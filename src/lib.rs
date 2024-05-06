use flate2::read::GzDecoder;
use indicatif::{ProgressBar, ProgressState, ProgressStyle};
use memchr::memmem;
use parasail_rs::{Aligner, Matrix, Profile};
use polars::prelude::*;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;
use seq_io::fastq::{Reader, Record};
use seq_io::parallel::parallel_fastq;
use std::collections::HashMap;
use std::fmt::Write;
use std::io::Read;
use std::{fs::File, path::Path};

/// Attempt to translate a DNA sequence to an amino acid sequence. If the sequence
/// length is not divisible by 3, the None variant is returned.
pub fn translate(seq: &[u8]) -> Option<String> {
    if seq.len() % 3 != 0 {
        return None;
    }

    let mut peptide = String::with_capacity(seq.len() / 3);

    'outer: for triplet in seq.chunks_exact(3) {
        for c in triplet {
            if !c.is_ascii() {
                peptide.push('X');
                continue 'outer;
            }
        }

        let c1 = ASCII_TO_INDEX[triplet[0] as usize];
        let c2 = ASCII_TO_INDEX[triplet[1] as usize];
        let c3 = ASCII_TO_INDEX[triplet[2] as usize];

        let amino_acid = if c1 == 4 || c2 == 4 || c3 == 4 {
            'X'
        } else {
            AA_TABLE_CANONICAL[c1][c2][c3]
        };

        peptide.push(amino_acid);
    }
    Some(peptide)
}

/// https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
/// U is equivalent to T here
///
/// The 1st index picks the 4x4 block
/// The 2nd index picks the row
/// the 3rd index picks the column
static AA_TABLE_CANONICAL: [[[char; 4]; 4]; 4] = [
    [
        ['K', 'N', 'K', 'N'], // AAA, AAC, AAG, AAU/AAT
        ['T', 'T', 'T', 'T'], // ACA, ACC, ACG, ACU/ACT
        ['R', 'S', 'R', 'S'], // AGA, AGC, AGG, AGU/AGT
        ['I', 'I', 'M', 'I'], // AUA/ATA, AUC/ATC, AUG/ATG, AUU/ATT
    ],
    [
        ['Q', 'H', 'Q', 'H'], // CAA, CAC, CAG, CAU/CAT
        ['P', 'P', 'P', 'P'], // CCA, CCC, CCG, CCU/CCT
        ['R', 'R', 'R', 'R'], // CGA, CGC, CGG, CGU/CGT
        ['L', 'L', 'L', 'L'], // CUA/CTA, CUC/CTC, CUG/CTG, CUU/CTT
    ],
    [
        ['E', 'D', 'E', 'D'], // GAA, GAC, GAG, GAU/GAT
        ['A', 'A', 'A', 'A'], // GCA, GCC, GCG, GCU/GCT
        ['G', 'G', 'G', 'G'], // GGA, GGC, GGG, GGU/GGT
        ['V', 'V', 'V', 'V'], // GUA/GTA, GUC/GTC, GUG/GTG, GUU/GTT
    ],
    [
        ['*', 'Y', '*', 'Y'], // UAA/TAA, UAC/TAC, UAG/TAG, UAU/TAT
        ['S', 'S', 'S', 'S'], // UCA/TCA, UCC/TCC, UCG/TCG, UCU/TCT
        ['*', 'C', 'W', 'C'], // UGA/TGA, UGC/TGC, UGG/TGG, UGU/TGT
        ['L', 'F', 'L', 'F'], // UUA/TTA, UUC/TTC, UUG/TTG, UUU/TTT
    ],
];

/// Maps an ASCII character to array index
///
/// A = 65, a = 97  => 0
/// C = 67, c = 99  => 1
/// G = 71, g = 103 => 2
/// T = 84, t = 116 => 3
/// U = 85, u = 117 => 3
static ASCII_TO_INDEX: [usize; 128] = [
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 0-15
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 16-31
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 32-47
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 48-63
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 64-79    (65 = A, 67 = C, 71 = G)
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 80-95    (84 = T, 85 = U)
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 96-111   (97 = a, 99 = c, 103 = g)
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 112-127  (116 = t, 117 = u)
];

/// Construct and return aligner with an adapter profile.
fn get_aligner(
    adapter: &[u8],
    scoring_matrix: &Matrix,
    gap_open_penalty: i32,
    gap_extend_penalty: i32,
) -> Aligner {
    let adapter_profile = Profile::new(adapter, true, scoring_matrix)
        .map_err(|e| PyValueError::new_err(format!("Error creating profile for adapter: {}", e)))
        .unwrap();

    Aligner::new()
        .profile(adapter_profile)
        .gap_open(gap_open_penalty)
        .gap_extend(gap_extend_penalty)
        .semi_global()
        .scan()
        .use_stats()
        .build()
}

/// Find where adapter matches in the read if it exists
fn find_adapter_match(
    seq: &[u8],
    adapter: &[u8],
    aligner: &Aligner,
    min_align_score: f64,
    is_prefix: bool,
    skip_alignment: bool,
) -> Option<usize> {
    let exact_match = memmem::find(seq, adapter);
    if let Some(exact_match) = exact_match {
        match is_prefix {
            true => Some(exact_match + adapter.len()),
            false => Some(exact_match),
        }
    } else {
        if skip_alignment {
            return None;
        }

        let alignment = aligner.align(None, seq).unwrap();
        let score = alignment.get_score();
        if score as f64 > min_align_score {
            match is_prefix {
                true => Some(alignment.get_length().unwrap() as usize),
                false => Some(seq.len() - alignment.get_length().unwrap() as usize),
            }
        } else {
            None
        }
    }
}

#[pyfunction]
#[pyo3(signature = (
    fq_path,
    adapters,
    match_score=3,
    mismatch_score=-2,
    gap_open_penalty=5,
    gap_extend_penalty=2,
    accept_prefix_alignment=0.75,
    accept_suffix_alignment=0.75,
    n_threads=3,
    queue_len=2,
    skip_translation=false,
    skip_alignment=false,
    show_progress=true,
))]
/// Find variable regions flanked by adapters in a FASTQ dataset.
///
/// Args:
///     fq_path (str): Path to FASTQ file
///     adapters (tuple(str, str)): Tuple of prefix and suffix adapters
///     match_score (int): Match score for alignment (default = 3)
///     mismatch_score (int): Mismatch score for alignment (default = -2)
///     gap_open_penalty (int): Gap open penalty for alignment (default = 5)
///     gap_extend_penalty (int): Gap extend penalty for alignment (default = 2)
///     accept_prefix_alignment (float): Threshold for accepting prefix alignment (default = 0.75)
///     accept_suffix_alignment (float): Threshold for accepting suffix alignment (default = 0.75)
///     n_threads (int): Number of threads (default = 3)
///     queue_len (int): Queue length (default = 2)
///     skip_translation (bool): Skip translation (default = False)
///     skip_alignment (bool): Skip alignments (default = False)
///     show_progress (bool): Show progress bar (default = True)
///
/// Returns:
///     polars.DataFrame: dataframe with columns 'sequence' and 'count'
pub fn find_variants(
    fq_path: String,
    adapters: (String, String),
    match_score: i32,
    mismatch_score: i32,
    gap_open_penalty: i32,
    gap_extend_penalty: i32,
    accept_prefix_alignment: f64,
    accept_suffix_alignment: f64,
    n_threads: u32,
    queue_len: usize,
    skip_translation: bool,
    skip_alignment: bool,
    show_progress: bool,
) -> PyResult<PyDataFrame> {
    if accept_prefix_alignment > 1.0 || accept_prefix_alignment <= 0.0 {
        return Err(PyValueError::new_err(
            "accept_prefix_alignment must be between 0 and 1",
        ));
    }

    if accept_suffix_alignment > 1.0 || accept_suffix_alignment <= 0.0 {
        return Err(PyValueError::new_err(
            "accept_suffix_alignment must be between 0 and 1",
        ));
    }

    assert!(accept_suffix_alignment > 0.0 && accept_suffix_alignment <= 1.0);
    if accept_prefix_alignment == 1.0 && accept_suffix_alignment == 1.0 {
        return Err(PyValueError::new_err(
            "Both accept_prefix_alignment and accept_suffix_alignment are set to 1.0. This will result in only accepting exact matches of adapters. Set skip_alignment=True if this is the intended behavior.",
        ));
    }

    let fq_path = Path::new(&fq_path);
    let fq_file = File::open(fq_path)?;
    let mut decoder = GzDecoder::new(&fq_file);
    let mut data = Vec::new();
    decoder.read_to_end(&mut data)?;
    let reader = Reader::new(data.as_slice());

    let scoring_matrix = Matrix::create(b"ATCG", match_score, mismatch_score)
        .expect("Error creating scoring matrix");

    let prefix = adapters.0.as_bytes();
    let prefix_aligner = get_aligner(
        prefix,
        &scoring_matrix,
        gap_open_penalty,
        gap_extend_penalty,
    );

    let suffix = adapters.1.as_bytes();
    let suffix_aligner = get_aligner(
        suffix,
        &scoring_matrix,
        gap_open_penalty,
        gap_extend_penalty,
    );

    // min score for accepted prefix and suffix alignment
    let min_prefix_score = accept_prefix_alignment * match_score as f64 * prefix.len() as f64;
    let min_suffix_score = accept_suffix_alignment * match_score as f64 * suffix.len() as f64;

    let mut variants: HashMap<String, u64> = HashMap::new();

    let pb = if show_progress {
        let pb = ProgressBar::new(data.len() as u64);
        pb.set_style(ProgressStyle::with_template("{spinner:.green} {msg} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {bytes}/{total_bytes} ({eta})")
            .unwrap()
            .with_key("eta", |state: &ProgressState, w: &mut dyn Write| write!(w, "{:.1}s", state.eta().as_secs_f64()).unwrap())
            .progress_chars("#>-"));
        pb
    } else {
        ProgressBar::hidden()
    };
    pb.set_message(format!("Processing {}", fq_path.display()));

    parallel_fastq(
        reader,
        n_threads,
        queue_len,
        |record, variant| {
            // find variable region in the read
            let seq = record.seq();
            let start = find_adapter_match(
                seq,
                prefix,
                &prefix_aligner,
                min_prefix_score,
                true,
                skip_alignment,
            );
            let end = find_adapter_match(
                seq,
                suffix,
                &suffix_aligner,
                min_suffix_score,
                false,
                skip_alignment,
            );

            if start.is_some() && end.is_some() && start.unwrap() < end.unwrap() {
                *variant = Some(seq[start.unwrap()..end.unwrap()].to_vec());
            }

            // determine number of bytes read and update progress
            let (id, desc) = record.id_desc_bytes();
            let n_bytes =
                record.seq().len() + record.qual().len() + id.len() + desc.unwrap_or(&[]).len();
            pb.inc(n_bytes as u64);
        },
        |_, variant| {
            if let Some(variant) = variant {
                if skip_translation {
                    if let Ok(variant) = String::from_utf8(variant.to_vec()) {
                        *variants.entry(variant).or_insert(0) += 1;
                    }
                } else {
                    // Attempt to translate and add to variants
                    if let Some(variant) = translate(variant) {
                        *variants.entry(variant).or_insert(0) += 1;
                    }
                }
            }
            None::<()>
        },
    )
    .unwrap();

    pb.finish_and_clear();

    let (seq, count): (Vec<String>, Vec<u64>) = variants.into_iter().unzip();
    let df = df!(
        "sequence" => seq,
        "count" => count,
    )
    .map_err(|e| PyValueError::new_err(format!("Error creating DataFrame: {}", e)))?;

    Ok(PyDataFrame(df))
}

/// vFind Python module
#[pymodule]
fn vfind(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(find_variants, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const MATCH_SCORE: i32 = 3;
    const MISMATCH_SCORE: i32 = -2;
    const GAP_OPEN_PENALTY: i32 = 5;
    const GAP_EXTEND_PENALTY: i32 = 2;
    const ACCEPT_ALIGNMENT: f64 = 0.75;

    #[test]
    fn test_good_prefix_alignment() {
        let prefix = b"GGGCCCAGCCGGCCGGAT";
        let seq = b"GGGCCCAGCCGGCGGGATATGGCGGGCATCTGTGCACTTCCGGAGGCGGAGGTTCAG";
        let good_alignment = true;
        _test_alignment_helper(prefix, seq, good_alignment);
    }

    #[test]
    fn test_bad_prefix_alignment() {
        let prefix = b"GGGCCCAGCCGGCCGGAT";
        let seq = b"GGGTCTAACCGGCGGGATATGGCGGGCATCTGTGCACTTCCGGAGGCGGAGGTTCAG";
        let good_alignment = false;
        _test_alignment_helper(prefix, seq, good_alignment)
    }

    #[test]
    fn test_good_suffix_alignment() {
        let suffix = b"CCGGAGGCGGAGGTTCAG";
        let seq = b"GGGCCCAGCCGGCCGGATATGGCGGGCATCTGTGCACTTCCGGAGGCGGAGGTTGAG";
        let good_alignment = true;
        _test_alignment_helper(suffix, seq, good_alignment);
    }

    #[test]
    fn test_bad_suffix_alignment() {
        let suffix = b"CCGGAGGCGGAGGTTCAG";
        let seq = b"GGGCCCAGCCGGCCGGATATGGCGGGCATCTGTGCACTTCCGGAGGCCGTGCTCCAC";
        let good_alignment = false;
        _test_alignment_helper(suffix, seq, good_alignment)
    }

    /// helper fn to test alignments
    fn _test_alignment_helper(query: &[u8], template: &[u8], pass: bool) {
        let min_align_score = ACCEPT_ALIGNMENT * MATCH_SCORE as f64 * query.len() as f64;

        let scoring_matrix = Matrix::create(b"ACTG", MATCH_SCORE, MISMATCH_SCORE).unwrap();

        let profile = Profile::new(query, true, &scoring_matrix).unwrap();
        let aligner = Aligner::new()
            .profile(profile)
            .gap_open(GAP_OPEN_PENALTY)
            .gap_extend(GAP_EXTEND_PENALTY)
            .semi_global()
            .scan()
            .use_stats()
            .build();

        let alignment = aligner.align(None, template).unwrap();
        let score = alignment.get_score();
        if pass {
            assert!(score as f64 > min_align_score);
        } else {
            assert!(min_align_score > score as f64);
        }
    }
}
