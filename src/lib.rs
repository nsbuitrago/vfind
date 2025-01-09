use std::{collections::HashMap, fs::File};

use memchr::memmem;
use parasail_rs::{Aligner, Matrix, Profile};
use polars::df;
use pyo3::{exceptions::PyValueError, prelude::*};
use pyo3_polars::PyDataFrame;
use seq_io::fastq::{Reader, Record};
use seq_io::parallel::parallel_fastq;

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

// Check that the accepted alignment thresholds are in the range (0, 1).
fn validate_align_thresholds(value: f32, name: &str) -> PyResult<()> {
    if value >= 1. || value <= 0. {
        Err(PyValueError::new_err(format!(
            "{} must be in range (0, 1)",
            name
        )))
    } else {
        Ok(())
    }
}

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
        .build()
}

/// Find where adapter matches in the read if it exists
fn find_adapter_match(
    seq: &[u8],
    adapter: &[u8],
    aligner: &Aligner,
    min_align_score: f32,
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
        if score as f32 > min_align_score {
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
    skip_translation=false,
    skip_alignment=false,
    n_threads=3,
    queue_len=2,
    show_progress=true,
))]
/// Find variable regions flanked by constant adapters in FASTQ dataset.
///
/// Args:
///     fq_path (str): Path to fastq file
///     adapters (tuple[str, str]): Tuple of adapters (prefix, suffix)
///     match_score (int): Match score for alignment (default = 3)
///     mismatch_score (int): Mismatch score for alignment (default = -2)
///     gap_open_penalty (int): Gap opening penalty for alignment (default = -2)
///     gap_extend_penalty (int): Gap extension penalty for alignment (default = -2)
///     accept_prefix_alignment (float): Accepted prefix alignment threshold (default = 0.75)
///     accept_prefix_alignment (float): Accepted suffix alignment threshold (default = 0.75)
///     skip_translation (bool): Skip translation to amino acid sequence (default = True)
///     skip_alignment (bool): Skip semi-global alignment step (default = False)
///     n_threads (int): Number of threads (default = 3)
///     queue_len (int): Queue length (default = 2)
///     show_progress (bool): Show progress bar (default = True)
pub fn find_variants(
    fq_path: String,
    adapters: (String, String),
    match_score: i32,
    mismatch_score: i32,
    gap_open_penalty: i32,
    gap_extend_penalty: i32,
    accept_prefix_alignment: f32,
    accept_suffix_alignment: f32,
    skip_translation: bool,
    skip_alignment: bool,
    n_threads: u32,
    queue_len: usize,
    show_progress: bool,
) -> PyResult<PyDataFrame> {
    // validate accept prefix/suffix accept thresholds
    validate_align_thresholds(accept_prefix_alignment, "accept_prefix_alignment")?;
    validate_align_thresholds(accept_suffix_alignment, "accept_suffix_alignment")?;
    // load data
    let file = File::open(fq_path)?;
    let mut buffer = buffer_redux::BufReader::new(file);
    buffer.read_into_buf().unwrap();
    let reader = Reader::new(flate2::bufread::GzDecoder::new(buffer));

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
    let min_prefix_score = accept_prefix_alignment * match_score as f32 * prefix.len() as f32;
    let min_suffix_score = accept_suffix_alignment * match_score as f32 * suffix.len() as f32;

    let mut variants: HashMap<String, u64> = HashMap::new();

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

    let (seq, count): (Vec<String>, Vec<u64>) = variants.into_iter().unzip();
    let df = df!(
        "sequence" => seq,
        "count" => count,
    )
    .map_err(|e| PyValueError::new_err(format!("Error creating DataFrame: {}", e)))?;

    Ok(PyDataFrame(df))
}

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
