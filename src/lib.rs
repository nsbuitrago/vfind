use std::{collections::HashMap, fs::File};

use memchr::memmem;
use parasail_rs::{Aligner, Matrix, Profile};
use polars::df;
use pyo3::{exceptions::PyValueError, prelude::*};
use pyo3_polars::PyDataFrame;
use seq_io::fastq::{Reader, Record};
use seq_io::parallel::parallel_fastq;

/// Attempt to translate a DNA sequence to an amino acid sequence. If the
/// sequence length is not divisible by 3, the None variant is returned.
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

/// Check that the accepted alignment thresholds are in the range (0, 1).
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

/// Where the adapter is located relative to the variable region.
enum Location {
    Prefix,
    Suffix,
}
/// Find any exact matches of the adapter on the sequence.
fn exact_adapter_match(seq: &[u8], adapter: &[u8], location: &Location) -> Option<usize> {
    if let Some(exact_match) = memmem::find(seq, adapter) {
        match location {
            Location::Prefix => Some(exact_match + adapter.len()),
            Location::Suffix => Some(exact_match),
        }
    } else {
        None
    }
}

/// Align adapter to sequence and check if produced score is greater than the
/// minimum score.
fn align_adapter(
    seq: &[u8],
    aligner: &Aligner,
    min_align_score: f32,
    orientation: &Location,
) -> Option<usize> {
    let alignment = aligner.align(None, seq).unwrap();
    let score = alignment.get_score();
    if score as f32 > min_align_score {
        match orientation {
            Location::Prefix => Some(alignment.get_length().unwrap() as usize),
            Location::Suffix => Some(seq.len() - alignment.get_length().unwrap() as usize),
        }
    } else {
        None
    }
}

/// Find where adapter matches in the read if it exists by exact match or
/// optional semi-global alignment.
fn find_adapter_match(
    seq: &[u8],
    adapter: &[u8],
    aligner: &Aligner,
    min_align_score: f32,
    orientation: &Location,
) -> Option<usize> {
    exact_adapter_match(seq, adapter, orientation)
        .or_else(|| align_adapter(seq, aligner, min_align_score, orientation))
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
                &Location::Prefix,
            );
            let end = find_adapter_match(
                seq,
                suffix,
                &suffix_aligner,
                min_suffix_score,
                &Location::Suffix,
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

    #[test]
    fn translate_dna_seq() {
        let dna_seq = b"TTCTTAATTATGGTCTCTCCTACTGCCTACCATCAAAATAAAGATGAATGCTGGCGTGGGTAA";
        let amino_acid_seq = "FLIMVSPTAYHQNKDECWRG*";

        let aa_from_dna = translate(dna_seq);

        assert!(aa_from_dna.is_some());
        assert_eq!(aa_from_dna.unwrap(), amino_acid_seq);
    }

    #[test]
    fn exact_adapter_matches() {
        let seq = b"ACGTNNNNTGCA";
        let fwd_adapter = b"ACGT";
        let rev_adapter = b"TGCA";

        let variant_start_idx = exact_adapter_match(seq, fwd_adapter, &Location::Prefix);
        let variant_end_idx = exact_adapter_match(seq, rev_adapter, &Location::Suffix);

        assert!(variant_start_idx.is_some());
        assert!(variant_end_idx.is_some());

        assert_eq!(variant_start_idx.unwrap(), 4);
        assert_eq!(variant_end_idx.unwrap(), 8);

        let variant = &seq[variant_start_idx.unwrap()..variant_end_idx.unwrap()];
        assert_eq!(variant, b"NNNN");
    }

    #[test]
    fn aligned_adapter_match() {
        let seq = b"ACCTNNNNTGGA"; // fwd and rev adapters have single nucleotide substitution
        let fwd_adapter = b"ACGT";
        let rev_adapter = b"TGCA";

        let exact_prefix_match = exact_adapter_match(seq, fwd_adapter, &Location::Prefix);
        let exact_suffix_match = exact_adapter_match(seq, rev_adapter, &Location::Suffix);

        assert!(exact_prefix_match.is_none());
        assert!(exact_suffix_match.is_none());

        // setup aligners
        let align_threshold = 0.5;
        let match_score = 3;
        let mismatch_score = -2;
        let gap_open_penalty = 5;
        let gap_extend_penalty = 2;
        let scoring_matrix = Matrix::create(b"ACGT", match_score, mismatch_score).unwrap();

        let fwd_aligner = get_aligner(
            fwd_adapter,
            &scoring_matrix,
            gap_open_penalty,
            gap_extend_penalty,
        );
        let rev_aligner = get_aligner(
            rev_adapter,
            &scoring_matrix,
            gap_open_penalty,
            gap_extend_penalty,
        );

        let min_prefix_align_score =
            fwd_adapter.len() as f32 * match_score as f32 * align_threshold;
        let min_suffix_align_score = min_prefix_align_score;

        let aligned_prefix_match =
            align_adapter(seq, &fwd_aligner, min_prefix_align_score, &Location::Prefix);

        let aligned_suffix_match =
            align_adapter(seq, &rev_aligner, min_suffix_align_score, &Location::Suffix);

        assert!(aligned_prefix_match.is_some());
        assert!(aligned_suffix_match.is_some());
        assert_eq!(aligned_prefix_match.unwrap(), 4);
        assert_eq!(aligned_suffix_match.unwrap(), 8);
    }
}
