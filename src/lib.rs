use flate2::read::GzDecoder;
use memchr::memmem;
use parasail_rs::{Aligner, Matrix, Profile};
use pyo3_polars::PyDataFrame;
use seq_io::fastq::{Reader, Record};
use seq_io::parallel::parallel_fastq;
use std::collections::HashMap;
use std::io::Read;
use std::{fs::File, path::Path};

use pyo3::{exceptions::PyIOError, prelude::*};

const N_THREADS: u32 = 4;
const QUEUE_LEN: usize = 2;

/// Collection of constant adapter sequences flanking the variable regions of interest.
#[derive(Debug)]
#[pyclass]
pub struct Adapters {
    prefix: Vec<u8>,
    suffix: Vec<u8>,
}

#[pymethods]
impl Adapters {
    /// Create a new collection of adapters
    #[new]
    pub fn new(prefix: String, suffix: String) -> PyResult<Self> {
        Ok(Adapters {
            prefix: prefix.into_bytes(),
            suffix: suffix.into_bytes(),
        })
    }
}

/// Alignment scoring parameters
#[pyclass]
pub struct AlignParams {
    #[pyo3(get, set)]
    match_score: i32,
    #[pyo3(get, set)]
    mismatch_score: i32,
    #[pyo3(get, set)]
    gap_open_penalty: i32,
    #[pyo3(get, set)]
    gap_extend_penalty: i32,
    #[pyo3(get, set)]
    accept_prefix_alignment: f64,
    #[pyo3(get, set)]
    accept_suffix_alignment: f64,
}

impl Default for AlignParams {
    fn default() -> Self {
        AlignParams {
            match_score: 3,
            mismatch_score: -2,
            gap_extend_penalty: 5,
            gap_open_penalty: 2,
            accept_prefix_alignment: 0.75,
            accept_suffix_alignment: 0.75,
        }
    }
}

pub fn translate(seq: &[u8]) -> String {
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
    peptide
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

#[pyfunction]
pub fn find_variants(
    fq_path: String,
    adapters: &Adapters,
    align_params: Option<&AlignParams>,
    n_threads: Option<u32>,
    queue_len: Option<usize>,
    skip_translation: Option<bool>,
) -> PyResult<PyDataFrame> {
    let fq_path = Path::new(&fq_path);
    if !fq_path.exists() {
        return Err(PyIOError::new_err("File {} does not exist."));
    }

    let fq_file = File::open(fq_path)?;
    let mut decoder = GzDecoder::new(&fq_file);
    let mut data = Vec::new();
    decoder.read_to_end(&mut data)?;

    let reader = Reader::new(data.as_slice());

    let n_threads = n_threads.unwrap_or(N_THREADS);
    let queue_len = queue_len.unwrap_or(QUEUE_LEN);

    let default_align_params = &AlignParams::default();
    let align_params = align_params.unwrap_or(default_align_params);
    let scoring_matrix = Matrix::create(
        b"ATCG",
        align_params.match_score,
        align_params.mismatch_score,
    )
    .expect("Error creating scoring matrix");

    let prefix_profile = Profile::new(&adapters.prefix, true, &scoring_matrix).unwrap();

    let prefix_aligner = Aligner::new()
        .profile(prefix_profile)
        .gap_open(align_params.gap_open_penalty)
        .gap_extend(align_params.gap_extend_penalty)
        .semi_global()
        .scan()
        .use_stats()
        .build();

    let suffix_profile = Profile::new(&adapters.suffix, true, &scoring_matrix).unwrap();

    let suffix_aligner = Aligner::new()
        .profile(suffix_profile)
        .gap_open(align_params.gap_open_penalty)
        .gap_extend(align_params.gap_extend_penalty)
        .semi_global()
        .scan()
        .use_stats()
        .build();

    let max_prefix_score = align_params.accept_prefix_alignment
        * align_params.match_score as f64
        * adapters.prefix.len() as f64;

    let max_suffix_score = align_params.accept_suffix_alignment
        * align_params.match_score as f64
        * adapters.suffix.len() as f64;

    let mut variants: HashMap<String, u32> = HashMap::new();

    parallel_fastq(
        reader,
        n_threads,
        queue_len,
        |record, variant| {
            let seq = record.seq();

            let start = {
                let exact_start = memmem::find(seq, &adapters.prefix);
                if let Some(exact_start) = exact_start {
                    Some(exact_start + adapters.prefix.len())
                } else {
                    let alignment = prefix_aligner.align(None, seq).unwrap();
                    let score = alignment.get_score();
                    if score as f64 > max_prefix_score {
                        Some(alignment.get_length().unwrap() as usize)
                    } else {
                        None
                    }
                }
            };

            let end = {
                let exact_end = memmem::find(seq, &adapters.suffix);
                if let Some(exact_end) = exact_end {
                    Some(exact_end)
                } else {
                    let alignment = suffix_aligner.align(None, seq).unwrap();
                    let score = alignment.get_score();
                    if score as f64 > max_suffix_score {
                        Some(seq.len() - alignment.get_length().unwrap() as usize)
                    } else {
                        None
                    }
                }
            };

            if start.is_some() && end.is_some() && start.unwrap() < end.unwrap() {
                *variant = Some(seq[start.unwrap()..end.unwrap()].to_vec());
            }
        },
        |_, variant| {
            if let Some(variant) = variant {
                if skip_translation.unwrap_or(false) {
                    let variant = String::from_utf8(variant.to_vec()).unwrap();
                    *variants.entry(variant).or_insert(0) += 1;
                } else {
                    // translate and add to variants
                    *variants.entry(translate(variant)).or_insert(0) += 1;
                }
            }

            None::<()>
        },
    )
    .unwrap();

    let df = polars::df!(
        "sequence" => variants.keys().cloned().collect::<Vec<String>>(),
        "count" => variants.values().cloned().collect::<Vec<u32>>(),
    )
    .unwrap();

    Ok(PyDataFrame(df))
}

/// vFind Python module
#[pymodule]
fn vfind(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Adapters>()?;
    m.add_class::<AlignParams>()?;
    m.add_function(wrap_pyfunction!(find_variants, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// sanity check
    #[test]
    fn test_adapter_construction() -> Result<(), PyErr> {
        let seq_a = String::from("ATG");
        let seq_b = String::from("TAA");

        let adapters = Adapters::new(seq_a, seq_b)?;
        assert_eq!(adapters.prefix, b"ATG");
        assert_eq!(adapters.suffix, b"TAA");

        Ok(())
    }

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
        let align_params = AlignParams::default();
        let accept_align_score = align_params.accept_prefix_alignment
            * align_params.match_score as f64
            * query.len() as f64;

        let scoring_matrix = Matrix::create(
            b"ACTG",
            align_params.match_score,
            align_params.mismatch_score,
        )
        .unwrap();

        let profile = Profile::new(query, true, &scoring_matrix).unwrap();
        let aligner = Aligner::new()
            .profile(profile)
            .gap_open(align_params.gap_open_penalty)
            .gap_extend(align_params.gap_extend_penalty)
            .semi_global()
            .scan()
            .use_stats()
            .build();

        let alignment = aligner.align(None, template).unwrap();
        let score = alignment.get_score();
        if pass {
            assert!(score as f64 > accept_align_score);
        } else {
            assert!(accept_align_score > score as f64);
        }
    }
}
