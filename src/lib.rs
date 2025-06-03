use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;

mod binseq;
mod fastq;
mod processor;
mod translate;

use crate::binseq::process_binseq;
use crate::fastq::process_fq;

#[pyfunction]
#[pyo3(signature = (
    filepath,
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
    show_progress=true,
))]
/// Find variable regions flanked by adapters in a FASTQ dataset.
///
/// Arguments
/// ----------
///
///     filepath : str
///         Path to FASTQ or Binseq file.
///     adapters : tuple(str, str)
///         Tuple of prefix and suffix adapters
///     match_score : int
///         Match score for alignment (Optional, default = 3)
///     mismatch_score : int
///         Mismatch score for alignment (Optional, default = -2)
///     gap_open_penalty : int
///         Gap open penalty for alignment (Optional, default = 5)
///     gap_extend_penalty : int
///         Gap extend penalty for alignment (Optional, default = 2)
///     accept_prefix_alignment : float
///         Threshold for accepting prefix alignment (Optional, default = 0.75)
///     accept_suffix_alignment : float
///         Threshold for accepting suffix alignment (Optional, default = 0.75)
///     n_threads : int
///         Number of threads (Optional, default = 3)
///     queue_len : int
///         Queue length (Optional, default = 2)
///     skip_translation : bool
///         Skip translation to amino-acid sequence (Optional, default = False)
///     skip_alignment : bool
///         Skip semi-global alignments (Optional, default = False)
///     show_progress : bool
///         Show progress bar (Optional, default = True)
///
/// Returns
/// -------
///
///     polars.DataFrame: dataframe with 'sequence' and 'count' columns
pub fn find_variants(
    filepath: String,
    adapters: (String, String),
    match_score: i32,
    mismatch_score: i32,
    gap_open_penalty: i32,
    gap_extend_penalty: i32,
    accept_prefix_alignment: f64,
    accept_suffix_alignment: f64,
    n_threads: usize,
    queue_len: usize,
    skip_translation: bool,
    show_progress: bool,
) -> PyResult<PyDataFrame> {
    // check if it is binseq or fastq
    if filepath.contains("fastq") || filepath.contains("fq") {
        process_fq(
            filepath,
            adapters,
            match_score,
            mismatch_score,
            gap_open_penalty,
            gap_extend_penalty,
            accept_prefix_alignment,
            accept_suffix_alignment,
            n_threads as u32,
            queue_len,
            skip_translation,
            show_progress,
        )
    } else if filepath.contains("bq") || filepath.contains("vbq") {
        process_binseq(
            filepath,
            adapters,
            match_score,
            mismatch_score,
            gap_open_penalty,
            gap_extend_penalty,
            accept_prefix_alignment,
            accept_suffix_alignment,
            n_threads,
            skip_translation,
            show_progress,
        )
    } else {
        panic!("Filetype could not be found");
    }
}

/// vFind Python module
#[pymodule]
fn vfind(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(find_variants, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use parasail_rs::{Aligner, Matrix, Profile};

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
