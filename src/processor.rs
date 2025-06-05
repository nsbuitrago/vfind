use std::{
    collections::HashMap,
    sync::{Arc, Mutex},
};

use memchr::memmem;
use parasail_rs::{Aligner, Matrix, Profile};
use pyo3::{PyResult, exceptions::PyValueError};

// Checks if the accepted alignment threshold is valid. Thresholds must
// be a float in the range (0, 1]. If the threshold is set to 1, alignments
// are skipped.
pub(crate) fn align_threshold_preflight(align_threshold: f64) -> PyResult<bool> {
    if align_threshold > 0. && align_threshold < 1. {
        return Ok(false);
    } else if align_threshold == 1. {
        return Ok(true);
    }

    Err(PyValueError::new_err(
        "Accept alignment threshold must be between 0 and 1.",
    ))
}

/// Construct and return aligner with an adapter profile.
pub(crate) fn get_aligner(
    adapter: &[u8],
    scoring_matrix: &Matrix,
    gap_open_penalty: i32,
    gap_extend_penalty: i32,
    skip_alignment: bool,
) -> Option<Aligner> {
    if skip_alignment {
        return None;
    }

    let adapter_profile = Profile::new(adapter, true, scoring_matrix)
        .map_err(|e| PyValueError::new_err(format!("Error creating profile for adapter: {}", e)))
        .unwrap();

    let aligner = Aligner::new()
        .profile(adapter_profile)
        .gap_open(gap_open_penalty)
        .gap_extend(gap_extend_penalty)
        .semi_global()
        .scan()
        .use_stats()
        .build();

    Some(aligner)
}

/// Find where adapter matches in the read if it exists
pub(crate) fn find_adapter_match(
    seq: &[u8],
    adapter: &[u8],
    aligner: Option<&Aligner>,
    min_align_score: f64,
    is_prefix: bool,
) -> Option<usize> {
    let exact_match = memmem::find(seq, adapter);
    if let Some(exact_match) = exact_match {
        match is_prefix {
            true => Some(exact_match + adapter.len()),
            false => Some(exact_match),
        }
    } else {
        let alignment = aligner?.align(None, seq).unwrap();
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

#[derive(Clone)]
pub(crate) struct Processor {
    pub(crate) prefix: Vec<u8>,
    pub(crate) suffix: Vec<u8>,
    pub(crate) prefix_aligner: Option<Aligner>,
    pub(crate) suffix_aligner: Option<Aligner>,
    pub(crate) min_prefix_score: f64,
    pub(crate) min_suffix_score: f64,
    pub(crate) skip_translation: bool,
    pub(crate) local_variants: HashMap<String, u64>,
    pub(crate) variants: Arc<Mutex<HashMap<String, u64>>>, // FIXME: does dashmap improve performance?
}

impl Processor {
    pub fn new(
        adapters: &(String, String),
        match_score: i32,
        mismatch_score: i32,
        gap_open_penalty: i32,
        gap_extend_penalty: i32,
        accept_prefix_alignment: f64,
        accept_suffix_alignment: f64,
        skip_translation: bool,
    ) -> PyResult<Self> {
        let scoring_matrix = Matrix::create(b"ACGT", match_score, mismatch_score)
            .expect("Error creating scoring matrix");

        let skip_prefix_alignment = align_threshold_preflight(accept_prefix_alignment)?;
        let prefix = adapters.0.as_bytes();
        let prefix_aligner = get_aligner(
            prefix,
            &scoring_matrix,
            gap_open_penalty,
            gap_extend_penalty,
            skip_prefix_alignment,
        );

        let skip_suffix_alignment = align_threshold_preflight(accept_suffix_alignment)?;
        let suffix = adapters.1.as_bytes();
        let suffix_aligner = get_aligner(
            suffix,
            &scoring_matrix,
            gap_open_penalty,
            gap_extend_penalty,
            skip_suffix_alignment,
        );

        // min score for accepted prefix and suffix alignment
        let min_prefix_score = accept_prefix_alignment * match_score as f64 * prefix.len() as f64;
        let min_suffix_score = accept_suffix_alignment * match_score as f64 * suffix.len() as f64;

        let local_variants: HashMap<String, u64> = HashMap::new();
        let variants: Arc<Mutex<HashMap<String, u64>>> = Arc::new(Mutex::new(HashMap::new()));

        Ok(Self {
            prefix: prefix.to_vec(),
            suffix: suffix.to_vec(),
            prefix_aligner,
            suffix_aligner,
            min_prefix_score,
            min_suffix_score,
            skip_translation,
            local_variants,
            variants,
        })
    }

    pub(crate) fn get_matches(&self, seq: &Vec<u8>) -> (Option<usize>, Option<usize>) {
        let start = find_adapter_match(
            seq,
            &self.prefix,
            self.prefix_aligner.as_ref(),
            self.min_prefix_score,
            true,
        );

        let end = find_adapter_match(
            seq,
            &self.suffix,
            self.suffix_aligner.as_ref(),
            self.min_suffix_score,
            false,
        );

        return (start, end);
    }

    // fn get_dataframe(self) -> PyResult<PyDataFrame> {
    //     let (seq, count): (Vec<String>, Vec<u64>) = self.variants.into_iter().unzip();
    //     let df = df!(
    //         "sequence" => seq,
    //         "count" => count,
    //     )
    //     .map_err(|e| PyValueError::new_err(format!("Error creating DataFrame: {}", e)))?;

    //     Ok(PyDataFrame(df))
    // }
}
