use std::{
    collections::HashMap,
    sync::{Arc, Mutex},
};

use binseq::{BinseqRecord, ParallelProcessor, Result as BinseqResult};
use memchr::memmem;
use parasail_rs::{Aligner, Matrix, Profile};
use pyo3::{PyResult, exceptions::PyValueError};

use crate::translate::translate;

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
pub struct Processor {
    prefix: Vec<u8>,
    suffix: Vec<u8>,
    prefix_aligner: Option<Aligner>,
    suffix_aligner: Option<Aligner>,
    min_prefix_score: f64,
    min_suffix_score: f64,
    skip_translation: bool,
    local_variants: Vec<String>,
    variants: Arc<Mutex<HashMap<String, u64>>>,
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
        variants: Arc<Mutex<HashMap<String, u64>>>,
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

        let mut local_variants: Vec<String> = Vec::new();
        // let mut variants: Arc<Mutex<HashMap<String, u64>>> = Arc::new(Mutex::new(HashMap::new()));

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

    fn get_matches(&self, seq: &Vec<u8>) -> (Option<usize>, Option<usize>) {
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

impl ParallelProcessor for Processor {
    fn process_record<B: BinseqRecord>(&mut self, record: B) -> BinseqResult<()> {
        // per record logic
        // get the seq
        let mut seq: Vec<u8> = Vec::new();
        record.decode_s(&mut seq)?;
        // let (start, end) = self.get_matches(seq);

        let mut variant: Option<Vec<u8>> = None;
        if let (Some(start), Some(end)) = self.get_matches(&seq) {
            if start < end {
                variant = Some(seq[start..end].to_vec());
            }
        }

        if let Some(variant) = variant {
            if self.skip_translation {
                if let Ok(variant) = String::from_utf8(variant) {
                    // *self.variants.entry(variant).or_insert(0) += 1;
                    self.local_variants.push(variant);
                }
            } else {
                // attempt to translate and add to variants
                if let Some(variant) = translate(&variant) {
                    // *self.variants.entry(variant).or_insert(0) += 1;
                    self.local_variants.push(variant);
                }
            }
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> BinseqResult<()> {
        {
            let mut variants = self.variants.lock().unwrap();
            self.local_variants.iter().for_each(|seq| {
                *variants.entry(seq.into()).or_insert(0) += 1;
            });
            self.local_variants.clear();
        }

        Ok(())
    }
}
