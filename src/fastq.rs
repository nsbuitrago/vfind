use flate2::read::MultiGzDecoder;
use indicatif::ProgressBar;
use polars::prelude::*;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;
use seq_io::fastq::{Reader, Record};
use seq_io::parallel::parallel_fastq;
use std::collections::HashMap;
use std::fs::File;

use crate::processor::Processor;
use crate::translate::translate;

pub(crate) fn process_fq(
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
    show_progress: bool,
) -> PyResult<PyDataFrame> {
    let gzdecoder = File::open(fq_path).map(MultiGzDecoder::new)?;
    let reader = Reader::new(gzdecoder);

    let processor = Processor::new(
        &adapters,
        match_score,
        mismatch_score,
        gap_open_penalty,
        gap_extend_penalty,
        accept_prefix_alignment,
        accept_suffix_alignment,
        skip_translation,
    )?;

    let mut variants: HashMap<String, u64> = HashMap::new();

    let pb = if show_progress {
        ProgressBar::new_spinner()
    } else {
        ProgressBar::hidden()
    };

    parallel_fastq(
        reader,
        n_threads,
        queue_len,
        |record, variant| {
            // find variable region in the read
            let seq = record.seq();
            // let mut variant: Option<Vec<u8>> = None;
            if let (Some(start), Some(end)) = processor.get_matches(&seq.to_vec()) {
                if start < end {
                    *variant = Some(seq[start..end].to_vec());
                }
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

    pb.finish_and_clear();

    let (seq, count): (Vec<String>, Vec<u64>) = variants.into_iter().unzip();
    let df = df!(
        "sequence" => seq,
        "count" => count,
    )
    .map_err(|e| PyValueError::new_err(format!("Error creating DataFrame: {}", e)))?;

    Ok(PyDataFrame(df))
}
