use std::{
    collections::HashMap,
    sync::{Arc, Mutex},
};

use crate::processor::Processor;
use binseq::{BinseqReader, ParallelReader};
use polars::df;
use pyo3::{PyResult, exceptions::PyValueError};
use pyo3_polars::PyDataFrame;

pub fn process_binseq(
    path: String,
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
    let reader = BinseqReader::new(path.as_ref()).unwrap();
    let variants: Arc<Mutex<HashMap<String, u64>>> = Arc::new(Mutex::new(HashMap::new()));
    let processor = Processor::new(
        &adapters,
        match_score,
        mismatch_score,
        gap_open_penalty,
        gap_extend_penalty,
        accept_prefix_alignment,
        accept_suffix_alignment,
        skip_translation,
        variants.clone(),
    )?;

    reader.process_parallel(processor, n_threads).unwrap();
    let mut variant_map = variants.lock().unwrap();
    let owned_map = std::mem::take(&mut *variant_map);
    let (seq, count): (Vec<String>, Vec<u64>) = owned_map.into_iter().unzip();
    let df = df!(
        "sequence" => seq,
        "count" => count,
    )
    .map_err(|e| PyValueError::new_err(format!("Error creating DataFrame: {}", e)))?;
    Ok(PyDataFrame(df))
}
