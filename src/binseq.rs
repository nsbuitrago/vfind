use std::{
    collections::HashMap,
    sync::{Arc, Mutex},
};

use crate::processor::Processor;
use crate::translate::*;

use binseq::{
    BinseqReader, BinseqRecord, ParallelProcessor, ParallelReader, Result as BinseqResult,
};
use polars::df;
use pyo3::{PyResult, exceptions::PyValueError};
use pyo3_polars::PyDataFrame;

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
                    *self.local_variants.entry(variant).or_insert(0) += 1;
                }
            } else {
                // attempt to translate and add to variants
                if let Some(variant) = translate(&variant) {
                    // *self.variants.entry(variant).or_insert(0) += 1;
                    *self.local_variants.entry(variant).or_insert(0) += 1;
                }
            }
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> BinseqResult<()> {
        {
            let mut variants = self.variants.lock().unwrap();
            self.local_variants
                .iter()
                .for_each(|(seq, count)| *variants.entry(seq.into()).or_insert(0) += count);
            self.local_variants.clear();
        }

        Ok(())
    }
}

pub(crate) fn process_binseq(
    path: String,
    adapters: (String, String),
    match_score: i32,
    mismatch_score: i32,
    gap_open_penalty: i32,
    gap_extend_penalty: i32,
    accept_prefix_alignment: f64,
    accept_suffix_alignment: f64,
    n_threads: usize,
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
    )?;

    // show the progress
    println!("{}", show_progress);

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
