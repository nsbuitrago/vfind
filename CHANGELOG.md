# Change Log

All notable changes are documented here under headings (YYYY-MM-DD) in reverse
chronological order.

## 0.3.1

### Update

bump deps:
    - flate2 1.1.0 -> 1.1.1
    - log 0.4.26 -> 0.4.27
    - parasail-rs 0.7.7 -> 0.7.8
    - polars 0.46.0 -> 0.48.1
    - pyo3 0.23.4 -> 0.24.2
    - pyo3-polars 0.20.0 -> 0.21.0
    - seq_io 0.3.3 -> 0.3.4

## 0.3.0 - 2025-05-25

### Changes

- Remove `skip_alignment` keyword argument in favor of setting `accept_*_alignment`
to 1.

## 0.2.0 - 2025-05-24

### Bug Fixes

- fix issue where only first stream in multi-gz fastq files were being processed
(#3)

### Improvements

- Updated docs

## 0.1.1 - 2024-05-01

### Features

- Added `skip_alignment` to optionally skip semi-global alignment and only
use exact adapter matches to find variants.
- `find_variants` now has visual feedback for progress since jobs are
    usually long. To turn this off, set `show_progress` to `False`.

### Bug Fixes

- DNA translation assumed coding sequences and would fail for sequences not
    divisible by 3. Non-coding sequences are skipped by default.
    To include non-coding sequences, set `skip_translation` to `False` to get
    all DNA sequences that are recovered before translation.

### Misc.

- Support x86_64 linux


## 0.1.0 - 2024-04-15

Initial release.

### Features

- Find variants flanked by constant adapter sequences in gzipped FASTQ files.
