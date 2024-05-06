# Change Log

All notable changes are documented here under headings (YYYY-MM-DD) in reverse
chronological order.

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



