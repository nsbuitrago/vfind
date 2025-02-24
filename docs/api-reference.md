# API Reference

```python
def vfind.find_variants(
    fq_path: str,
    adapters: tuple(str, str),
    match_score: int = 3,
    mismatch_score: int = -2,
    gap_open_penalty: int = 5,
    gap_extend_penalty: int = 2,
    accept_prefix_alignment: float = 0.75,
    accept_suffix_alignment: float = 0.75,
    skip_translation: bool = False,
    skip_alignment: bool = False,
    show_progress: bool = True,
    n_threads: int = 3,
    queue_len: int = 2,
)
```

**Find variable regions flanked by constant adapters.**

## Arguments

**fq_path: str \[Required\]:**  Path to fastq file.

**adapters: tuple(str, str) \[Required\]:** Constant adapter sequences flanking
variable region.

**match_score: int \[Optional, Default = 3\]:** Match score used for alignment.

**mismatch_score: int \[Optional, Default = -2\]:** Mismatch score used for
alignment.

**gap_open_penalty: int \[Optional, Default = 5\]:** Gap open penalty used for
alignment. Note that this is given as a positive integer.

**gap_extend_penalty: int \[Optional, Default = 2\]:** Gap extension penalty
used for alignment. Note that this is given as a positive integer.

**accept_prefix_alignment: float (0, 1] \[Optional, Default = 0.75\]:** Threshold for
accepting alignments between read and 5' adapter (prefix) sequences.

**accept_suffix_alignment: float (0, 1] \[Optional, Default = 0.75\]:** Threshold for
accepting alignments between read and 3' adapter (suffix) sequences.

Note that if you use an accept threshold of 1, you should set `skip_alignment` to `True`.

**skip_translation: bool \[Optional, Default = False\]:** Whether to skip
in-frame translation and return DNA sequences.

**skip_alignment: bool \[Optional, Default = False\]:** Whether to skip
semi-global alignment and only search for exact adapter matches.

**show_progress: bool \[Optional, Default = True\]:** Whether to show
visual progress bar.*

**n_threads: int \[Optional, Default = 3\]:** Number of threads to use for
processing fastq.

**queue_len: int \[Optional, Default = 2\]:** Queue length used for processing
fastq.

## Returns

**polars.DataFrame:** Dataframe with `sequence` and `count` columns
corresponding to recovered sequences mapped to raw read counts.

