# vFind

[![PyPI - Version](https://img.shields.io/pypi/v/vfind)](https://pypi.org/project/vfind/)

*A simple variant finder for NGS data.*

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Examples](#examples)
4. [Contributing](#contributing)
5. [License](#license)

## Introduction

vFind is unlike a traditional [*variant caller*](https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/).
It is actually using a simpler algorithm which is *usually* sufficient for
screening experiments. The main use case is finding variants from a library that
has constant adapter sequences flanking a variable region.

This simple algorithm is summarized as:

1. Define a pair of adapter sequences that flank the variable region.
2. For each fastq read, search for exact matches of these adapters.
3. If both adapters are found exactly, recover the variable region.
4. For each adapter without an exact match, perform semi-global alignment between the given adapter and read (see the [alignment parameters](#using-custom-alignment-parameters) section
for more details).
5. If the alignment score meets a set threshold, that adapter is considered to match.
6. If both adapters are exactly or partially matched, recover the variable region.
7. For exact matches of both adapters, recover the variable region. Otherwise, continue to the next read.
8. Finally, translate the variable region to its amino acid sequence and filter out any sequences with partial codons (see the [miscellaneous](#miscellaneous) section for more details).

> [!WARNING]
> Note that vFind doesn't do any kind of preprocessing. For initial quality
> filtering, merging, and other common preprocessing operations, you might be
> interested in something like [fastp](https://github.com/OpenGene/fastp) or
> [ngmerge](https://github.com/jsh58/NGmerge). We generally recommend using
> fastp for quality filtering and merging fastq files before using vFind.

Installation details and usage examples are given below. For more usage details,
please see the [API reference](docs/api-reference.md)

## Installation

The package is available on [PyPI](https://pypi.org/project/vfind) and can be installed via pip (or alternatives like [uv](https://github.com/astral-sh/uv)).

### PyPI (Recommended)

Below is an example using uv to initialize a project and add vfind as a dependency.

```bash
uv init
uv add vfind
```

and with pip after creating and activating a new virtual environment

```bash
python3 -m venv .venv
source .venv/bin/activate

python3 -m pip install vfind

```

### Source

vFind is developed using PyO3 and Rust. You will need to make sure you have
a Rust toolchain installed as well as standard C tooling to build some
dependencies (i.e., parasail-rs crate).

1. Clone the repository

```bash
git clone https://github.com/nsbuitrago/vfind.git
cd vfind
```

2. Inside the vfind directory, sync dependencies with uv and build vfind

```bash
uv sync

# this will build and install vfind in the virtual env
uv run maturin develop --uv # or `make dev`
```

## Examples

### Basic Usage

```python
from vfind import find_variants
import polars as pl # variants are returned in a polars dataframe

adapters = ("GGG", "CCC") # define the adapters
fq_path = "./path/to/your/fastq/file.fq.gz" # path to fq file

variants = find_variants(fq_path, adapters)

# print the number of unique sequences
print(variants.n_unique())
```

`find_variants` returns a polars dataframe with `sequence` and `count` columns.
`sequence` contains the amino acid sequence of the variable regions and
`count` contains the frequency of those variant.

We can then use [dataframe methods](https://docs.pola.rs/py-polars/html/reference/dataframe/index.html)
to further analyze the recovered variants. Some examples are shown below.

```python
# Get the top 5 most frequent variants
variants.sort("count", descending=True) # sort by the counts in descending order
print(variants.head(5)) # print the first 5 (most frequent) variants

# filter out sequences with less than 10 read counts
# also any sequences that have a pre-mature stop codon (i.e., * before the last residue)

filtered_variants = variants.filter(
    variants["count"] > 10,
    ~variants["sequence"][::-2].str.contains("*")
)

# write the filtered variants to a csv file
filtered_variants.write_csv("filtered_variants.csv")
```

### Using Custom Alignment Parameters

By default, vFind uses semi-global alignment with the following parameters:

- match score = 3
- mismatch score = -2
- gap open penalty = 5
- gap extend penalty = 2

Note that the gap penalties are represented as positive integers. This is largely due to how the underlying
alignment library works.

To adjust these alignment parameters, use the `match_score`, `mismatch_score`,
`gap_open_penalty`, and `gap_extend_penalty` keyword arguments:

```python
from vfind import find_variants

# ... define adapters and fq_path

# use identity scoring with no gap penalties for alignments
variants = find_variants(
    fq_path,
    adapters,
    match_score = 1,
    mismatch_score = -1,
    gap_open_penalty: 0,
    gap_extend_penalty: 0,
)
```

Alignments are accepted if they produce a score above a set threshold. The
threshold for considering an acceptable alignment can be adjusted with the
`accept_prefix_alignment` and `accept_suffix_alignment` arguments. By default,
both thresholds are set to 0.75.

The thresholds are represent a percentage of the maximum alignment score. So, a
value of 0.75 means alignments producing scores that are greater than 75% the
maximum theoretical score will be accepted. Thus, valid values are between 0 and
1.

Either an exact match or partial match (accepted alignment) must be made for
both adapter sequences to recover a variant.  In order to skip alignment and
only look for exact matches, set the alignment thresolds to 1 (i.e.,
`accept_suffix_alignment = 1` to only allow perfect matches of the suffix
adapter).

### Miscellaneous

**Q:** I don't need the amino acid sequence. Can I just get the DNA sequence?

**A:** Yes. Just set `skip_translation` to True.

```python
# ...
dna_seqs = find_variants(fq_path, adapters, skip_translation=True)
```

---

**Q:** I don't want to use polars. Can I use pandas instead?

**A:** Yes. Use the [`to_pandas`](https://docs.pola.rs/py-polars/html/reference/dataframe/api/polars.DataFrame.to_pandas.html#polars.DataFrame.to_pandas) method on the dataframe.

---

**Q:** I have a lot of data and `find_variants` is slow. Is there anything I can do to speed it up?

**A:** Maybe. Try changing the number of threads or queue length the function uses.

```python
# ...
variants = find_variants(fq_path, adapters, n_threads=6, queue_len=4)
```

For more usage details, see the [API reference](./docs/api-reference.md).

## Contributing

Feedback is a gift and contributions are more than welcome. Please submit an
issue or pull request for any bugs, suggestions, or feedback. Please see the
[developing](./docs/developing.md) guide for more details on how to work on vFind.

## License

vFind is licensed under the [MIT license](./LICENSE)
