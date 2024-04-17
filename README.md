# vFind

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

The workflow can be described generally by these steps:

1. define a set of constant adapter sequences that immediately flank a variable region of interest.
2. for each fastq read, search for exact matches of these adapters.
3. in case of no exact match for either adapter, perform a semi-global alignment of
the adapter sequence on the fastq read.
4. if the alignment is acceptable, recover the variable region.
To adjust the thresholds for accepting alignments, see the [alignment parameters](#using-custom-alignment-parameters) section

> [!WARNING]
> Note that vFind doesn't do any kind of fastq preprocessing. For initial quality
> filtering, merging, and other common preprocessing operations, you might be
> interested in something like [fastp]() or [ngmerge](). We generally recommend
> using fastp for quality filtering and merging fastq files before using vFind.

## Installation

### PyPI (Recommended for most)

vFind is available on [PyPI](https://pypi.org/project/vfind/0.1.0/) and can be installed via pip (or alternatives like
[uv](https://github.com/astral-sh/uv)).

Below is an example using pip with Python3 in a new project.

```bash
# create a new virtual env
python3 -m venv .venv # create a new virtual env if haven't already
source .venv/bin/activate # activate the virtual env

python3 -m pip install vfind # install vfind
```

### Nix

vFind is also available as a Python package on [NixPkgs](https://search.nixos.org/packages?). You can declare new
development enviroments using [nix flakes](https://wiki.nixos.org/wiki/Flakes).

For something quick, you can use nix-shell,

```bash
nix-shell -p python311 python3Packages.vfind
```

## Examples

### Basic Usage

```python
from vfind import find_variants
import polars as pl # variants are returned in a polars dataframe

adapters = ("GGG", "CCC") # define the Adapters
fq_path = "./path/to/your/fastq/file.fq.gz" # path fo fq file

# now, let's find some variants
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
    ~variants["counts"] < 10,
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

Note that the gap penalties are represented as positive integers. This is largely due to the underlying
alignment library.

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

Alignments are accepted if they produce a score above a set threshold. The threshold
for considering an acceptable alignment can be adjusted with the `accept_prefix_alignment`
and `accept_suffix_alignment` arguments. By default, both thresholds are set to 0.75.

The thresholds represent a fraction of the maximum alignment score. So, a value of 0.75
means alignments producing scores that are greater than 75% the maximum theoretical score
will be accepted. If an alignment produces a score below the threshold, it will
not be accepted and no sequence will be recovered from that given read.

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

**Q:** `find_variants` is slow. Is there anything I can do to speed it up?

**A:** Maybe? Try changing the number of threads or queue length the function uses.

```python
# ...
variants = find_variants(fq_path, adapters, n_threads=6, queue_len=4)
```

## Contributing

Contributions are welcome. Please submit an issue or pull request for any bugs,
suggestions, or feedback.

### Developing

vFind is written in Rust and uses [PyO3](https://pyo3.rs/v0.21.1/) and [Maturin](https://github.com/PyO3/maturin)
to generate the Python module. To get started, you will need to have the
[rust toolchain](https://www.rust-lang.org/tools/install) and [Python >= 3.10](https://www.python.org/downloads/).

Below are some general steps to get up and running. Note that these examples
use [uv]). However, you could do this with standard pip or your preferred method.

1. clone the repository to your machine

```bash
git clone git@github.com:nsbuitrago/vfind.git
cd vfind
```

2. Create a new Python virtual environment and install dev dependencies.

```bash
uv venv
source .venv/bin/activate

# this will intsall maturin, polars, and any other required packages.
uv pip install -r dev-requirements.txt
```

3. Build and install the vFind package in your virtual environment with maturin. 

```bash
# in the root project directory
maturin develop
```

4. From here, you can make changes to the Rust lib and run `maturin` develop
to rebuild the package with those changes.

## License

vFind is licensed under the [MIT license](LICENSE)

