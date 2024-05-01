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
3. For exact matches of both adapters, recover the variable region.
3. in case of no exact match for either adapter, perform a semi-global alignment of
the adapter sequence on the fastq read (optional, see the [alignment parameters](#using-custom-alignment-parameters) section).
4. if the alignment is produces a score that is above a set threshold, recover the variable region.
To adjust the thresholds for accepting alignments, see the [alignment parameters](#using-custom-alignment-parameters) section
5. translate the variable region to an amino acid sequence (optional, see the [miscellaneous](#miscellaneous) section).
Note that the translation is done in-frame and assumes the sequence is a coding sequence.

> [!WARNING]
> Note that vFind doesn't do any kind of preprocessing. For initial quality
> filtering, merging, and other common preprocessing operations, you might be
> interested in something like [fastp](https://github.com/OpenGene/fastp) or
> [ngmerge](https://github.com/jsh58/NGmerge). We generally recommend using
> fastp for quality filtering and merging fastq files before using vFind.

Installation details and usage examples are given below. For more usage details,
please see the [API reference](docs/api-reference.md)

## Installation

### PyPI (Recommended for most)

vFind is available on [PyPI](https://pypi.org/) and can be installed via pip (or alternatives like
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
enviroments using [nix flakes](https://wiki.nixos.org/wiki/Flakes).

For something quick, you can use nix-shell. For example, the following will
create a new shell with Python 3.11, vFind, and polars installed.

```bash
nix-shell -p python311 python3Packages.vfind python3Packages.polars
```

## Examples

### Basic Usage

```python
from vfind import find_variants
import polars as pl # variants are returned in a polars dataframe

adapters = ("GGG", "CCC") # define the adapters
fq_path = "./path/to/your/fastq/file.fq.gz" # path to fq file

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

Only alignments that produce scores meeting a threshold will be considered accepted. 
The threshold for considering an acceptable alignment can be adjusted with the
`accept_prefix_alignment` and `accept_suffix_alignment` arguments. By default,
both thresholds are set to 0.75.

The thresholds represent a fraction of the maximum alignment score. So, a value of 0.75
means alignments producing scores that are greater than 75% the maximum theoretical score
will be accepted.

Either an exact match or partial match (accepted alignment) must be made for both adapter sequences to recover a variant. 
In order to skip alignment and only look for exact matches, set the `skip_alignment` argument to `True`.

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

For more usage details, see the [API reference](docs/api-reference.md).

## Contributing

Feedback is a gift and contributions are more than welcome. Please submit an
issue or pull request for any bugs, suggestions, or feedback. Please see the 
[developing](docs/developing) guide for more details on how to work on vFind.

## License

vFind is licensed under the [MIT license](LICENSE)

