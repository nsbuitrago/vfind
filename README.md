# vFind

A simple variant finder for NGS data.

## Examples

### Basic Usage

```python
from vfind import find_variants, Adapters, AlignParams
import polars as pl # variants are returned in a polars dataframe

# we want to find all variable regions that are flanked by constant adapters

# define the Adapters
adapters = Adapters("GGG", "CCC")

# now, let's find some variants
fq_path = "./path/to/your/fastq/file.fq.gz" # we support gzipped fastq files
variants = find_variants(fq_path, adapters)

# see the top 5 most frequent variants
print(variants.head(5))
```

variants is a polars dataframe with `sequence` and `count` columns.
`sequence` contains the amino acid sequence of the variable region and
`count` contains the frequency of that variant. 

We can then use [dataframe] methods to perform some additional filtering.

