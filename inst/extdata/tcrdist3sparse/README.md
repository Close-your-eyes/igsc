# CDR3-only TCRdist

Install dependencies:

```sh
python3 -m pip install pwseqdist pandas
```

Run from the folder containing the script and example:

```sh
python3 cdr3_tcrdist.py example_cdr3.csv -o distances.csv
```

For a different CDR3 column:

```sh
python3 cdr3_tcrdist.py input.csv --cdr3-column cdr3_b_aa -o distances.csv
```

The input needs a header and a column of CDR3 amino-acid sequences. No V or J
genes are required. Either TRA or TRB sequences work; analyze each chain
separately for within-chain comparisons. Only the 20 standard amino-acid
letters are accepted. Surrounding whitespace is removed and lowercase letters
are converted to uppercase. Missing or invalid sequences produce an error.

The output is a dense CSV containing only distance values, without a header row
or record ID column. Rows and columns follow the input record order. Duplicate
sequences remain separate records. Read this output with `pd.read_csv(path,
header=None)` in pandas.

To include labels, add `--include-labels`. Labels then come from `--id-column`
when supplied; otherwise they are `row_1`, `row_2`, and so on, numbering data
records after the header. IDs used as labels must be nonempty and unique.

Defaults match the earlier CDR3-only example: substitution weight 1, gap penalty
4, N-terminal trim 3, C-terminal trim 2, and optimized gap placement. No additional
CDR3 weight of three is applied. These explicit settings differ from the bare
`nb_vector_tcrdist` function defaults. To score the entire CDR3, add
`--ntrim 0 --ctrim 0`. Sequences must have at least five residues and retain at
least one residue after trimming. Use `--help` to see all scoring options.

The script directly calls `pwseqdist.metrics.nb_vector_tcrdist` through
`pwseqdist.apply_pairwise_rect`, which is also used by tcrdist3's `_pw` wrapper.
It does not need the full tcrdist3 package. Scores depend on the selected scoring
settings; changing the substitution weight does not scale the gap penalty.

Dense mode computes and stores every pairwise distance in memory. CSV output can
also become very large. The first calculation can take longer while Numba
compiles the metric.

## Sparse output for larger repertoires

```sh
python3 cdr3_tcrdist.py input.csv --cdr3-column cdr3 \
  --sparse --radius 12 --chunk-size 100 --cpus 4 -o distances.npz
```

`12` is an example cutoff, not a validated biological threshold. Choose a cutoff
appropriate to your data and scoring settings. `--radius` is required in sparse
mode. Only distances at or below that cutoff are retained. The script calculates
blocks of at most `--chunk-size` rows against all input records, then stores each
block sparsely. This also works with one CPU. Smaller chunks reduce temporary
memory. More workers can increase memory use. All pairs are still evaluated,
so sparse output reduces storage and peak memory but does not remove the growth
in computation with repertoire size. The retained sparse matrix must still fit
in RAM; a large radius can make it effectively dense. Numba worker startup may
add overhead for each block.

No IDs or labels are written in sparse mode. Input row order and duplicate
records are preserved. `--include-labels` is only available for dense CSV output.

The sparse encoding follows tcrdist3:

- Positive stored values are actual distances within the cutoff.
- Stored `-1` means an actual distance of zero, including the diagonal and
  identical sequences.
- Unstored entries, read numerically as zero, mean the distance exceeds the
  cutoff. They are not identical sequences.

Load a Python sparse output with:

```python
from scipy.sparse import load_npz

distances = load_npz("distances.npz")
```

For R, choose the Matrix Market format by using an `.mtx` output filename:

```sh
python3 cdr3_tcrdist.py input.csv --sparse --radius 12 \
  --chunk-size 100 --cpus 4 -o distances.mtx
```

```r
distances <- Matrix::readMM("distances.mtx")
```

This is a thresholded neighbor matrix, not a complete distance matrix. Do not
pass it directly to `as.dist()`, `hclust()`, or `uwot::umap()`: absent entries and
the `-1` encoding need appropriate handling. A UMAP neighbor-list workflow must
decode zero distances, include self as the first neighbor, and provide enough
known neighbors per record. If a record has too few retained neighbors, rerun
with a larger radius or calculate its nearest neighbors explicitly. Keep dense
mode when you need all exact distances for hierarchical clustering.

The script can also be imported:

```python
from cdr3_tcrdist import calculate_distances

distances = calculate_distances(["CASSLGQETQYF", "CASSLAQETQYF"])
```

Upstream references:
- https://tcrdist3.readthedocs.io/en/latest/tcrdistances.html
- https://github.com/agartland/pwseqdist/blob/master/pwseqdist/nb_metrics.py
