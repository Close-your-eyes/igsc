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

This computes and stores every pairwise distance in memory. Large repertoires
need a sparse or chunked workflow instead; CSV output can also become very large.
The first calculation can take longer while Numba compiles the metric.

The script can also be imported:

```python
from cdr3_tcrdist import calculate_distances

distances = calculate_distances(["CASSLGQETQYF", "CASSLAQETQYF"])
```

Upstream references:
- https://tcrdist3.readthedocs.io/en/latest/tcrdistances.html
- https://github.com/agartland/pwseqdist/blob/master/pwseqdist/nb_metrics.py
