# Native 10x TCR annotation on macOS

This standalone FASTA annotator calls the original Cell Ranger annotation and
feature functions. It does not assemble reads, infer clonotypes, or call cells.
It runs locally without uploading sequences.

The included executable was built and tested on Apple Silicon macOS using
Rust 1.88.0. No Rust installation is needed to run it. Keep the executable and
the adjacent `reference/` directory together.

## Run

```sh
./tenx-tcr-annotate \
  --input /Users/chris/Documents/R_packages/igsc/inst/extdata/tcr_test.fasta \
  --output ../tcr_10x_annotations.csv \
  --json ../tcr_10x_annotations.json
```

Run this command from this directory, or call the executable by its absolute
path. The default reference is resolved relative to the executable, independently
of your working directory. Use `--reference /path/to/fasta/regions.fa` to supply
another 10x V(D)J reference. `--help` lists arguments.

The CSV uses the 46 columns in the supplied reference CSV, with four additional
columns: `sequence_id`, `sequence`, `reverse_complemented`, `annotation_notes`.
The JSON exposes the detailed gene alignments and productivity criteria.

Coordinates are zero-based and end-exclusive on the oriented output `sequence`.
If the reverse complement produces the better annotation, the output contains
that oriented sequence and sets `reverse_complemented` to `true`.

The CDR3 sequence includes the conserved flanking residues according to the
10x output convention. Framework four includes the trailing partial codon in
its nucleotide sequence, while its amino-acid sequence includes complete codons
only, consistent with the supplied Cell Ranger reference annotations.

`reads` and `umis` are blank because they cannot be inferred from FASTA. Existing
clonotype and consensus labels are parsed from identifiers matching
`clonotypeN_consensus_M`; this does not infer a biological clonotype.
Other FASTA identifiers are retained in `sequence_id`.

## Validation

The supplied seven TCR sequences were compared with `tcr_annotations_ref.csv`:

- All 44 reference columns other than `reads` and `umis` matched for every row:
  308 matching values.
- All seven regions matched in amino-acid sequence, nucleotide sequence, start,
  and end: 196 matching values.
- All seven reverse-complemented inputs recovered the same oriented sequences
  and annotations.
- Ambiguous-base input was rejected without overwriting an existing output.
- A sequence with no receptor alignment produced an unannotated result.

The machine-readable results are in `validation.json`. To repeat validation:

```sh
python3 validate.py \
  --binary ./tenx-tcr-annotate \
  --fasta /Users/chris/Documents/R_packages/igsc/inst/extdata/tcr_test.fasta \
  --reference-csv /Users/chris/Documents/R_packages/igsc/inst/extdata/tcr_annotations_ref.csv
```

Matching these examples is not a claim that this sequence-only wrapper reproduces
the complete Cell Ranger pipeline on every dataset or every software release.

## Implementation and scope

`vendor/` contains ten unchanged 10x crates from Cell Ranger 10.1.0, commit
`669395e208db7ce03354091e271d893074294d28`:

https://github.com/10XGenomics/cellranger/tree/669395e208db7ce03354091e271d893074294d28/lib/rust

`app/src/main.rs` is the wrapper. It uses:

- `RefData::from_fasta` to load a standard 10x reference.
- `ContigAnnotation::from_seq` for gene annotation, CDR3 and productivity.
- The `fr1_start`, `cdr1_start`, `fr2_start`, `cdr2_start`, and `fr3_start` functions
  on translated V references, as in `enclone_process/src/populate_features.rs`.
- The same region slicing convention as `enclone_proto/src/types.rs`.

The wrapper performs a direct feature transfer for V alignments without indels.
When the V alignment contains an insertion or deletion, upstream framework/CDR1/
CDR2 features are left blank with an explanation in `annotation_notes`. Cell
Ranger's full clonotype-consensus and indel-processing stages are not included.
Feature extraction is limited to TRA and TRB. Input must be unambiguous A/C/G/T
DNA; ambiguous bases are rejected because of the upstream packed-DNA behavior.

Read counts, UMIs, cell confidence, sequencing quality and clonotype assignments
are not inferred. Placeholder values required by the upstream API are excluded
from the JSON output.

The included human reference is 10x V(D)J reference 7.1.0. Its source is recorded
in `reference/SOURCE.txt`. SHA-256 of `reference/regions.fa`:

`5b9d03db44c1cdd66212942f334929ba3babd62e7797d93086779e84a0181e57`

## Rebuild

With a Rust toolchain supporting edition 2024 installed:

```sh
./build.sh
```

The build downloads standard Rust dependencies on its first run. `Cargo.lock`
pins the versions tested here. An Intel Mac requires rebuilding the executable
for that architecture. The upstream source files are unchanged; the workspace
manifest includes only the required crates and dependencies.

## Upstream license

The 10x license is retained in `vendor/LICENSE.10x`. It allows internal use in
connection with 10x products or associated data and prohibits redistribution or
sublicensing. This local project retains the upstream notices and is not a
publicly redistributed package. 10x does not support this wrapper.
