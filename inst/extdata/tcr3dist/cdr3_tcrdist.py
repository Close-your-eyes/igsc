#!/usr/bin/env python3
"""Calculate CDR3-only pairwise distances from a CSV, without V/J genes.

Install: python3 -m pip install pwseqdist pandas
Run: python3 cdr3_tcrdist.py sequences.csv --cdr3-column cdr3 -o distances.csv

Defaults match the raw CDR3 calculation in the accompanying chat example:
dist_weight=1, gap_penalty=4, ntrim=3, ctrim=2, fixed_gappos=False.
No additional CDR3 weight of three is applied.
"""

import argparse
from pathlib import Path
import sys

import pandas as pd
import pwseqdist as pw


AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def calculate_distances(
    sequences, *, cpus=1, ntrim=3, ctrim=2, gap_penalty=4,
    dist_weight=1, fixed_gappos=False,
):
    """Return a dense matrix in input order, retaining duplicate sequences.

    Accepts an iterable of CDR3 amino-acid strings from either TRA or TRB.
    Whitespace around sequences is removed and letters are uppercased.
    """
    if cpus < 1 or dist_weight < 1:
        raise ValueError("cpus and dist_weight must be positive integers.")
    if min(ntrim, ctrim, gap_penalty) < 0:
        raise ValueError("Trimming values and gap_penalty cannot be negative.")
    cleaned = []
    for position, sequence in enumerate(sequences, start=1):
        if not isinstance(sequence, str) or not sequence.strip():
            raise ValueError(f"Record {position}: CDR3 is missing or empty.")
        sequence = sequence.strip().upper()
        invalid = set(sequence) - AMINO_ACIDS
        if invalid:
            raise ValueError(
                f"Record {position}: unsupported amino-acid characters: "
                + ", ".join(sorted(invalid))
            )
        if len(sequence) < 5 or len(sequence) <= ntrim + ctrim:
            raise ValueError(
                f"Record {position}: CDR3 must have at least five residues "
                "and retain at least one residue after trimming."
            )
        cleaned.append(sequence)
    if not cleaned:
        raise ValueError("No CDR3 sequences were provided.")

    # The upstream metric stores scores in signed 16-bit integers.
    max_length = max(map(len, cleaned))
    max_substitution = int(pw.matrices.tcr_nb_distance_matrix.max())
    if max_length * (max_substitution * dist_weight + gap_penalty) > 32767:
        raise ValueError("Sequence lengths or scoring settings risk integer overflow.")

    return pw.apply_pairwise_rect(
        metric=pw.metrics.nb_vector_tcrdist,
        seqs1=pd.Series(cleaned).to_numpy(),
        ncpus=cpus,
        use_numba=True,
        uniqify=True,
        reexpand=True,
        distance_matrix=pw.matrices.tcr_nb_distance_matrix,
        dist_weight=dist_weight,
        gap_penalty=gap_penalty,
        ntrim=ntrim,
        ctrim=ctrim,
        fixed_gappos=fixed_gappos,
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_csv", type=Path, help="CSV with a header row")
    parser.add_argument("--cdr3-column", default="cdr3", help="CDR3 column (default: cdr3)")
    parser.add_argument("--id-column", help="Optional column of unique record identifiers for labeled output")
    parser.add_argument("--include-labels", action="store_true",
                        help="Include column headers and row identifiers (default: values only)")
    parser.add_argument("-o", "--output", type=Path, help="Output matrix CSV")
    parser.add_argument("--cpus", type=int, default=1, help="Worker count (default: 1)")
    parser.add_argument("--ntrim", type=int, default=3, help="N-terminal trim (default: 3)")
    parser.add_argument("--ctrim", type=int, default=2, help="C-terminal trim (default: 2)")
    parser.add_argument("--gap-penalty", type=int, default=4, help="Gap penalty (default: 4)")
    parser.add_argument("--dist-weight", type=int, default=1,
                        help="Substitution weight before gap penalties (default: 1)")
    parser.add_argument("--fixed-gappos", action="store_true",
                        help="Use fixed gap placement instead of optimized placement")
    args = parser.parse_args(argv)
    output = args.output or args.input_csv.with_name(args.input_csv.stem + ".tcrdist.csv")

    try:
        if output.resolve() == args.input_csv.resolve():
            raise ValueError("Output must be different from the input file.")
        frame = pd.read_csv(args.input_csv, dtype=str, keep_default_na=False)
        if args.cdr3_column not in frame.columns:
            raise ValueError(f"CDR3 column {args.cdr3_column!r} was not found.")
        labels = None
        if args.include_labels and args.id_column:
            if args.id_column not in frame.columns:
                raise ValueError(f"ID column {args.id_column!r} was not found.")
            labels = frame[args.id_column].str.strip().tolist()
            if any(not label for label in labels) or len(set(labels)) != len(labels):
                raise ValueError("Record identifiers must be nonempty and unique.")
        elif args.include_labels:
            labels = [f"row_{i}" for i in range(1, len(frame) + 1)]

        print(f"Calculating pairwise distances for {len(frame)} records...", file=sys.stderr)
        distances = calculate_distances(
            frame[args.cdr3_column], cpus=args.cpus, ntrim=args.ntrim,
            ctrim=args.ctrim, gap_penalty=args.gap_penalty,
            dist_weight=args.dist_weight, fixed_gappos=args.fixed_gappos,
        )
        pd.DataFrame(distances, index=labels, columns=labels).to_csv(
            output, index=args.include_labels, header=args.include_labels,
            index_label="record_id" if args.include_labels else None,
        )
    except (OSError, ValueError, pd.errors.ParserError) as error:
        parser.exit(2, f"Error: {error}\n")
    print(f"Saved {len(frame)} by {len(frame)} distance matrix to {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
