#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 14:09:27 2026

@author: chris
"""

# pip install tcrdist3 pandas

import argparse

import pandas as pd
import os
from pathlib import Path
from tcrdist.repertoire import TCRrep


parser = argparse.ArgumentParser(
    description="Calculate TRA or TRB TCR distance matrices."
)
parser.add_argument(
    # UserWarning: db_file must be 'alphabeta_gammadelta_db.tsv' or 'alphabeta_db.tsv' or 'gammadelta_db.tsv' unless you have built tcrdist3 from scratch self._validate_db_file()
    "--db-file",
    default="alphabeta_gammadelta_db.tsv",
    choices=["alphabeta_db.tsv", " alphabeta_gammadelta_db.tsv", "gammadelta_db"],
    help="Path to the TCRdist database file (default: %(default)s)",
)
parser.add_argument(
    "--chain",
    choices=["TRA", "TRB"],
    default="TRB",
    type=str.upper,
    help="T-cell receptor chain to process (default: %(default)s)",
)
parser.add_argument(
    "--input",
    default="tcrdist3_ex.csv",
    help="Path to the input CSV file (default: %(default)s)",
)
parser.add_argument(
    "--species",
    choices=["human", "mouse"],
    default="human",
    help="Path to the input CSV file (default: %(default)s)",
)
parser.add_argument(
    "--suffix",
    default="",
    help="Optional suffix added to output filenames, before .csv",
)
parser.add_argument(
    "--workdir",
    type=Path,
    default=Path.cwd(),
    help="Directory for input and output files (default: current directory)",
)

args = parser.parse_args()

workdir = args.workdir.expanduser().resolve()

if not workdir.is_dir():
    parser.error(f"Working directory does not exist: {workdir}")

try:
    os.chdir(workdir)
except OSError as error:
    parser.error(f"Cannot use working directory {workdir}: {error}")


# Configure column and attribute names for the selected chain.
if args.chain == "TRA":
    chain = "alpha"
    cdr3_column = "cdr3_a_aa"
    v_column = "v_a_gene"
    j_column = "j_a_gene"
    distance_attribute = "pw_alpha"
    cdr3_distance_attribute = "pw_cdr3_a_aa"
    output_prefix = "alpha"
else:
    chain = "beta"
    cdr3_column = "cdr3_b_aa"
    v_column = "v_b_gene"
    j_column = "j_b_gene"
    distance_attribute = "pw_beta"
    cdr3_distance_attribute = "pw_cdr3_b_aa"
    output_prefix = "beta"


# Load the CSV.
df = pd.read_csv(args.input)

# Check and clean the chain-specific columns.
required = [cdr3_column, v_column, j_column]
missing = [column for column in required if column not in df.columns]

if missing:
    parser.error(
        f"{args.chain} requires these missing CSV columns: "
        + ", ".join(missing)
    )

df = df.dropna(subset=required).copy()

df[cdr3_column] = df[cdr3_column].str.strip().str.upper()
df[v_column] = df[v_column].str.strip()
df[j_column] = df[j_column].str.strip()

if "count" not in df.columns:
    df["count"] = 1


# Construct the repertoire and calculate distances.
tr = TCRrep(
    cell_df=df,
    organism=args.species,
    chains=[chain],
    db_file=args.db_file,
)

labels = tr.clone_df[cdr3_column]

full_distance = pd.DataFrame(
    getattr(tr, distance_attribute),
    index=labels,
    columns=labels,
)

cdr3_distance = pd.DataFrame(
    getattr(tr, cdr3_distance_attribute),
    index=labels,
    columns=labels,
)

suffix = f"_{args.suffix}" if args.suffix else ""

full_distance.to_csv(
    f"tcrdist3_full_distmat_{output_prefix}{suffix}.csv"
)
cdr3_distance.to_csv(
    f"tcrdist3_cdr3_distmat_{output_prefix}{suffix}.csv"
)

# print(full_distance)