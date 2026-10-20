#!/usr/bin/env python3
"""Rename paired gzip-compressed FASTQ reads and replace R2 qualities.

Both mates receive the same name, starting with read_1 (zero-padded to the
number of digits in the total pair count). R1 qualities are preserved.
Inputs are read twice to determine padding and validate them before writing.
Only the Python standard library is required. Wrapped FASTQ is supported;
output records always use four lines and a bare '+' separator.
"""

import argparse
import gzip
import os
from pathlib import Path
import sys
import tempfile
from itertools import zip_longest


def read_fastq(path):
    """Yield sequence and quality strings, validating each FASTQ record."""
    with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
        record = 0
        while True:
            header = handle.readline()
            if not header:
                return
            record += 1
            context = f"{path}: record {record}"
            if not header.startswith("@"):
                raise ValueError(f"{context}: expected a header starting with '@'")

            sequence_parts = []
            while True:
                line = handle.readline()
                if not line:
                    raise ValueError(f"{context}: missing '+' separator")
                if line.startswith("+"):
                    break
                sequence_parts.append(line.rstrip("\r\n"))
            if not sequence_parts:
                raise ValueError(f"{context}: missing sequence line")
            sequence = "".join(sequence_parts)

            quality_parts = []
            quality_length = 0
            while not quality_parts or quality_length < len(sequence):
                line = handle.readline()
                if not line:
                    raise ValueError(f"{context}: truncated quality values")
                quality = line.rstrip("\r\n")
                quality_parts.append(quality)
                quality_length += len(quality)
            if quality_length != len(sequence):
                raise ValueError(f"{context}: sequence and quality lengths differ")
            yield sequence, "".join(quality_parts)


def paired_records(r1, r2):
    for index, pair in enumerate(zip_longest(read_fastq(r1), read_fastq(r2)), 1):
        if any(record is None for record in pair):
            raise ValueError(f"R1 and R2 have different read counts (pair {index})")
        yield pair


def default_output(path):
    name = path.name
    for suffix in (".fastq.gz", ".fq.gz", ".gz"):
        if name.lower().endswith(suffix):
            name = name[:-len(suffix)]
            break
    return path.with_name(name + ".renamed.fastq.gz")


def run(args):
    r1, r2 = args.r1.resolve(), args.r2.resolve()
    out1 = (args.out_r1 or default_output(args.r1)).resolve()
    out2 = (args.out_r2 or default_output(args.r2)).resolve()
    if r1 == r2:
        raise ValueError("R1 and R2 must be different input files")
    if out1 == out2 or out1 in (r1, r2) or out2 in (r1, r2):
        raise ValueError("Output paths must be distinct and different from input paths")
    for path in (out1, out2):
        if path.exists():
            raise ValueError(f"Refusing to overwrite existing output: {path}")
        if not path.parent.is_dir():
            raise ValueError(f"Output directory does not exist: {path.parent}")

    print("Counting and validating read pairs...", file=sys.stderr)
    total = sum(1 for _ in paired_records(r1, r2))
    digits = len(str(total))
    temporary_paths = []
    try:
        for destination in (out1, out2):
            with tempfile.NamedTemporaryFile(
                prefix=".fastq-", suffix=".gz", dir=destination.parent, delete=False
            ) as temporary:
                temporary_paths.append(Path(temporary.name))

        with gzip.open(
            temporary_paths[0], "wt", encoding="ascii", newline="\n",
            compresslevel=args.compression_level,
        ) as f1, gzip.open(
            temporary_paths[1], "wt", encoding="ascii", newline="\n",
            compresslevel=args.compression_level,
        ) as f2:
            written = 0
            for index, ((seq1, qual1), (seq2, _)) in enumerate(paired_records(r1, r2), 1):
                name = f"read_{index:0{digits}d}"
                f1.write(f"@{name}\n{seq1}\n+\n{qual1}\n")
                f2.write(f"@{name}\n{seq2}\n+\n{args.quality * len(seq2)}\n")
                written = index
            if written != total:
                raise ValueError("Input read count changed between passes; retry with stable files")

        for temporary, destination in zip(temporary_paths, (out1, out2)):
            os.replace(temporary, destination)
    finally:
        for path in temporary_paths:
            path.unlink(missing_ok=True)

    print(f"Wrote {total:,} read pairs:\n  {out1}\n  {out2}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("r1", type=Path, help="Input R1 FASTQ.gz")
    parser.add_argument("r2", type=Path, help="Input R2 FASTQ.gz; paired by record order")
    parser.add_argument("--out-r1", type=Path, help="Output R1 (default: INPUT.renamed.fastq.gz)")
    parser.add_argument("--out-r2", type=Path, help="Output R2 (default: INPUT.renamed.fastq.gz)")
    parser.add_argument("--quality", default="I", help="Single R2 quality character (default: I)")
    parser.add_argument(
        "--compression-level", type=int, choices=range(10), default=9,
        help="Gzip level for both outputs: 0 disables compression, 1 is fastest, "
             "9 compresses most (default: 9)",
    )
    args = parser.parse_args()
    if len(args.quality) != 1 or not 33 <= ord(args.quality) <= 126:
        parser.error("--quality must be one printable ASCII character between '!' and '~'")
    try:
        run(args)
    except (OSError, EOFError, UnicodeError, ValueError) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    main()
