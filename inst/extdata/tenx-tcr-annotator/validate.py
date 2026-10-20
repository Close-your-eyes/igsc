#!/usr/bin/env python3
"""Integration checks against an independently supplied Cell Ranger CSV."""
import argparse
import csv
import json
import subprocess
import tempfile
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--binary', type=Path, required=True)
parser.add_argument('--fasta', type=Path, required=True)
parser.add_argument('--reference-csv', type=Path, required=True)
parser.add_argument('--report', type=Path)
args = parser.parse_args()
binary = args.binary.resolve()

def read_csv(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))

def run(fasta, output):
    return subprocess.run([str(binary), '--input', str(fasta), '--output', str(output)], capture_output=True, text=True)

with tempfile.TemporaryDirectory() as folder:
    temp = Path(folder)
    actual_csv = temp / 'forward.csv'
    result = run(args.fasta.resolve(), actual_csv)
    assert result.returncode == 0, result.stderr
    actual = read_csv(actual_csv)
    expected = {row['consensus_id']: row for row in read_csv(args.reference_csv)}
    assert len(actual) == len(expected), 'Different row counts'
    assert len({row['consensus_id'] for row in actual}) == len(actual), 'Duplicate output identifiers'
    keys = [key for key in next(iter(expected.values())) if key not in ('reads', 'umis')]
    discrepancies = []
    for row in actual:
        for key in keys:
            if row[key] != expected[row['consensus_id']][key]:
                discrepancies.append({'sequence_id': row['sequence_id'], 'column': key, 'expected': expected[row['consensus_id']][key], 'actual': row[key]})
        assert row['reads'] == row['umis'] == '', 'Inferred unsupported read/UMI metadata'
        seq = row['sequence']
        for region in ('fwr1', 'cdr1', 'fwr2', 'cdr2', 'fwr3', 'cdr3', 'fwr4'):
            if row[region + '_start']:
                assert seq[int(row[region + '_start']):int(row[region + '_end'])] == row[region + '_nt'], (row['sequence_id'], region)
    assert not discrepancies, json.dumps(discrepancies, indent=2)

    # Reverse-complementing complete inputs must preserve annotations and
    # coordinates on the restored forward-oriented sequence.
    reverse_fasta = temp / 'reverse.fasta'
    complement = str.maketrans('ACGT', 'TGCA')
    reverse_fasta.write_text(''.join('>' + row['sequence_id'] + '\n' + row['sequence'].translate(complement)[::-1] + '\n' for row in actual))
    reverse_csv = temp / 'reverse.csv'
    result = run(reverse_fasta, reverse_csv)
    assert result.returncode == 0, result.stderr
    reverse = read_csv(reverse_csv)
    assert len(reverse) == len(actual)
    for forward_row, reverse_row in zip(actual, reverse):
        assert reverse_row['reverse_complemented'] == 'true'
        for key in keys + ['sequence']:
            assert forward_row[key] == reverse_row[key], (forward_row['sequence_id'], key)

    # A rejected ambiguous-base record must not overwrite an existing output.
    invalid = temp / 'ambiguous.fasta'
    invalid.write_text('>ambiguous\n' + 'ACGT' * 12 + 'N\n')
    sentinel = temp / 'preserved.csv'
    sentinel.write_text('do not overwrite\n')
    result = run(invalid, sentinel)
    assert result.returncode != 0
    assert sentinel.read_text() == 'do not overwrite\n'

    nohit = temp / 'nohit.fasta'
    nohit.write_text('>nohit\n' + 'A' * 100 + '\n')
    nohit_csv = temp / 'nohit.csv'
    result = run(nohit, nohit_csv)
    assert result.returncode == 0, result.stderr
    nohit_row = read_csv(nohit_csv)[0]
    assert nohit_row['cdr3'] == '' and nohit_row['full_length'] == 'false'

    report = {
        'sequence_count': len(actual),
        'reference_columns_checked_per_sequence': len(keys),
        'reference_values_matched': len(actual) * len(keys),
        'region_sequence_and_coordinate_values_matched': len(actual) * 28,
        'reverse_complement_checks_passed': len(actual),
        'ambiguous_input_preserves_existing_output': True,
        'non_receptor_sequence_handled': True,
        'excluded_reference_fields': ['reads', 'umis'],
        'discrepancies': discrepancies,
    }
    print(json.dumps(report, indent=2))
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + '\n')
