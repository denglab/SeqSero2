#!/usr/bin/env python3
"""Combine SeqSero_result.tsv files using their named columns."""
import argparse
import csv
from pathlib import Path


def summarize(paths, output):
    rows = []
    columns = []
    for source in paths:
        path = Path(source)
        if path.is_dir():
            path = path / 'SeqSero_result.tsv'
        with path.open(newline='', encoding='utf-8') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if not reader.fieldnames or 'Predicted serotype' not in reader.fieldnames:
                raise ValueError('missing TSV header (do not use -s): ' + str(path))
            for column in reader.fieldnames:
                if column not in columns:
                    columns.append(column)
            for row in reader:
                if None in row or any(value is None for value in row.values()):
                    raise ValueError('malformed TSV row: ' + str(path))
                rows.append(row)
    writer = csv.DictWriter(output, fieldnames=columns, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    return len(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reports', nargs='+', help='TSV reports or their output directories')
    parser.add_argument('-o', '--output', required=True, help='new summary TSV (will not overwrite)')
    args = parser.parse_args()
    # Validate and combine before creating the destination, so malformed inputs
    # cannot leave a partial summary that looks complete.
    import io
    buffer = io.StringIO(newline='')
    summarize(args.reports, buffer)
    with open(args.output, 'x', newline='', encoding='utf-8') as output:
        output.write(buffer.getvalue())


if __name__ == '__main__':
    main()
