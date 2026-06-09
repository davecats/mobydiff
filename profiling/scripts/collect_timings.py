#!/usr/bin/env python3
"""Collect per-case summary.csv files into one CSV table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", nargs="?", default="profiling/results")
    args = parser.parse_args()

    root = Path(args.results)
    rows = []
    fieldnames = None
    for summary in sorted(root.rglob("summary.csv")):
        with summary.open(newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                row["summary_file"] = str(summary)
                rows.append(row)
                if fieldnames is None:
                    fieldnames = list(reader.fieldnames or []) + ["summary_file"]

    if not rows:
        return

    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)


if __name__ == "__main__":
    main()
