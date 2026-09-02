#!/usr/bin/env python3
"""Build the id -> structure catalogue the decoy pool is drawn from.

    python3 tools/build_catalog.py chembl_36_chemreps.txt.gz catalog.sqlite

A similarity search returns identifiers and scores; it has no idea what a molecule is.
This is the lookup that turns those integers back into structures, kept as SQLite so a
2.8M-row catalogue costs nothing to open and does not have to be held in memory.

Vendored here so `dhfrcamp decoys` has a working setup path in one repository.
"""

from __future__ import annotations

import argparse
import gzip
import sqlite3
import sys
import time


def rows(path: str):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        try:
            id_col = header.index("chembl_id")
            smiles_col = header.index("canonical_smiles")
        except ValueError:
            sys.exit(f"{path}: expected chembl_id and canonical_smiles, got {header}")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > smiles_col and parts[smiles_col]:
                yield int(parts[id_col].removeprefix("CHEMBL")), parts[smiles_col]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="chembl_NN_chemreps.txt[.gz]")
    parser.add_argument("output", help="destination .sqlite file")
    args = parser.parse_args()

    started = time.time()
    connection = sqlite3.connect(args.output)
    connection.execute("PRAGMA journal_mode = OFF")
    connection.execute("PRAGMA synchronous = OFF")
    connection.execute("DROP TABLE IF EXISTS molecule")
    # INTEGER PRIMARY KEY aliases the rowid, so the lookup is a b-tree seek with no
    # secondary index to build or carry.
    connection.execute("CREATE TABLE molecule (id INTEGER PRIMARY KEY, smiles TEXT NOT NULL)")
    connection.executemany("INSERT OR REPLACE INTO molecule VALUES (?, ?)", rows(args.input))
    connection.commit()
    count = connection.execute("SELECT count(*) FROM molecule").fetchone()[0]
    connection.close()
    print(f"wrote {count:,} molecules to {args.output} in {time.time() - started:.1f}s")


if __name__ == "__main__":
    main()
