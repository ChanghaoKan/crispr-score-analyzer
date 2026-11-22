"""Build data/gene_aliases.json from an HGNC table and the bundled DepMap matrix.

Download the HGNC complete set TSV, then run:

    python scripts/build_gene_aliases.py --hgnc /path/to/hgnc_complete_set.txt

Only previous symbols and aliases that map to exactly one approved symbol in
the bundled matrix are kept. A token that is itself a current matrix symbol
is not redirected.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from datetime import date
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def matrix_symbols(path: Path) -> set[str]:
    genes = np.load(path, allow_pickle=False)["genes"].astype(str)
    return {gene.split(" (", 1)[0].strip().upper() for gene in genes if str(gene).strip()}


def alias_map(hgnc_path: Path, symbols: set[str]) -> dict[str, str]:
    targets = defaultdict(set)
    with hgnc_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            approved = (row.get("symbol") or "").strip().upper()
            if approved not in symbols:
                continue
            for field in ("prev_symbol", "alias_symbol"):
                for token in (row.get(field) or "").split("|"):
                    token = token.strip().upper()
                    if token and token != approved:
                        targets[token].add(approved)
    return {key: next(iter(values)) for key, values in sorted(targets.items())
            if len(values) == 1 and key not in symbols}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hgnc", type=Path, required=True)
    parser.add_argument("--matrix", type=Path, default=ROOT / "data" / "depmap_chronos.npz")
    parser.add_argument("--out", type=Path, default=ROOT / "data" / "gene_aliases.json")
    args = parser.parse_args()
    payload = {
        "source": "HGNC complete set",
        "source_url": "https://www.genenames.org/",
        "retrieved": date.today().isoformat(),
        "note": "Previous symbols and aliases that map to exactly one approved symbol present in the bundled DepMap matrix. Current DepMap symbols are not redirected.",
        "map": alias_map(args.hgnc, matrix_symbols(args.matrix)),
    }
    args.out.write_text(json.dumps(payload, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    print(f"{len(payload['map'])} aliases -> {args.out}")


if __name__ == "__main__":
    main()
