"""Build the bundled DepMap Chronos matrix from the official release files.

The portal files are not committed. Download DepMap Public CRISPRGeneEffect.csv
and the matching Model.csv, then run this script. The app reads the resulting
float32 matrix from disk and does not download it when a page opens.

DepMap's own terms cover the matrix. The MIT license covers the application code.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from analysis_core import GENE_COLUMN_PATTERN, builtin_provenance_path, file_sha256


def checksums(path: Path) -> dict:
    md5 = hashlib.md5()
    sha = hashlib.sha256()
    size = 0
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            md5.update(chunk)
            sha.update(chunk)
            size += len(chunk)
    return {"bytes": size, "md5": md5.hexdigest(), "sha256": sha.hexdigest()}


def build(effect_csv: Path, model_csv: Path, out: Path, release: str) -> dict:
    effect_info = checksums(effect_csv)
    model_info = checksums(model_csv)
    model = pd.read_csv(model_csv, low_memory=False)
    required = {"ModelID", "OncotreeLineage"}
    missing = required.difference(model.columns)
    if missing:
        raise SystemExit(f"Model.csv is missing columns: {', '.join(sorted(missing))}")
    effect = pd.read_csv(effect_csv, index_col=0, low_memory=False)
    if effect.index.duplicated().any():
        raise SystemExit("CRISPRGeneEffect.csv has duplicate ModelID rows")
    if effect.columns.duplicated().any():
        raise SystemExit("CRISPRGeneEffect.csv has duplicate gene columns")
    sample = effect.iloc[:, : min(100, effect.shape[1])].to_numpy(dtype=np.float64)
    scores = effect.to_numpy(dtype=np.float32)
    reconstruction = scores[:, : sample.shape[1]].astype(np.float64)
    difference = np.abs(sample - reconstruction)
    comparable = np.isfinite(sample) & np.isfinite(reconstruction)
    difference[~comparable] = 0
    max_abs_error = float(difference.max()) if difference.size else 0.0
    genes = [str(column).strip() for column in effect.columns]
    unmatched = [gene for gene in genes if GENE_COLUMN_PATTERN.fullmatch(gene) is None]
    if unmatched:
        raise SystemExit(f"{len(unmatched)} gene headers do not match SYMBOL (Entrez ID); first: {unmatched[:5]}")
    ids = effect.index.astype(str)
    meta = model.drop_duplicates("ModelID").set_index("ModelID").reindex(ids)
    if "CellLineName" in meta.columns:
        names = meta["CellLineName"]
    elif "StrippedCellLineName" in meta.columns:
        names = meta["StrippedCellLineName"]
    else:
        names = pd.Series(index=meta.index, dtype="string")
    names = names.astype("string").fillna("").str.strip()
    names = names.mask(names.eq(""), pd.Series(ids.tolist(), index=names.index, dtype="string"))
    lineage = meta["OncotreeLineage"].astype("string").fillna("").str.strip()
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out,
        scores=np.ascontiguousarray(scores, dtype=np.float32),
        genes=np.array(genes, dtype="U128"),
        depmap_id=np.array(ids.tolist(), dtype="U32"),
        cell_line_display_name=np.array(names.tolist(), dtype="U128"),
        lineage=np.array(lineage.tolist(), dtype="U64"),
    )
    provenance = {
        "release": release,
        "score_type": "Chronos Gene Effect",
        "source": "https://depmap.org/portal/",
        "retrieved": "2026-11-21",
        "files": {
            "CRISPRGeneEffect.csv": effect_info,
            "Model.csv": model_info,
        },
        "matrix_sha256": file_sha256(out),
        "n_models": int(scores.shape[0]),
        "n_genes": int(scores.shape[1]),
        "score_dtype": "float32",
        "float32_max_abs_error_first_100_genes": max_abs_error,
        "lineage_column": "OncotreeLineage",
        "name_column": "CellLineName" if "CellLineName" in model.columns else "ModelID",
        "terms": "DepMap portal terms cover this matrix. The MIT license covers only the application code.",
    }
    builtin_provenance_path(out).write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    return provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--effect", type=Path, required=True, help="Official CRISPRGeneEffect.csv")
    parser.add_argument("--model", type=Path, required=True, help="Official Model.csv from the same release")
    parser.add_argument("--out", type=Path, default=Path("data/depmap_chronos.npz"))
    parser.add_argument("--release", default="DepMap Public 26Q3")
    args = parser.parse_args()
    provenance = build(args.effect, args.model, args.out, args.release)
    print(json.dumps({
        "release": provenance["release"],
        "n_models": provenance["n_models"],
        "n_genes": provenance["n_genes"],
        "matrix_bytes": args.out.stat().st_size,
        "matrix_sha256": provenance["matrix_sha256"],
        "source_md5": provenance["files"]["CRISPRGeneEffect.csv"]["md5"],
        "float32_max_abs_error_first_100_genes": provenance["float32_max_abs_error_first_100_genes"],
    }, indent=2))


if __name__ == "__main__":
    main()
