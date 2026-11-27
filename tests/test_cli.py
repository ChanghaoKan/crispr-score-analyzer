"""Tests for the crispr-score-analyzer command-line interface on synthetic CSVs."""

from __future__ import annotations

import json
from pathlib import Path
import zipfile

import pandas as pd

from crispr_score_analyzer.cli import main


def _write_synthetic_csv(path: Path) -> Path:
    path.write_text(
        "ModelID,cell_line_display_name,lineage,MYC (4609),PTEN (5728),"
        "TEAD1 (7003),ERBB2 (2064),BRAF (673),MAPK1 (5594)\n"
        "ACH-000001,Skin1,Skin,-1.95,0.28,-1.40,-0.12,-1.65,-1.42\n"
        "ACH-000002,Skin2,Skin,-2.05,0.31,-1.35,-0.10,-1.58,-1.35\n"
        "ACH-000003,Skin3,Skin,-1.88,0.24,-1.28,-0.14,-1.50,-1.28\n"
        "ACH-000004,Skin4,Skin,-1.92,0.29,-1.30,-0.09,-1.48,-1.22\n"
        "ACH-000005,Skin5,Skin,-2.01,0.33,-1.25,-0.11,-1.62,-1.39\n"
        "ACH-000006,Breast1,Breast,-1.90,0.35,-0.20,-1.85,-0.08,-0.12\n"
        "ACH-000007,Breast2,Breast,-2.02,0.26,-0.18,-1.72,-0.15,-0.10\n"
        "ACH-000008,Breast3,Breast,-1.85,0.38,-0.15,-1.90,-0.05,-0.07\n"
        "ACH-000009,Breast4,Breast,-1.94,0.30,-0.22,-1.68,-0.11,-0.09\n"
        "ACH-000010,Breast5,Breast,-1.99,0.27,-0.19,-1.80,-0.09,-0.11\n",
        encoding="utf-8",
    )
    return path


def test_cli_info_reports_synthetic_dataset(tmp_path, capsys):
    csv_path = _write_synthetic_csv(tmp_path / "scores.csv")
    rc = main(["info", "--input", str(csv_path)])
    assert rc == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["package"] == "crispr-score-analyzer"
    assert payload["dataset"]["n_cell_lines"] == 10
    assert payload["dataset"]["n_ranked_genes"] == 6
    assert payload["dataset"]["n_lineages"] == 2


def test_cli_rank_with_cohort_alias_and_bundle(tmp_path, capsys):
    csv_path = _write_synthetic_csv(tmp_path / "scores.csv")
    out_dir = tmp_path / "rank_out"
    rc = main([
        "rank",
        "--input", str(csv_path),
        "--genes", "TEAD1,HER2,BRAF",
        "--cohort", "Skin",
        "--output-dir", str(out_dir),
        "--bundle",
    ])
    assert rc == 0
    captured = capsys.readouterr()
    assert "TEAD1" in captured.out
    assert "ERBB2" in captured.out
    assert "HER2 -> ERBB2" in captured.err

    selected = pd.read_csv(out_dir / "selected_results.csv")
    assert selected["gene"].tolist() == ["TEAD1", "ERBB2", "BRAF"]
    assert (selected["n_valid"] == 5).all()
    assert (out_dir / "full_rankings.csv").is_file()
    assert (out_dir / "lineage_preference.csv").is_file()
    assert (out_dir / "diagnostics.csv").is_file()

    zip_path = out_dir / "analysis_bundle.zip"
    assert zip_path.is_file()
    with zipfile.ZipFile(zip_path) as zf:
        names = set(zf.namelist())
        assert {"selected_results.csv", "full_rankings.csv", "plot_data.csv", "metadata.json", "diagnostics.csv"} <= names


def test_cli_codep_scan_and_matrix(tmp_path, capsys):
    csv_path = _write_synthetic_csv(tmp_path / "scores.csv")
    scan_csv = tmp_path / "braf_scan.csv"
    rc = main([
        "codep",
        "--input", str(csv_path),
        "--query", "BRAF",
        "--method", "pearson",
        "--min-shared", "5",
        "--top", "3",
        "--output", str(scan_csv),
    ])
    assert rc == 0
    out = capsys.readouterr().out
    assert "MAPK1" in out
    scan_df = pd.read_csv(scan_csv)
    assert scan_df.iloc[0]["gene"] == "MAPK1"

    corr_csv = tmp_path / "panel_corr.csv"
    rc = main([
        "codep",
        "--input", str(csv_path),
        "--genes", "TEAD1,BRAF,MAPK1",
        "--method", "spearman",
        "--output", str(corr_csv),
    ])
    assert rc == 0
    corr_df = pd.read_csv(corr_csv, index_col=0)
    assert list(corr_df.columns) == ["TEAD1", "BRAF", "MAPK1"]


def test_cli_returns_nonzero_for_missing_file(tmp_path, capsys):
    rc = main(["rank", "--input", str(tmp_path / "does_not_exist.csv"), "--genes", "MYC"])
    assert rc == 2
    assert "Error:" in capsys.readouterr().err
