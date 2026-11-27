"""Tests for the importable Streamlit-free crispr_score_analyzer Python library."""

from __future__ import annotations

import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

import crispr_score_analyzer as csa
from crispr_score_analyzer import (
    analyze_dataset,
    create_correlation_heatmap,
    create_lineage_boxplot,
    create_preference_heatmap,
    create_preference_plot,
    create_rank_plot,
    filter_gene_suggestions,
    gene_correlation_matrix,
    genomewide_codependency,
    get_lineage_frame,
    lineage_preference,
    load_gene_aliases,
    read_score_csv,
    resolve_gene_list,
)


def _synthetic_csv_bytes() -> bytes:
    return (
        b"ModelID,cell_line_display_name,lineage,MYC (4609),PTEN (5728),"
        b"TEAD1 (7003),TEAD2 (8463),ERBB2 (2064),BRAF (673),MAPK1 (5594)\n"
        b"ACH-000001,Skin1,Skin,-1.95,0.28,-1.40,-0.40,-0.12,-1.65,-1.42\n"
        b"ACH-000002,Skin2,Skin,-2.05,0.31,-1.35,-0.38,-0.10,-1.58,-1.35\n"
        b"ACH-000003,Skin3,Skin,-1.88,0.24,-1.28,-0.41,-0.14,-1.50,-1.28\n"
        b"ACH-000004,Skin4,Skin,-1.92,0.29,-1.30,-0.39,-0.09,-1.48,-1.22\n"
        b"ACH-000005,Skin5,Skin,-2.01,0.33,-1.25,-0.35,-0.11,-1.62,-1.39\n"
        b"ACH-000006,Breast1,Breast,-1.90,0.35,-0.20,-0.10,-1.85,-0.08,-0.12\n"
        b"ACH-000007,Breast2,Breast,-2.02,0.26,-0.18,-0.12,-1.72,-0.15,-0.10\n"
        b"ACH-000008,Breast3,Breast,-1.85,0.38,-0.15,-0.09,-1.90,-0.05,-0.07\n"
        b"ACH-000009,Breast4,Breast,-1.94,0.30,-0.22,-0.08,-1.68,-0.11,-0.09\n"
        b"ACH-000010,Breast5,Breast,-1.99,0.27,-0.19,-0.11,-1.80,-0.09,-0.11\n"
    )


def test_core_library_imports_without_streamlit():
    """Verify the package never imports streamlit during library use."""
    code = (
        "import sys; import crispr_score_analyzer as csa; "
        "assert 'streamlit' not in sys.modules; "
        "assert hasattr(csa, 'create_rank_plot'); "
        "assert hasattr(csa, 'genomewide_codependency')"
    )
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=False)
    assert proc.returncode == 0, proc.stderr


def test_public_version_and_aliases_loading():
    assert isinstance(csa.__version__, str) and len(csa.__version__) >= 5
    aliases = load_gene_aliases()
    assert len(aliases) > 10_000
    assert aliases.get("HER2") == "ERBB2"
    assert aliases.get("TEF-1") == "TEAD1"


def test_filter_gene_suggestions_exact_prefix_and_alias_ranking():
    symbols = ["A1BG", "ERBB2", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "STEAP1", "MYC"]
    aliases = {"HER2": "ERBB2", "TEF-1": "TEAD1", "TCF13": "TEAD1"}

    # Exact numbered symbol TEAD1 must rank first and never collapse to A1BG or TEAD2.
    exact = filter_gene_suggestions("TEAD1", symbols, aliases)
    assert [item["gene"] for item in exact] == ["TEAD1"]
    assert exact[0]["matched_via"] is None

    # Prefix TEAD returns TEAD1..TEAD4 in order without fuzzy noise.
    prefix = filter_gene_suggestions("tead", symbols, aliases)
    assert [item["gene"] for item in prefix] == ["TEAD1", "TEAD2", "TEAD3", "TEAD4"]

    # Entrez suffix paste also resolves cleanly.
    entrez = filter_gene_suggestions("TEAD1 (7003)", symbols, aliases)
    assert [item["gene"] for item in entrez] == ["TEAD1"]

    # HGNC aliases resolve cleanly with matched_via provenance.
    alias_hit = filter_gene_suggestions("HER2", symbols, aliases)
    assert alias_hit == [{"gene": "ERBB2", "matched_via": "HER2"}]

    tef_hit = filter_gene_suggestions("TEF-1", symbols, aliases)
    assert tef_hit[0] == {"gene": "TEAD1", "matched_via": "TEF-1"}


def test_end_to_end_streamlit_free_pipeline():
    df = read_score_csv(_synthetic_csv_bytes())
    rankings, diagnostics = analyze_dataset(df)
    assert len(rankings) == 7
    assert diagnostics.empty is False

    aliases = load_gene_aliases()
    matched, missing, resolved = resolve_gene_list(rankings, ["TEAD1", "HER2", "BRAF"], aliases)
    assert matched == ["TEAD1", "ERBB2", "BRAF"]
    assert missing == []
    assert resolved == [{"input": "HER2", "gene": "ERBB2"}]

    plot_df = get_lineage_frame(df, matched, rankings, "lineage")
    pref = lineage_preference(plot_df)
    assert set(pref["gene"]) == {"TEAD1", "ERBB2", "BRAF"}
    for theme in ("light", "dark"):
        rank_fig = create_rank_plot(rankings, matched, n_cell_lines=len(df), theme=theme)
        assert len(rank_fig.data) >= 2
        box_fig = create_lineage_boxplot(plot_df, matched, theme=theme)
        assert len(box_fig.data) >= 2
        pref_fig = create_preference_plot(pref, theme=theme)
        assert pref_fig is not None
        hm_fig = create_preference_heatmap(pref, theme=theme)
        assert hm_fig is not None

    scan = genomewide_codependency(df, "BRAF", rankings, method="pearson", min_shared=5)
    assert scan.iloc[0]["gene"] == "MAPK1"
    assert scan.iloc[0]["correlation"] > 0.9

    corr, counts = gene_correlation_matrix(df, ["TEAD1", "BRAF", "MAPK1"], rankings, method="spearman")
    corr_fig = create_correlation_heatmap(corr, counts, method="spearman", theme="dark")
    assert corr_fig.data[0].type == "heatmap"
    assert np.isfinite(corr.loc["BRAF", "MAPK1"])
