"""Unit tests for genome-wide co-dependency scanning and pairwise scatter views."""

import unittest

import numpy as np
import pandas as pd

from analysis_core import analyze_dataset
from codependency import (
    SCAN_COLUMNS,
    SCATTER_COLUMNS,
    create_codependency_scan_plot,
    create_gene_pair_scatter,
    gene_pair_scatter_frame,
    genomewide_codependency,
)
from figure_export import prepare_export_figure


def _sample_matrix():
    """Synthetic CRISPR matrix with positive, negative, missing, and constant profiles."""
    rng = np.random.default_rng(42)
    n = 24
    base = np.linspace(-1.5, 0.2, n)
    pos = base * 0.9 + rng.normal(0, 0.05, n)
    neg = -base * 0.8 + rng.normal(0, 0.05, n)
    indep = np.sin(np.arange(n) * 1.7) * 0.4
    const = np.full(n, -0.3)

    # Introduce NaN and Inf edge cases.
    base_with_nan = base.copy()
    base_with_nan[0] = np.nan
    pos_with_nan = pos.copy()
    pos_with_nan[1] = np.nan
    pos_with_nan[2] = np.inf
    neg_with_ties = neg.copy()
    neg_with_ties[3:6] = -0.5

    df = pd.DataFrame({
        "ModelID": [f"ACH-{i:06d}" for i in range(1, n + 1)],
        "cell_line_display_name": [f"Line {i}" for i in range(1, n + 1)],
        "lineage": ["Melanoma"] * 8 + ["Lung"] * 8 + ["Breast"] * 8,
        "BRAF (673)": base_with_nan,
        "MAPK1 (5594)": pos_with_nan,
        "PTEN (5728)": neg_with_ties,
        "INDEP (1001)": indep,
        "CONST (1002)": const,
    })
    rankings, _ = analyze_dataset(df)
    return df, rankings


class GenomewideCodependencyTests(unittest.TestCase):
    def test_pearson_and_spearman_match_pandas_exact_pairwise(self):
        df, rankings = _sample_matrix()
        for method in ("pearson", "spearman"):
            with self.subTest(method=method):
                scan = genomewide_codependency(df, "BRAF", rankings, method=method, min_shared=5)
                self.assertEqual(list(scan.columns), SCAN_COLUMNS)
                self.assertNotIn("BRAF", set(scan["gene"]))
                self.assertNotIn("CONST", set(scan["gene"]))
                self.assertEqual(scan.iloc[0]["gene"], "MAPK1")
                self.assertEqual(scan.iloc[-1]["gene"], "PTEN")

                x = pd.to_numeric(df["BRAF (673)"], errors="coerce").replace([np.inf, -np.inf], np.nan)
                for _, row in scan.iterrows():
                    raw_col = rankings.loc[rankings["gene"] == row["gene"], "gene_raw"].iloc[0]
                    y = pd.to_numeric(df[raw_col], errors="coerce").replace([np.inf, -np.inf], np.nan)
                    both = x.notna() & y.notna()
                    if method == "pearson":
                        expected_r = float(x[both].corr(y[both], method="pearson"))
                    else:
                        expected_r = float(x[both].rank(method="average").corr(y[both].rank(method="average"), method="pearson"))
                    self.assertEqual(int(row["n_shared"]), int(both.sum()))
                    self.assertAlmostEqual(float(row["correlation"]), expected_r, places=6)

    def test_object_dtype_fallback_matches_numeric_fast_path(self):
        df, rankings = _sample_matrix()
        fast = genomewide_codependency(df, "BRAF", rankings, method="pearson", min_shared=5)

        obj_df = df.copy()
        for col in ["BRAF (673)", "MAPK1 (5594)", "PTEN (5728)", "INDEP (1001)", "CONST (1002)"]:
            obj_df[col] = obj_df[col].map(lambda v: "" if pd.isna(v) else ("inf" if np.isinf(v) else f"{v:.8f}"))
        obj_rankings, _ = analyze_dataset(obj_df)
        slow = genomewide_codependency(obj_df, "BRAF", obj_rankings, method="pearson", min_shared=5)

        self.assertEqual(list(fast["gene"]), list(slow["gene"]))
        np.testing.assert_allclose(fast["correlation"].to_numpy(), slow["correlation"].to_numpy(), atol=1e-5)
        self.assertEqual(list(fast["n_shared"]), list(slow["n_shared"]))

    def test_min_shared_and_error_validation(self):
        df, rankings = _sample_matrix()
        empty = genomewide_codependency(df, "BRAF", rankings, method="pearson", min_shared=50)
        self.assertTrue(empty.empty)
        self.assertEqual(list(empty.columns), SCAN_COLUMNS)

        with self.assertRaises(ValueError):
            genomewide_codependency(df, "BRAF", rankings, method="kendall")
        with self.assertRaises(ValueError):
            genomewide_codependency(df, "UNKNOWN_GENE", rankings, method="pearson")

    def test_cohort_scoped_scan_reflects_subset_correlation(self):
        df = pd.DataFrame({
            "ModelID": [f"ACH-{i:06d}" for i in range(1, 13)],
            "lineage": ["Skin"] * 6 + ["Lung"] * 6,
            "GENE_A (1)": [-1.2, -1.0, -0.8, -0.6, -0.4, -0.2, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2],
            "GENE_B (2)": [-1.1, -0.9, -0.7, -0.5, -0.3, -0.1, -0.1, -0.3, -0.5, -0.7, -0.9, -1.1],
        })
        skin_df = df[df["lineage"] == "Skin"]
        lung_df = df[df["lineage"] == "Lung"]
        skin_rankings, _ = analyze_dataset(skin_df)
        lung_rankings, _ = analyze_dataset(lung_df)

        skin_scan = genomewide_codependency(skin_df, "GENE_A", skin_rankings, min_shared=4)
        lung_scan = genomewide_codependency(lung_df, "GENE_A", lung_rankings, min_shared=4)
        self.assertGreater(float(skin_scan.iloc[0]["correlation"]), 0.99)
        self.assertLess(float(lung_scan.iloc[0]["correlation"]), -0.99)


class GenePairScatterAndPlotTests(unittest.TestCase):
    def test_scatter_frame_and_stats(self):
        df, rankings = _sample_matrix()
        scatter_df, stats = gene_pair_scatter_frame(df, "BRAF", "MAPK1", rankings)
        self.assertEqual(list(scatter_df.columns), SCATTER_COLUMNS)
        # Rows 0 (BRAF NaN), 1 (MAPK1 NaN), 2 (MAPK1 Inf) excluded -> 21 valid pairs.
        self.assertEqual(stats["n_shared"], 21)
        self.assertEqual(len(scatter_df), 21)
        self.assertGreater(stats["pearson_r"], 0.95)
        self.assertGreater(stats["spearman_rho"], 0.95)
        self.assertGreater(stats["slope"], 0.5)
        self.assertTrue(np.isfinite(stats["intercept"]))
        self.assertEqual(set(scatter_df["lineage"]), {"Melanoma", "Lung", "Breast"})

    def test_scatter_missing_gene_returns_empty(self):
        df, rankings = _sample_matrix()
        scatter_df, stats = gene_pair_scatter_frame(df, "BRAF", "MISSING", rankings)
        self.assertTrue(scatter_df.empty)
        self.assertEqual(stats["n_shared"], 0)
        self.assertTrue(np.isnan(stats["pearson_r"]))

    def test_scan_plot_and_scatter_plot_render_and_export_cleanly(self):
        df, rankings = _sample_matrix()
        scan = genomewide_codependency(df, "BRAF", rankings, method="pearson", min_shared=5)
        scatter_df, stats = gene_pair_scatter_frame(df, "BRAF", "MAPK1", rankings)
        dark_theme = {
            "plot_bg": "#151d26",
            "plot_text": "#ebf0f5",
            "plot_axis": "#b6c3cf",
            "plot_reference": "rgba(182,195,207,0.34)",
            "plot_grid": "rgba(182,195,207,0.09)",
            "accent": "#67a4e0",
        }
        scan_fig = create_codependency_scan_plot(scan, "BRAF", method="pearson", top_n=5, th=dark_theme, lang="zh")
        self.assertIsNotNone(scan_fig)
        self.assertGreaterEqual(len(scan_fig.data), 2)
        self.assertIn("BRAF", scan_fig.layout.title.text)

        scatter_fig = create_gene_pair_scatter(scatter_df, stats, "BRAF", "MAPK1", th=dark_theme, lang="en")
        self.assertIsNotNone(scatter_fig)
        self.assertGreaterEqual(len(scatter_fig.data), 4)  # 3 lineages + linear trendline
        exported = prepare_export_figure(scatter_fig)
        self.assertEqual(exported.layout.paper_bgcolor, "white")
        # Dark marker halos matching dark plot_bg should be reset to white for clean publication export.
        self.assertEqual(exported.data[0].marker.line.color, "white")

        self.assertIsNone(create_codependency_scan_plot(pd.DataFrame(columns=SCAN_COLUMNS), "BRAF"))
        self.assertIsNone(create_gene_pair_scatter(pd.DataFrame(columns=SCATTER_COLUMNS), stats, "BRAF", "MAPK1"))


if __name__ == "__main__":
    unittest.main()
