#!/usr/bin/env python3
"""Reproduce all quantitative statistics reported in the CRISPR Score Analyzer preprint.

Run from the repository root:
    python preprint/reproduce_stats.py

Uses only NumPy, pandas, and the Python standard library (no SciPy) together with
the exact analysis functions from analysis_core.py and codependency.py on
data/depmap_chronos.npz.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis_core import (  # noqa: E402
    analyze_dataset,
    file_sha256,
    get_lineage_frame,
    lineage_preference,
    read_builtin_matrix,
)
from codependency import _pair_correlation, gene_pair_scatter_frame, genomewide_codependency  # noqa: E402


def _betacf(a: float, b: float, x: float) -> float:
    """Continued fraction evaluation of the incomplete beta function (Lentz)."""
    maxit, eps, fpmin = 200, 3e-14, 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, maxit + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def pearson_two_sided_p(r: float, n: int) -> float:
    """Exact two-sided Student-t p-value for Pearson correlation r with df = n - 2."""
    df_ = n - 2
    if df_ <= 0:
        return float("nan")
    if abs(r) >= 1.0:
        return 0.0
    x = 1.0 - r * r
    a, b = df_ / 2.0, 0.5
    bt = math.exp(
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
        + a * math.log(x) + b * math.log(1.0 - x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def ols_partial_correlation(x: pd.Series, y: pd.Series, z: pd.Series) -> tuple[float, int, float]:
    """Compute partial Pearson r(x, y | z) on complete finite cases via OLS residuals."""
    sub = pd.DataFrame({"x": x.astype(float), "y": y.astype(float), "z": z.astype(float)})
    sub = sub.replace([np.inf, -np.inf], np.nan).dropna()
    n_complete = len(sub)
    xv = sub["x"].to_numpy(dtype=np.float64)
    yv = sub["y"].to_numpy(dtype=np.float64)
    zv = sub["z"].to_numpy(dtype=np.float64)
    design = np.column_stack([np.ones(n_complete, dtype=np.float64), zv])
    rx = xv - design @ np.linalg.lstsq(design, xv, rcond=None)[0]
    ry = yv - design @ np.linalg.lstsq(design, yv, rcond=None)[0]
    raw_r = float(np.corrcoef(xv, yv)[0, 1])
    partial_r = float(np.corrcoef(rx, ry)[0, 1])
    return partial_r, n_complete, raw_r


def main() -> None:
    matrix_path = ROOT / "data" / "depmap_chronos.npz"
    prov_path = ROOT / "data" / "depmap_chronos.provenance.json"
    prov = json.loads(prov_path.read_text(encoding="utf-8"))
    actual_sha = file_sha256(matrix_path)
    assert actual_sha == prov["matrix_sha256"], (
        f"Matrix SHA-256 mismatch: {actual_sha} != {prov['matrix_sha256']}"
    )

    df = read_builtin_matrix(matrix_path)
    rankings, _ = analyze_dataset(df)
    col_map = dict(zip(rankings["gene"], rankings["gene_raw"]))
    non_empty_lin = (df["lineage"].astype(str).str.strip() != "").sum()
    n_models = len(df)
    n_complete = int((rankings["n_valid"] == n_models).sum())
    n_partial = int((rankings["n_valid"] < n_models).sum())
    n_lt_1000 = int((rankings["n_valid"] < 1000).sum())
    n_lt_500 = int((rankings["n_valid"] < 500).sum())

    print("=== 1. Matrix Provenance, Coverage & Missingness Audit (Section 3.1) ===")
    print(f"Matrix SHA-256 verified : {actual_sha}")
    print(f"Models x Genes          : {n_models:,} x {len(rankings):,}")
    print(f"Matched OncotreeLineage : {non_empty_lin:,}/{n_models:,} across {df['lineage'].nunique()} lineages")
    print(f"Complete genes (n={n_models}): {n_complete:,} ({100.0 * n_complete / len(rankings):.2f}%)")
    print(f"Partial genes (n<{n_models}) : {n_partial:,} ({100.0 * n_partial / len(rankings):.2f}%)")
    print(f"Min n_valid             : {int(rankings['n_valid'].min())} | n<1000: {n_lt_1000} | n<500: {n_lt_500}")
    print(f"Genes with mean < -0.5  : {(rankings['mean_score'] < -0.5).sum():,}")
    print(f"Mean Chronos range      : {rankings['mean_score'].min():+.3f} to {rankings['mean_score'].max():+.3f}")

    print("\n=== 2. Numerical Equivalence of Chunked Correlation Kernels (Section 3.1) ===")
    mdm2_scan = genomewide_codependency(df, "MDM2", rankings, method="pearson")
    mdm2_p = mdm2_scan.set_index("gene")
    mdm2_sp = genomewide_codependency(df, "MDM2", rankings, method="spearman").set_index("gene")
    mdm2_vec = df[col_map["MDM2"]].to_numpy(dtype=np.float64)
    sample_genes = list(dict.fromkeys(
        list(mdm2_p.head(50).index)
        + list(mdm2_p.tail(50).index)
        + list(rankings.loc[rankings["n_valid"] < n_models, "gene"].head(100))
    ))
    sample_genes = [g for g in sample_genes if g in mdm2_p.index]
    max_diff_p = max(
        abs(float(mdm2_p.loc[g, "correlation"]) - float(df[col_map["MDM2"]].astype(np.float64).corr(df[col_map[g]].astype(np.float64))))
        for g in sample_genes
    )
    sp_subset = sample_genes[:100]
    max_diff_s = 0.0
    for g in sp_subset:
        yv = df[col_map[g]].to_numpy(dtype=np.float64)
        mask = np.isfinite(mdm2_vec) & np.isfinite(yv)
        xr = pd.Series(mdm2_vec[mask]).rank(method="average").to_numpy(dtype=np.float64)
        yr = pd.Series(yv[mask]).rank(method="average").to_numpy(dtype=np.float64)
        ref_rho = float(np.corrcoef(xr, yr)[0, 1])
        max_diff_s = max(max_diff_s, abs(float(mdm2_sp.loc[g, "correlation"]) - ref_rho))
    print(f"Pearson max |delta r|  vs pandas.Series.corr ({len(sample_genes)} genes): {max_diff_p:.3e}")
    print(f"Spearman max |delta rho| vs exact rank corr  ({len(sp_subset)} genes): {max_diff_s:.3e}")

    print("\n=== 3. Pan-Cancer Benchmark Rankings (Section 3.2 / Fig. 2A-B) ===")
    for gene in ["MYC", "PTEN", "KRAS", "SOX10", "BRAF", "ESR1", "TP53", "PLK1", "AURKA", "CDK1", "CCNE1"]:
        row = rankings[rankings["gene"] == gene].iloc[0]
        print(
            f"{gene:8s} mean={row['mean_score']:+.3f}  rank={int(row['rank']):5d}  "
            f"pct={row['percentile']:6.2f}%  n_valid={int(row['n_valid'])}"
        )

    print("\n=== 4. Lineage Preference & Top Cell Lines (Section 3.2 / Fig. 2C-D) ===")
    lin_df = get_lineage_frame(df, ["SOX10", "ESR1"], rankings, "lineage")
    pref_df = lineage_preference(lin_df, threshold=-0.5)
    for gene, target_lins in [("SOX10", ["Skin", "Eye"]), ("ESR1", ["Breast", "Ampulla of Vater"])]:
        for lin in target_lins:
            row = pref_df[(pref_df["gene"] == gene) & (pref_df["lineage"] == lin)].iloc[0]
            print(
                f"{gene:6s} | {lin:18s} n={int(row['n_valid']):2d}  "
                f"mean={row['mean']:+.3f}  delta_mean={row['difference']:+.3f}  "
                f"delta_median={row['median_difference']:+.3f}  pct<-0.5={100.0 * row['dependent_fraction']:.1f}%"
            )

    eye_sox = df[df["lineage"] == "Eye"][["depmap_id", "cell_line_display_name", col_map["SOX10"]]].sort_values(col_map["SOX10"]).head(3)
    print("Top Eye SOX10 models:", [(r["cell_line_display_name"], r["depmap_id"], round(float(r[col_map["SOX10"]]), 3)) for _, r in eye_sox.iterrows()])
    breast_esr = df[df["lineage"] == "Breast"][["depmap_id", "cell_line_display_name", col_map["ESR1"]]].sort_values(col_map["ESR1"]).head(3)
    print("Top Breast ESR1 models:", [(r["cell_line_display_name"], r["depmap_id"], round(float(r[col_map["ESR1"]]), 3)) for _, r in breast_esr.iterrows()])

    print("\n=== 5. Systematic Co-Dependency Decomposition Across Canonical Axes (Section 3.3 / Table 2) ===")
    canonical_pairs = [
        ("MDM2", "TP53", "Lung"),
        ("MDM2", "PPM1D", "Lung"),
        ("BRAF", "MAPK1", "Skin"),
        ("BRAF", "MAP2K1", "Skin"),
        ("BRAF", "SOX10", "Skin"),
        ("SOX10", "MITF", "Skin"),
        ("ESR1", "FOXA1", "Breast"),
        ("ESR1", "GATA3", "Breast"),
        ("CDK4", "CCND1", "Breast"),
        ("KRAS", "RAF1", "Bowel"),
    ]
    print(f"{'Pair':13s} {'Lineage':8s} {'n_pan':>5s} {'r_pan':>7s} {'rho_pan':>7s} {'r_within':>8s} {'n_lin':>5s} {'r_lin':>7s} {'rho_lin':>7s}")
    for g1, g2, lin in canonical_pairs:
        sub = pd.DataFrame({
            "lineage": df["lineage"],
            "x": df[col_map[g1]].astype(float),
            "y": df[col_map[g2]].astype(float),
        }).dropna()
        xv = sub["x"].to_numpy(dtype=np.float64)
        yv = sub["y"].to_numpy(dtype=np.float64)
        r_pan = _pair_correlation(xv, yv, "pearson")
        rho_pan = _pair_correlation(xv, yv, "spearman")
        xr = (sub["x"] - sub.groupby("lineage")["x"].transform("mean")).to_numpy(dtype=np.float64)
        yr = (sub["y"] - sub.groupby("lineage")["y"].transform("mean")).to_numpy(dtype=np.float64)
        r_within = _pair_correlation(xr, yr, "pearson")
        sl = sub[sub["lineage"] == lin]
        r_lin = _pair_correlation(sl["x"].to_numpy(dtype=np.float64), sl["y"].to_numpy(dtype=np.float64), "pearson")
        rho_lin = _pair_correlation(sl["x"].to_numpy(dtype=np.float64), sl["y"].to_numpy(dtype=np.float64), "spearman")
        print(f"{g1 + '-' + g2:13s} {lin:8s} {len(sub):5d} {r_pan:+7.3f} {rho_pan:+7.3f} {r_within:+8.3f} {len(sl):5d} {r_lin:+7.3f} {rho_lin:+7.3f}")

    print("\n=== 6. Pan-Cancer MDM2 Scan & MDM2-TP53 Lineage Stratification (Section 3.3 / Fig. 3A-C) ===")
    for label, sub in [("Top pos", mdm2_scan.head(7)), ("Top neg", mdm2_scan.tail(7).iloc[::-1])]:
        items = ", ".join(f"{r['gene']} ({r['correlation']:+.3f}, n={int(r['n_shared'])})" for _, r in sub.iterrows())
        print(f"{label}: {items}")

    _, pair_stats = gene_pair_scatter_frame(df, "MDM2", "TP53", rankings, "lineage")
    print(
        f"Pooled MDM2 vs TP53: r={pair_stats['pearson_r']:+.3f}, "
        f"rho={pair_stats['spearman_rho']:+.3f}, n={pair_stats['n_shared']}"
    )

    pair = pd.DataFrame({
        "lineage": df["lineage"],
        "MDM2": df[col_map["MDM2"]].astype(float),
        "TP53": df[col_map["TP53"]].astype(float),
    }).dropna()
    m_res = pair["MDM2"] - pair.groupby("lineage")["MDM2"].transform("mean")
    t_res = pair["TP53"] - pair.groupby("lineage")["TP53"].transform("mean")
    print(f"Within-lineage demeaned pooled r (30 lineages, n={len(pair)}): {m_res.corr(t_res):+.3f}")

    lin_corrs = []
    for lin, grp in pair.groupby("lineage"):
        if len(grp) >= 15:
            r_val = float(grp["MDM2"].corr(grp["TP53"]))
            lin_corrs.append((lin, len(grp), r_val))
    lin_corrs.sort(key=lambda x: -x[1])
    neg_n = sum(1 for _, _, rv in lin_corrs if rv < 0)
    by_r = sorted(lin_corrs, key=lambda x: x[2])
    print(
        f"Lineages with n >= 15: {len(lin_corrs)} (negative r: {neg_n}/{len(lin_corrs)}, "
        f"median r = {np.median([x[2] for x in lin_corrs]):+.3f}, "
        f"range {by_r[0][2]:+.3f} [{by_r[0][0]}, n={by_r[0][1]}] to "
        f"{by_r[-1][2]:+.3f} [{by_r[-1][0]}, n={by_r[-1][1]}])"
    )
    for lin, cnt, rv in lin_corrs[:10]:
        print(f"  {lin:26s} n={cnt:3d}  r={rv:+.3f}")

    print("\n=== 7. BRAF Pan-Cancer vs Skin Scan, SOX10 Partial r, and Exact Two-Sided t-Tests (Section 3.3 / Fig. 3D) ===")
    braf_pan = genomewide_codependency(df, "BRAF", rankings, method="pearson").set_index("gene")
    print("Pan-cancer BRAF partners:", {g: round(float(braf_pan.loc[g, "correlation"]), 3) for g in ["MAPK1", "MAP2K1", "MAP2K2", "SOX10", "DUSP4", "SOX5"]})

    skin_df = df[df["lineage"] == "Skin"]
    skin_ranks, _ = analyze_dataset(skin_df)
    braf_skin = genomewide_codependency(skin_df, "BRAF", skin_ranks, method="pearson").copy()
    m_tests = len(braf_skin)
    excluded = len(skin_ranks) - 1 - m_tests
    print(f"Skin cohort: n={len(skin_df)} lines, {m_tests:,} partner genes (n_shared >= 10), {excluded} excluded")

    braf_skin["p_two_sided"] = [
        pearson_two_sided_p(float(r["correlation"]), int(r["n_shared"]))
        for _, r in braf_skin.iterrows()
    ]
    order = np.argsort(braf_skin["p_two_sided"].to_numpy())
    p_sorted = braf_skin["p_two_sided"].to_numpy()[order]
    bh_sorted = np.minimum.accumulate((p_sorted * m_tests / np.arange(1, m_tests + 1))[::-1])[::-1]
    bh_fdr = np.empty_like(bh_sorted)
    bh_fdr[order] = np.clip(bh_sorted, 0.0, 1.0)
    braf_skin["bh_fdr"] = bh_fdr
    braf_skin["bonferroni"] = np.clip(braf_skin["p_two_sided"] * m_tests, 0.0, 1.0)

    skin_indexed = braf_skin.set_index("gene")
    braf_s = skin_df[col_map["BRAF"]]
    sox10_s = skin_df[col_map["SOX10"]]
    print("Partner    raw_r   n_shared  p_two_sided  Bonferroni  BH_FDR   partial_r|SOX10  n_complete")
    for gene in ["MAPK1", "SOX10", "SOX5", "MAP2K1", "MAP2K2", "UBQLN1", "AKIRIN1", "STYXL1"]:
        row = skin_indexed.loc[gene]
        if gene == "SOX10":
            pr_str, nc_str = "    N/A", "   75"
        else:
            pr, nc, _ = ols_partial_correlation(braf_s, skin_df[col_map[gene]], sox10_s)
            pr_str, nc_str = f"{pr:+7.3f}", f"{nc:5d}"
        print(
            f"{gene:9s} {row['correlation']:+7.3f}  {int(row['n_shared']):8d}  "
            f"{row['p_two_sided']:11.2e}  {row['bonferroni']:10.4f}  {row['bh_fdr']:7.4f}  "
            f"{pr_str}  {nc_str}"
        )


if __name__ == "__main__":
    main()
