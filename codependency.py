"""Co-dependency analysis helpers and Plotly builders for CRISPR Score Analyzer.

Provides genome-wide single-gene co-dependency scanning (including cohort-scoped
scans unavailable on the DepMap portal), pairwise gene scatter data with
cell-line identity, and theme-aware Plotly visualizations.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from analysis_core import (
    _cell_line_ids,
    _cell_line_names,
    _numeric_values,
    find_lineage_column,
    match_genes,
)

FONT_FAMILY = "Inter, Helvetica Neue, Arial, sans-serif"
ESSENTIALITY_THRESHOLD = -0.5

SCAN_COLUMNS = ["gene", "correlation", "n_shared", "mean_score", "rank", "percentile"]
SCATTER_COLUMNS = ["cell_line", "cell_line_name", "lineage", "score_x", "score_y"]

# Refined, colorblind-conscious editorial palette for scatter-plot lineages.
LINEAGE_PALETTE = [
    "#1b5e91", "#c4583b", "#1d8a68", "#8f5585", "#c98a2c", "#4d90b8",
    "#4b4e8f", "#a3475c", "#3e8c84", "#6b8e4e", "#948447", "#74a5c7",
]


def _default_theme():
    return {
        "plot_bg": "#ffffff",
        "plot_text": "#182028",
        "plot_axis": "#36414d",
        "plot_reference": "rgba(54,65,77,0.28)",
        "plot_grid": "rgba(54,65,77,0.06)",
        "accent": "#1b4d78",
    }


def _ui(lang: str, en: str, zh: str) -> str:
    return zh if lang == "zh" else en


def _pair_correlation(x: np.ndarray, y: np.ndarray, method: str = "pearson") -> float:
    """Exact Pearson or Spearman correlation for finite 1D arrays without scipy."""
    if len(x) < 2:
        return float("nan")
    xv = np.asarray(x, dtype=np.float64)
    yv = np.asarray(y, dtype=np.float64)
    if xv.min() == xv.max() or yv.min() == yv.max():
        return float("nan")
    if method == "spearman":
        xv = pd.Series(xv).rank(method="average").to_numpy(dtype=np.float64)
        yv = pd.Series(yv).rank(method="average").to_numpy(dtype=np.float64)
    dx = xv - xv.mean()
    dy = yv - yv.mean()
    denom = float(np.sqrt(np.dot(dx, dx) * np.dot(dy, dy)))
    if denom <= 0.0:
        return float("nan")
    return float(np.clip(np.dot(dx, dy) / denom, -1.0, 1.0))


def _pearson_chunk(x: np.ndarray, valid_x: np.ndarray, values_sub: np.ndarray, min_shared: int):
    """Exact pairwise-complete Pearson r between ``x[valid_x]`` and ``values_sub``."""
    x_v = x[valid_x]
    n_rows, width = values_sub.shape
    finite = np.isfinite(values_sub)
    n = finite.sum(axis=0)
    r = np.full(width, np.nan, dtype=np.float64)

    complete = n == n_rows
    if complete.any() and n_rows >= min_shared and x_v.min() < x_v.max():
        y_c = values_sub[:, complete]
        nonconst = y_c.max(axis=0) > y_c.min(axis=0)
        dx = x_v - x_v.mean()
        dx_ss = float(np.dot(dx, dx))
        if dx_ss > 0 and nonconst.any():
            dy = y_c - y_c.mean(axis=0, keepdims=True)
            cov = dx @ dy
            denom = np.sqrt(dx_ss * np.einsum("ij,ij->j", dy, dy))
            r_c = np.full(y_c.shape[1], np.nan, dtype=np.float64)
            valid_denom = nonconst & (denom > 0)
            np.divide(cov, denom, out=r_c, where=valid_denom)
            np.clip(r_c, -1.0, 1.0, out=r_c)
            r[complete] = r_c

    inc = (~complete) & (n >= min_shared)
    if inc.any():
        y_i = values_sub[:, inc]
        m_i = finite[:, inc]
        n_i = n[inc]
        x_masked_hi = np.where(m_i, x_v[:, None], -np.inf).max(axis=0)
        x_masked_lo = np.where(m_i, x_v[:, None], np.inf).min(axis=0)
        y_masked_hi = np.where(m_i, y_i, -np.inf).max(axis=0)
        y_masked_lo = np.where(m_i, y_i, np.inf).min(axis=0)
        nonconst = (x_masked_hi > x_masked_lo) & (y_masked_hi > y_masked_lo)
        x_m = np.where(m_i, x_v[:, None], 0.0)
        y_m = np.where(m_i, y_i, 0.0)
        mx = x_m.sum(axis=0) / n_i
        my = y_m.sum(axis=0) / n_i
        dx = np.where(m_i, x_v[:, None] - mx[None, :], 0.0)
        dy = np.where(m_i, y_i - my[None, :], 0.0)
        cov = (dx * dy).sum(axis=0)
        denom = np.sqrt((dx * dx).sum(axis=0) * (dy * dy).sum(axis=0))
        r_i = np.full(y_i.shape[1], np.nan, dtype=np.float64)
        valid_denom = nonconst & (denom > 0)
        np.divide(cov, denom, out=r_i, where=valid_denom)
        np.clip(r_i, -1.0, 1.0, out=r_i)
        r[inc] = r_i
    return r, n.astype(int)


def _spearman_chunk(x_valid: np.ndarray, x_ranks: np.ndarray, values_sub: np.ndarray, min_shared: int):
    """Exact pairwise-complete Spearman rho between ``x_valid`` and ``values_sub``."""
    n_rows, width = values_sub.shape
    finite = np.isfinite(values_sub)
    n = finite.sum(axis=0).astype(int)
    r = np.full(width, np.nan, dtype=np.float64)
    complete = n == n_rows
    if complete.any() and n_rows >= min_shared and x_ranks.min() < x_ranks.max():
        y_c = values_sub[:, complete]
        nonconst = y_c.max(axis=0) > y_c.min(axis=0)
        y_ranks = pd.DataFrame(y_c).rank(axis=0, method="average").to_numpy(dtype=np.float64)
        dx = x_ranks - x_ranks.mean()
        dx_ss = float(np.dot(dx, dx))
        if dx_ss > 0 and nonconst.any():
            dy = y_ranks - y_ranks.mean(axis=0, keepdims=True)
            denom = np.sqrt(dx_ss * np.einsum("ij,ij->j", dy, dy))
            cov = dx @ dy
            r_complete = np.full(y_ranks.shape[1], np.nan, dtype=np.float64)
            valid_denom = nonconst & (denom > 0)
            np.divide(cov, denom, out=r_complete, where=valid_denom)
            np.clip(r_complete, -1.0, 1.0, out=r_complete)
            r[complete] = r_complete
    incomplete_indices = np.flatnonzero((~complete) & (n >= min_shared))
    for col_idx in incomplete_indices:
        col_mask = finite[:, col_idx]
        r[col_idx] = _pair_correlation(x_valid[col_mask], values_sub[col_mask, col_idx], method="spearman")
    return r, n


def genomewide_codependency(
    df: pd.DataFrame,
    query_gene: str,
    rankings: pd.DataFrame,
    method: str = "pearson",
    min_shared: int | None = None,
) -> pd.DataFrame:
    """Correlate one query gene against every other ranked gene in ``df``.

    Uses exact pairwise-complete Pearson or Spearman correlations over finite
    observations. Returns genes ordered from strongest positive co-dependency
    to strongest negative correlation, excluding ``query_gene`` itself.
    """
    method = str(method).lower()
    if method not in {"pearson", "spearman"}:
        raise ValueError("Correlation method must be 'pearson' or 'spearman'")
    matched, _ = match_genes(rankings, [query_gene])
    if not matched:
        raise ValueError(f"Query gene not found in current rankings: {query_gene}")
    canonical = matched[0]
    query_upper = canonical.upper()
    query_row = rankings.loc[rankings["gene_upper"] == query_upper].iloc[0]
    x_series, _ = _numeric_values(df[query_row["gene_raw"]])
    x = x_series.to_numpy(dtype=np.float64)
    valid_x = np.isfinite(x)
    n_valid_x = int(valid_x.sum())
    if min_shared is None:
        min_shared = min(10, max(3, n_valid_x))
    if n_valid_x < min_shared:
        return pd.DataFrame(columns=SCAN_COLUMNS)

    targets = rankings.loc[rankings["gene_upper"] != query_upper].reset_index(drop=True)
    if targets.empty:
        return pd.DataFrame(columns=SCAN_COLUMNS)

    col_positions = {col: idx for idx, col in enumerate(df.columns)}
    raw_positions = np.array([col_positions[raw] for raw in targets["gene_raw"]], dtype=int)
    pos_order = np.argsort(raw_positions, kind="stable")
    ordered_positions = raw_positions[pos_order].tolist()
    unique_dtypes = set(df.dtypes.iloc[ordered_positions])
    use_fast = all(
        pd.api.types.is_numeric_dtype(dt) and not pd.api.types.is_bool_dtype(dt)
        for dt in unique_dtypes
    )

    width = len(ordered_positions)
    corrs_ordered = np.full(width, np.nan, dtype=np.float64)
    shared_ordered = np.zeros(width, dtype=int)

    if use_fast:
        x_valid = x[valid_x]
        x_ranks = (
            pd.Series(x_valid).rank(method="average").to_numpy(dtype=np.float64)
            if method == "spearman" else None
        )
        all_valid_x = bool(valid_x.all())
        valid_rows = None if all_valid_x else np.flatnonzero(valid_x)
        for start in range(0, width, 2048):
            stop = min(start + 2048, width)
            chunk = ordered_positions[start:stop]
            contiguous = chunk[-1] - chunk[0] + 1 == len(chunk)
            if contiguous:
                block = df.iloc[:, chunk[0]:chunk[-1] + 1] if all_valid_x else df.iloc[valid_rows, chunk[0]:chunk[-1] + 1]
            else:
                block = df.iloc[:, chunk] if all_valid_x else df.iloc[valid_rows, chunk]
            values_sub = block.to_numpy(dtype=np.float64)
            if method == "pearson":
                r_chunk, n_chunk = _pearson_chunk(x, valid_x, values_sub, min_shared)
            else:
                r_chunk, n_chunk = _spearman_chunk(x_valid, x_ranks, values_sub, min_shared)
            corrs_ordered[start:stop] = r_chunk
            shared_ordered[start:stop] = n_chunk
    else:
        for idx, pos in enumerate(ordered_positions):
            y_series, _ = _numeric_values(df.iloc[:, pos])
            y_arr = y_series.to_numpy(dtype=np.float64)
            both = valid_x & np.isfinite(y_arr)
            count = int(both.sum())
            shared_ordered[idx] = count
            if count >= min_shared:
                corrs_ordered[idx] = _pair_correlation(x[both], y_arr[both], method=method)

    corrs = np.empty(width, dtype=np.float64)
    shared = np.empty(width, dtype=int)
    corrs[pos_order] = corrs_ordered
    shared[pos_order] = shared_ordered

    result = pd.DataFrame({
        "gene": targets["gene"].to_numpy(),
        "correlation": corrs,
        "n_shared": shared,
        "mean_score": targets["mean_score"].to_numpy(dtype=float),
        "rank": targets["rank"].to_numpy(dtype=int),
        "percentile": targets["percentile"].to_numpy(dtype=float),
    })
    result = result.dropna(subset=["correlation"])
    if result.empty:
        return pd.DataFrame(columns=SCAN_COLUMNS)
    result = result.sort_values(
        ["correlation", "gene"], ascending=[False, True], kind="stable"
    ).reset_index(drop=True)
    return result[SCAN_COLUMNS]


def gene_pair_scatter_frame(
    df: pd.DataFrame,
    gene_x: str,
    gene_y: str,
    rankings: pd.DataFrame,
    lineage_col=None,
) -> tuple[pd.DataFrame, dict]:
    """Build pairwise finite cell-line observations and correlation statistics."""
    lineage_col = find_lineage_column(df) if lineage_col is None else lineage_col
    lookup = rankings.set_index("gene_upper")["gene_raw"].to_dict()
    mx, _ = match_genes(rankings, [gene_x])
    my, _ = match_genes(rankings, [gene_y])
    if not mx or not my:
        return pd.DataFrame(columns=SCATTER_COLUMNS), {
            "n_shared": 0, "pearson_r": np.nan, "spearman_rho": np.nan,
            "slope": np.nan, "intercept": np.nan,
        }
    sx, _ = _numeric_values(df[lookup[mx[0].upper()]])
    sy, _ = _numeric_values(df[lookup[my[0].upper()]])
    lineage = (
        df[lineage_col].astype("string").str.strip().replace("", pd.NA).fillna("Unknown")
        if lineage_col is not None and lineage_col in df.columns
        else pd.Series("Unknown", index=df.index, dtype="string")
    )
    ids = _cell_line_ids(df)
    names = _cell_line_names(df, ids)
    frame = pd.DataFrame({
        "cell_line": ids.to_numpy(),
        "cell_line_name": names.to_numpy(),
        "lineage": lineage.to_numpy(),
        "score_x": sx.to_numpy(dtype=float),
        "score_y": sy.to_numpy(dtype=float),
    }).dropna(subset=["score_x", "score_y"]).reset_index(drop=True)

    n = len(frame)
    xv = frame["score_x"].to_numpy(dtype=float)
    yv = frame["score_y"].to_numpy(dtype=float)
    stats = {
        "n_shared": n,
        "pearson_r": _pair_correlation(xv, yv, method="pearson") if n >= 2 else np.nan,
        "spearman_rho": _pair_correlation(xv, yv, method="spearman") if n >= 2 else np.nan,
        "slope": np.nan,
        "intercept": np.nan,
    }
    if n >= 2 and xv.min() < xv.max():
        var_x = float(np.var(xv))
        if var_x > 0:
            slope = float(np.cov(xv, yv, ddof=0)[0, 1] / var_x)
            stats["slope"] = slope
            stats["intercept"] = float(yv.mean() - slope * xv.mean())
    return frame[SCATTER_COLUMNS], stats


def _apply_theme(fig: go.Figure, th: dict) -> go.Figure:
    fig.update_layout(
        plot_bgcolor=th["plot_bg"],
        paper_bgcolor=th["plot_bg"],
        font=dict(family=FONT_FAMILY, size=13, color=th["plot_text"]),
        hoverlabel=dict(font=dict(family=FONT_FAMILY, size=13)),
        legend_font_color=th["plot_text"],
        title_font_color=th["plot_text"],
    )
    fig.update_xaxes(
        linecolor=th["plot_axis"], linewidth=1.1,
        ticks="outside", ticklen=4, tickwidth=1, tickcolor=th["plot_axis"],
        tickfont=dict(size=12, color=th["plot_text"]),
        title_font=dict(size=13, color=th["plot_text"]),
    )
    fig.update_yaxes(
        linecolor=th["plot_axis"], linewidth=1.1,
        ticks="outside", ticklen=4, tickwidth=1, tickcolor=th["plot_axis"],
        tickfont=dict(size=12, color=th["plot_text"]),
        title_font=dict(size=13, color=th["plot_text"]),
    )
    fig.update_annotations(font_color=th["plot_text"])
    return fig


def create_codependency_scan_plot(
    scan_df: pd.DataFrame,
    query_gene: str,
    method: str = "pearson",
    top_n: int = 12,
    th: dict | None = None,
    lang: str = "en",
) -> go.Figure | None:
    """Diverging horizontal lollipop of top positive and negative co-dependencies."""
    if scan_df.empty:
        return None
    th = th or _default_theme()
    is_dark = th["plot_bg"] != "#ffffff"
    pos = scan_df[scan_df["correlation"] > 0].head(top_n)
    neg = scan_df[scan_df["correlation"] < 0].tail(top_n)
    shown = pd.concat([neg, pos], ignore_index=True).drop_duplicates(subset=["gene"])
    if shown.empty:
        shown = scan_df.head(top_n)
    shown = shown.sort_values(["correlation", "gene"], ascending=[True, True])
    order = shown["gene"].tolist()
    vals = shown["correlation"].to_numpy(dtype=float)
    limit = max(float(np.nanmax(np.abs(vals))), 0.15)
    pos_col = "#5598db" if is_dark else "#1b5e91"
    neg_col = "#d67458" if is_dark else "#c4583b"
    neutral = "#1d2836" if is_dark else "#f7f4ef"
    metric_label = "Pearson r" if method == "pearson" else "Spearman ρ"

    stems = {pos_col: ([], []), neg_col: ([], [])}
    for gene, val in zip(order, vals):
        xs, ys = stems[pos_col if val >= 0 else neg_col]
        xs.extend((0, val, None))
        ys.extend((gene, gene, None))
    fig = go.Figure()
    for color, (xs, ys) in stems.items():
        if not xs:
            continue
        fig.add_trace(go.Scatter(
            x=xs, y=ys, mode="lines",
            line=dict(color=color, width=1.8),
            hoverinfo="skip", showlegend=False,
        ))
    fig.add_trace(go.Scatter(
        x=vals, y=order, mode="markers+text",
        text=[f"  {v:+.2f}" if v >= 0 else f"{v:+.2f}  " for v in vals],
        textposition=["middle right" if v >= 0 else "middle left" for v in vals],
        textfont=dict(size=12, color=th["plot_text"]),
        cliponaxis=False,
        marker=dict(
            size=10.5, color=vals,
            colorscale=[[0.0, neg_col], [0.5, neutral], [1.0, pos_col]],
            cmid=0, cmin=-limit, cmax=limit, showscale=False, line=dict(width=0),
        ),
        customdata=np.column_stack((shown["n_shared"], shown["mean_score"], shown["rank"])),
        hovertemplate=(
            f"<b>%{{y}}</b> × {query_gene}<br>"
            + metric_label + ": %{x:+.3f}<br>"
            + _ui(lang, "Shared cell lines", "共有有效细胞系") + ": %{customdata[0]:.0f}<br>"
            + _ui(lang, "Partner mean score", "伙伴基因平均分") + ": %{customdata[1]:.3f}<br>"
            + _ui(lang, "Partner rank", "伙伴基因排名") + ": %{customdata[2]:,.0f}<extra></extra>"
        ),
        showlegend=False,
    ))
    fig.add_vline(x=0, line_width=1.3, line_color=th["plot_axis"])
    quiet_grid = "rgba(182,195,207,0.12)" if is_dark else "rgba(54,65,77,0.08)"
    fig.update_yaxes(
        automargin=True, categoryorder="array", categoryarray=order,
        showgrid=True, gridcolor=quiet_grid, gridwidth=1,
    )
    fig.update_xaxes(
        range=[-limit * 1.42, limit * 1.42], zeroline=False, showgrid=False,
        title_text=_ui(
            lang,
            f"← anti-correlated  ·  {metric_label} with {query_gene}  ·  co-dependent →",
            f"← 负相关  ·  与 {query_gene} 的 {metric_label}  ·  正共依赖 →",
        ),
    )
    fig.update_layout(
        title=dict(
            text=_ui(
                lang,
                f"<b>{query_gene}</b>  ·  Top co-dependencies ({len(scan_df):,} genes scanned)",
                f"<b>{query_gene}</b>  ·  全基因组共依赖排名（{len(scan_df):,} 基因）",
            ),
            font=dict(size=13.5), x=0, xanchor="left",
        ),
        height=max(420, len(order) * 26 + 120),
        margin=dict(l=24, r=38, t=50, b=58),
        showlegend=False,
    )
    themed = _apply_theme(fig, th)
    themed.update_yaxes(showline=False, ticks="")
    return themed


def create_gene_pair_scatter(
    scatter_df: pd.DataFrame,
    stats: dict,
    gene_x: str,
    gene_y: str,
    th: dict | None = None,
    lang: str = "en",
) -> go.Figure | None:
    """Cell-line scatter plot comparing two genes, colored by cancer lineage."""
    if scatter_df.empty:
        return None
    th = th or _default_theme()
    fig = go.Figure()
    counts = scatter_df["lineage"].value_counts()
    top_lineages = counts.head(10).index.tolist()
    has_other = len(counts) > 10
    hover = (
        "<b>%{customdata[0]}</b> (%{customdata[1]})<br>"
        + _ui(lang, "Cancer type", "癌种") + ": %{customdata[2]}<br>"
        + f"{gene_x}: %{{x:.3f}}<br>"
        + f"{gene_y}: %{{y:.3f}}<extra></extra>"
    )
    # Draw unhighlighted background lineages first so top lineages stay crisp in foreground.
    if has_other:
        other = scatter_df[~scatter_df["lineage"].isin(top_lineages)]
        fig.add_trace(go.Scatter(
            x=other["score_x"], y=other["score_y"], mode="markers",
            name=_ui(lang, f"Other lineages ({len(other)})", f"其他癌种 ({len(other)})"),
            legendrank=100,
            marker=dict(size=5.2, opacity=0.28, color="#9aa4ae", line=dict(width=0.4, color=th["plot_bg"])),
            customdata=other[["cell_line_name", "cell_line", "lineage"]].values,
            hovertemplate=hover,
        ))
    for idx, lineage in enumerate(top_lineages):
        sub = scatter_df[scatter_df["lineage"] == lineage]
        fig.add_trace(go.Scatter(
            x=sub["score_x"], y=sub["score_y"], mode="markers",
            name=f"{lineage} ({len(sub)})",
            legendrank=idx + 1,
            marker=dict(
                size=6.8, opacity=0.82,
                color=LINEAGE_PALETTE[idx % len(LINEAGE_PALETTE)],
                line=dict(width=0.6, color=th["plot_bg"]),
            ),
            customdata=sub[["cell_line_name", "cell_line", "lineage"]].values,
            hovertemplate=hover,
        ))

    if np.isfinite(stats.get("slope", np.nan)) and np.isfinite(stats.get("intercept", np.nan)):
        xs = np.array([float(scatter_df["score_x"].min()), float(scatter_df["score_x"].max())])
        ys = stats["slope"] * xs + stats["intercept"]
        fig.add_trace(go.Scatter(
            x=xs, y=ys, mode="lines",
            name=_ui(lang, "Linear trend", "线性趋势线"),
            legendrank=200,
            line=dict(color=th["plot_axis"], width=1.6, dash="dash"),
            hoverinfo="skip",
        ))

    fig.add_vline(x=0, line=dict(color=th["plot_reference"], width=0.7))
    fig.add_hline(y=0, line=dict(color=th["plot_reference"], width=0.7))
    fig.add_vline(x=ESSENTIALITY_THRESHOLD, line=dict(dash="dot", color="#9aa4ae", width=0.9))
    fig.add_hline(y=ESSENTIALITY_THRESHOLD, line=dict(dash="dot", color="#9aa4ae", width=0.9))

    r_str = f"{stats['pearson_r']:+.3f}" if np.isfinite(stats.get("pearson_r", np.nan)) else "NA"
    rho_str = f"{stats['spearman_rho']:+.3f}" if np.isfinite(stats.get("spearman_rho", np.nan)) else "NA"
    title = (
        f"<b>{gene_x} vs {gene_y}</b>  ·  "
        f"Pearson r = {r_str}  ·  Spearman ρ = {rho_str}  ·  n = {stats['n_shared']:,}"
    )
    fig.update_layout(
        title=dict(text=title, font=dict(size=13), x=0, xanchor="left"),
        xaxis=dict(title=f"{gene_x}  ·  {_ui(lang, 'CRISPR score', 'CRISPR 分数')}",
                   showgrid=True, gridcolor=th["plot_grid"], automargin=True),
        yaxis=dict(title=f"{gene_y}  ·  {_ui(lang, 'CRISPR score', 'CRISPR 分数')}",
                   showgrid=True, gridcolor=th["plot_grid"], automargin=True),
        legend=dict(font=dict(size=11), itemsizing="constant"),
        width=780, height=520,
        margin=dict(l=64, r=24, t=48, b=56),
    )
    return _apply_theme(fig, th)
