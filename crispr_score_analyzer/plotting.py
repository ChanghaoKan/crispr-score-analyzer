"""Streamlit-free CNS/publication Plotly builders for CRISPR Score Analyzer.

Provides deterministic, theme-aware (neutral light and charcoal dark) Plotly
figures using colorblind-safe Nature Publishing Group (NPG) categorical hues
and a 7-stop ColorBrewer RdBu diverging scale. All figures can be used directly
from Python scripts, Jupyter notebooks, the CLI, or the Streamlit web interface.
"""

from __future__ import annotations

from html import escape

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


ESSENTIALITY_THRESHOLD = -0.5
PREFERENCE_MIN_N = 5
PREFERENCE_MAX_LOLLIPOP_GENES = 6

PLOT_COLORS = {
    # Colorblind-safe Nature Publishing Group (NPG) / CNS publication palette
    "essential": "#E64B35",
    "nonessential": "#00A087",
    "interest": "#3C5488",
    "threshold": "#7F8C8D",
    "boxplot_fill": "#3C5488",
}

FONT_FAMILY = "-apple-system, BlinkMacSystemFont, Inter, Helvetica Neue, Arial, sans-serif"

PLOT_THEMES = {
    "light": {
        "text": "#1d1d1f",
        "plot_bg": "#ffffff",
        "plot_text": "#1d1d1f",
        "plot_axis": "#424245",
        "plot_reference": "rgba(29,29,31,0.24)",
        "plot_grid": "rgba(29,29,31,0.06)",
        "plot_scatter_bg": "rgba(110,110,115,0.22)",
    },
    "dark": {
        "text": "#f5f5f7",
        "plot_bg": "#1c1c1e",
        "plot_text": "#f5f5f7",
        "plot_axis": "#aeaeb2",
        "plot_reference": "rgba(245,245,247,0.24)",
        "plot_grid": "rgba(255,255,255,0.07)",
        "plot_scatter_bg": "rgba(174,174,178,0.34)",
    },
}


def _ui(lang: str, en: str, zh: str) -> str:
    return zh if lang == "zh" else en


def get_plot_theme(theme: str | dict = "light") -> dict:
    """Resolve a theme name ('light' or 'dark') or custom token dict."""
    if isinstance(theme, dict):
        return theme
    return PLOT_THEMES.get(str(theme).lower(), PLOT_THEMES["light"])


def _is_dark_theme(theme: str | dict) -> bool:
    if isinstance(theme, dict):
        return str(theme.get("plot_bg", "")).lower() not in {"", "#ffffff", "white"}
    return str(theme).lower() == "dark"


def preference_colorscale(theme: str | dict = "light") -> list:
    """ColorBrewer 7-stop RdBu scale: deep steel blue (< 0) through neutral to crimson (> 0)."""
    if _is_dark_theme(theme):
        return [
            [0.0, "#4393C3"], [0.22, "#67A9CF"], [0.38, "#3B5970"],
            [0.5, "#24262B"],
            [0.62, "#6E3F3B"], [0.78, "#D6604D"], [1.0, "#E05A5A"],
        ]
    return [
        [0.0, "#2166AC"], [0.17, "#4393C3"], [0.33, "#92C5DE"],
        [0.5, "#F7F7F7"],
        [0.67, "#F4A582"], [0.83, "#D6604D"], [1.0, "#B2182B"],
    ]


def apply_theme_to_fig(fig: go.Figure, theme: str | dict = "light") -> go.Figure:
    """Apply publication typography, axis ticks, and light/dark surfaces to a figure."""
    th = get_plot_theme(theme)
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


def _rank_label_annotations(points, x_range, y_range, plot_width=552, plot_height=366):
    """Space labels in a conservative viewport, keeping arrows on the true data.

    A small spatial grid bounds collision checks even for a large submitted set.
    Positions use data coordinates so labels and leader lines also remain vector
    objects in exports and scale with the responsive chart.
    """
    x_span, y_span = x_range[1] - x_range[0], y_range[1] - y_range[0]
    occupied = {}
    annotations = []

    def cells(rect):
        left, top, right, bottom = rect
        return [
            (col, row)
            for col in range(int(left // 64), int(right // 64) + 1)
            for row in range(int(top // 64), int(bottom // 64) + 1)
        ]

    offsets = [
        (dx, dy)
        for distance in (24, 48, 72, 96, 120)
        for dx in (0, 48, -48, 96, -96)
        for dy in (-distance, distance)
    ]
    for point in sorted(points, key=lambda item: (item["x"], item["y"], item["gene"])):
        px = (point["x"] - x_range[0]) / x_span * plot_width
        py = (y_range[1] - point["y"]) / y_span * plot_height
        half_width = max(24, len(point["gene"]) * 3.7 + 7)
        half_height = 12
        best = None
        for dx, dy in offsets:
            cx = min(plot_width - half_width - 4, max(half_width + 4, px + dx))
            cy = min(plot_height - half_height - 4, max(half_height + 4, py + dy))
            rect = (cx - half_width, cy - half_height, cx + half_width, cy + half_height)
            nearby = {other for cell in cells(rect) for other in occupied.get(cell, [])}
            overlap = sum(
                max(0, min(rect[2], other[2]) - max(rect[0], other[0]))
                * max(0, min(rect[3], other[3]) - max(rect[1], other[1]))
                for other in nearby
            )
            if best is None or overlap < best[0]:
                best = (overlap, cx, cy, rect)
            if overlap == 0:
                break
        _, cx, cy, rect = best
        for cell in cells(rect):
            occupied.setdefault(cell, []).append(rect)
        annotations.append(dict(
            x=point["x"], y=point["y"], xref="x", yref="y",
            ax=x_range[0] + cx / plot_width * x_span,
            ay=y_range[1] - cy / plot_height * y_span, axref="x", ayref="y",
            text=escape(point["gene"]), showarrow=True, arrowhead=0,
            arrowwidth=0.8, arrowcolor=point["color"], standoff=5,
            xanchor="center", yanchor="middle", borderpad=2,
            bgcolor="rgba(0,0,0,0)", borderwidth=0,
            font=dict(size=13, color=point["color"], family=FONT_FAMILY),
        ))
    return annotations


def _ranking_fingerprint(rankings: pd.DataFrame):
    """Identify a ranking curve without hashing every gene name."""
    if len(rankings) == 0:
        return (0, "", "", 0.0)
    scores = rankings["mean_score"].to_numpy(dtype=np.float64, copy=False)
    return (
        len(rankings),
        str(rankings["gene"].iloc[0]),
        str(rankings["gene"].iloc[-1]),
        round(float(scores.sum()), 5),
    )


def _cached_background_trace(rankings, excluded, color, name, hover):
    """Reuse the full background curve when only labels or marker size change."""
    cache = getattr(_cached_background_trace, "cache", None)
    if cache is None:
        cache = {}
        _cached_background_trace.cache = cache
    key = (_ranking_fingerprint(rankings), tuple(sorted(excluded)), color, name, hover)
    cached = cache.get(key)
    if cached is None:
        upper = rankings["gene_upper"].to_numpy()
        keep = np.ones(len(rankings), dtype=bool) if not excluded else ~np.isin(upper, list(excluded))
        cached = go.Scattergl(
            x=rankings["rank"].to_numpy()[keep],
            y=np.round(rankings["mean_score"].to_numpy(dtype=np.float64)[keep], 4),
            mode="markers",
            marker=dict(size=2.5, color=color),
            name=name,
            text=rankings["gene"].to_numpy()[keep],
            customdata=np.column_stack((
                rankings["n_valid"].to_numpy()[keep],
                np.round(rankings["percentile"].to_numpy(dtype=np.float64)[keep], 2),
            )),
            hovertemplate=hover,
        )
        if len(cache) >= 8:
            cache.clear()
        cache[key] = cached
    return cached


def build_rank_figure(
    rankings: pd.DataFrame,
    layers: list[dict],
    references: list[str],
    n_cell_lines: int,
    point_size: int | float = 4,
    show_reference_labels: bool = True,
    width: int = 640,
    height: int = 480,
    theme: str | dict = "light",
    lang: str = "en",
) -> go.Figure:
    """Build a cohort-wide essentiality ranking scatter plot with layered highlights."""
    th = get_plot_theme(theme)
    fig = go.Figure()
    highlighted = {g.upper() for layer in layers for g in layer["genes"]}
    refs = list(dict.fromkeys(g.upper() for g in references if g))
    hover = (
        "<b>%{text}</b><br>" + _ui(lang, "Rank", "排名") + ": %{x:,}<br>"
        + _ui(lang, "Mean score", "平均分") + ": %{y:.4f}<br>"
        + _ui(lang, "Valid cell lines", "有效细胞系数") + ": %{customdata[0]}<br>"
        + _ui(lang, "Rank percentile (lower = stronger)", "排名百分位（越低依赖越强）")
        + ": %{customdata[1]:.2f}%<extra></extra>"
    )
    fig.add_trace(_cached_background_trace(
        rankings, highlighted | set(refs), th["plot_scatter_bg"], _ui(lang, "All genes", "全部基因"), hover,
    ))
    label_points = []

    def collect_labels(part, color):
        label_points.extend(
            {"x": record.rank, "y": record.mean_score, "gene": record.gene, "color": color}
            for record in part.itertuples(index=False)
        )

    for i, gene in enumerate(refs):
        if gene in highlighted:
            continue
        part = rankings[rankings.gene_upper == gene]
        if part.empty:
            continue
        color = PLOT_COLORS["essential"] if i == 0 else PLOT_COLORS["nonessential"]
        fig.add_trace(go.Scatter(
            x=part["rank"], y=part.mean_score, mode="markers",
            marker=dict(
                size=point_size * 2.6, color=color,
                symbol="diamond" if i == 0 else "square",
                line=dict(width=1, color=th["plot_bg"]),
            ),
            text=part.gene, name=_ui(lang, f"Reference · {gene}", f"参考基因 · {gene}"),
            customdata=part[["n_valid", "percentile"]].values, hovertemplate=hover,
        ))
        if show_reference_labels:
            collect_labels(part, color)

    for layer in layers:
        part = rankings[rankings.gene_upper.isin({g.upper() for g in layer["genes"]})].sort_values("rank")
        if part.empty:
            continue
        fig.add_trace(go.Scatter(
            x=part["rank"], y=part.mean_score,
            mode="markers",
            marker=dict(
                size=point_size * layer.get("size", 2.5), color=layer["color"],
                symbol=layer.get("symbol", "circle"),
                line=dict(width=1, color=th["plot_bg"]),
            ),
            text=part.gene,
            name=layer["name"],
            customdata=part[["n_valid", "percentile"]].values, hovertemplate=hover,
        ))
        if layer["labels"]:
            collect_labels(part, layer["color"])

    fig.add_hline(y=ESSENTIALITY_THRESHOLD, line=dict(dash="dash", color=PLOT_COLORS["threshold"], width=1))
    fig.add_hline(y=0, line=dict(color=th["plot_reference"], width=0.7))
    lo, hi = rankings.mean_score.min(), rankings.mean_score.max()
    span = max(hi - lo, 0.2)
    x_range = [0, max(len(rankings) * 1.02, 1)]
    y_range = [lo - span * 0.12, hi + span * 0.18]
    margin = dict(l=64, r=24, t=42, b=72)
    fig.update_layout(
        xaxis=dict(
            title=_ui(lang, "Gene rank", "基因排名"), tickformat=",d", range=x_range,
            showgrid=False, automargin=True,
        ),
        yaxis=dict(
            title=_ui(lang, "Mean CRISPR score", "平均 CRISPR 分数"), range=y_range,
            showgrid=True, gridcolor=th["plot_grid"], automargin=True,
        ),
        annotations=_rank_label_annotations(
            label_points, x_range, y_range,
            plot_width=max(240, width - margin["l"] - margin["r"]),
            plot_height=max(200, height - margin["t"] - margin["b"]),
        ),
        title=dict(
            text=_ui(
                lang,
                f"{n_cell_lines:,} cell lines  ·  −0.5 reference",
                f"{n_cell_lines:,} 个细胞系  ·  −0.5 参考线",
            ),
            font=dict(size=12.5), xref="paper", x=0, xanchor="left",
        ),
        legend=dict(
            orientation="h", yanchor="top", y=-0.16, xanchor="center", x=0.5,
            font=dict(size=12), itemsizing="constant",
        ),
        width=width, height=height, margin=margin,
    )
    return apply_theme_to_fig(fig, theme=th)


def create_rank_plot(
    gene_rank_df: pd.DataFrame,
    genes_of_interest: list[str],
    essential_gene: str = "MYC",
    nonessential_gene: str = "PTEN",
    n_cell_lines: int = 0,
    show_labels: bool = True,
    point_size: int | float = 4,
    width: int = 640,
    height: int = 480,
    theme: str | dict = "light",
    lang: str = "en",
) -> go.Figure:
    """Create a gene essentiality ranking curve highlighting ``genes_of_interest``."""
    interest_color = "#67A9CF" if _is_dark_theme(theme) else PLOT_COLORS["interest"]
    return build_rank_figure(
        gene_rank_df,
        [{"genes": genes_of_interest, "labels": show_labels, "color": interest_color,
          "name": _ui(lang, "Genes of interest", "目标基因")}],
        [essential_gene, nonessential_gene],
        n_cell_lines,
        point_size,
        show_reference_labels=show_labels,
        width=width,
        height=height,
        theme=theme,
        lang=lang,
    )


def create_multilayer_rank_plot(
    gene_rank_df: pd.DataFrame,
    background_genes: list[str],
    highlight_genes: list[str],
    bg_color: str = "#8E99A2",
    hl_color: str = "#3C5488",
    essential_gene: str = "MYC",
    nonessential_gene: str = "PTEN",
    n_cell_lines: int = 0,
    show_labels: bool = True,
    point_size: int | float | None = 4,
    width: int = 640,
    height: int = 480,
    theme: str | dict = "light",
    lang: str = "en",
) -> go.Figure:
    """Create a two-layer ranking curve overlaying highlight genes onto a background gene set."""
    bg_only = [g for g in background_genes if g not in highlight_genes]
    layers = [
        {"genes": bg_only, "labels": False, "color": bg_color, "size": 1.7,
         "name": _ui(lang, "Background gene set", "背景基因集")},
        {"genes": highlight_genes, "labels": show_labels, "color": hl_color, "symbol": "circle",
         "name": _ui(lang, "Highlight genes", "高亮基因")},
    ]
    if point_size is None:
        point_size = 4
    return build_rank_figure(
        gene_rank_df,
        layers,
        [essential_gene, nonessential_gene],
        n_cell_lines,
        point_size,
        show_reference_labels=show_labels,
        width=width,
        height=height,
        theme=theme,
        lang=lang,
    )


def create_lineage_boxplot(
    lineage_data: pd.DataFrame,
    genes: list[str],
    sort_mode: str = "median",
    show_points: bool = False,
    columns: int = 1,
    compare_mode: str = "by_lineage",
    theme: str | dict = "light",
    lang: str = "en",
) -> go.Figure:
    """Create horizontal per-lineage or cross-gene CRISPR score boxplots."""
    th = get_plot_theme(theme)
    # Horizontal groups keep long cancer-type names readable on narrow screens.
    genes = [g for g in genes if g in set(lineage_data.gene)]
    n_genes = len(genes)
    score_min, score_max = lineage_data.crispr_score.min(), lineage_data.crispr_score.max()
    pad = max((score_max - score_min) * 0.06, 0.15)
    is_dark = _is_dark_theme(theme)

    def box_palette(median_val):
        if median_val < ESSENTIALITY_THRESHOLD:
            return ("#67A9CF", "rgba(103,169,207,0.24)") if is_dark else ("#2166AC", "rgba(33,102,172,0.20)")
        if median_val < 0:
            return ("#92C5DE", "rgba(146,197,222,0.16)") if is_dark else ("#4E84B4", "rgba(78,132,180,0.16)")
        return ("#8E95A5", "rgba(142,149,165,0.14)") if is_dark else ("#8491B4", "rgba(132,145,180,0.14)")

    if compare_mode == "by_gene":
        fig = go.Figure()
        summary = (
            lineage_data[lineage_data.gene.isin(genes)]
            .groupby("gene").crispr_score.agg(["median", "count"]).reset_index()
        )
        if sort_mode == "median":
            summary = summary.sort_values(["median", "gene"])
        else:
            order_map = {gene: idx for idx, gene in enumerate(genes)}
            summary = summary.sort_values("gene", key=lambda s: s.map(order_map))
        labels = []
        for record in summary.itertuples(index=False):
            label = f"{record.gene} (n={record.count})"
            labels.append(label)
            points = lineage_data[lineage_data.gene == record.gene]
            stroke, fill = box_palette(record.median)
            fig.add_trace(go.Box(
                x=points.crispr_score, y=[label] * len(points), orientation="h",
                name=record.gene, showlegend=False,
                marker=dict(color=stroke, size=4, opacity=0.55),
                line=dict(color=stroke, width=1.25),
                fillcolor=fill,
                boxpoints="all" if show_points else "outliers", jitter=0.35, pointpos=0,
                customdata=points[["cell_line", "cell_line_name", "lineage"]].values,
                hovertemplate=(
                    _ui(lang, "Cell line", "细胞系") + ": %{customdata[1]} (%{customdata[0]})<br>"
                    + _ui(lang, "Cancer type", "癌种") + ": %{customdata[2]}<br>"
                    + _ui(lang, "Score", "分数") + ": %{x:.4f}<extra>%{y}</extra>"
                ),
            ))
        fig.update_yaxes(
            categoryorder="array", categoryarray=labels, autorange="reversed",
            automargin=True, showgrid=False,
        )
        fig.update_xaxes(
            title_text=_ui(lang, "CRISPR score", "CRISPR 分数"),
            range=[score_min - pad, score_max + pad],
            showgrid=True, gridcolor=th["plot_grid"],
        )
        fig.add_vline(x=0, line=dict(color=th["plot_reference"], width=0.7))
        fig.add_vline(x=ESSENTIALITY_THRESHOLD, line=dict(dash="dash", color=PLOT_COLORS["threshold"], width=1))
        fig.update_layout(
            height=max(320, len(labels) * 46 + 110),
            margin=dict(l=120, r=25, t=36, b=55), showlegend=False,
        )
        return apply_theme_to_fig(fig, theme=th)

    n_cols = 2 if columns == 2 and n_genes >= 2 else 1
    n_rows = (n_genes + n_cols - 1) // n_cols
    v_space = min(0.08, 0.35 / n_rows) if n_rows > 1 else 0.08
    fig = make_subplots(
        rows=n_rows, cols=n_cols, vertical_spacing=v_space,
        horizontal_spacing=0.18 if n_cols == 2 else 0.05,
        subplot_titles=genes,
    )
    group_max = 1
    for idx, gene in enumerate(genes):
        row = idx // n_cols + 1
        col = idx % n_cols + 1
        sub = lineage_data[lineage_data.gene == gene]
        summary = sub.groupby("lineage").crispr_score.agg(["median", "count"]).reset_index()
        summary = summary.sort_values(["median", "lineage"] if sort_mode == "median" else ["lineage"])
        group_max = max(group_max, len(summary))
        labels = []
        for record in summary.itertuples(index=False):
            label = f"{record.lineage} (n={record.count})"
            labels.append(label)
            points = sub[sub.lineage == record.lineage]
            stroke, fill = box_palette(record.median)
            fig.add_trace(go.Box(
                x=points.crispr_score, y=[label] * len(points), orientation="h",
                name=record.lineage, showlegend=False,
                marker=dict(color=stroke, size=4, opacity=0.5),
                line=dict(color=stroke, width=1.2),
                fillcolor=fill,
                boxpoints="all" if show_points else "outliers", jitter=0.35, pointpos=0,
                customdata=points[["cell_line", "cell_line_name"]].values,
                hovertemplate=(
                    _ui(lang, "Cell line", "细胞系") + ": %{customdata[1]} (%{customdata[0]})<br>"
                    + _ui(lang, "Score", "分数") + ": %{x:.4f}<extra>%{y}</extra>"
                ),
            ), row=row, col=col)
        fig.update_yaxes(
            categoryorder="array", categoryarray=labels, autorange="reversed",
            automargin=True, showgrid=False, row=row, col=col,
        )
        fig.update_xaxes(
            title_text=_ui(lang, "CRISPR score", "CRISPR 分数"),
            range=[score_min - pad, score_max + pad],
            showgrid=True, gridcolor=th["plot_grid"], row=row, col=col,
        )
        fig.add_vline(x=0, line=dict(color=th["plot_reference"], width=0.7), row=row, col=col)
        fig.add_vline(
            x=ESSENTIALITY_THRESHOLD,
            line=dict(dash="dash", color=PLOT_COLORS["threshold"], width=1),
            row=row, col=col,
        )
    fig.update_annotations(font=dict(size=15, family=FONT_FAMILY, color=th["plot_text"]))
    fig.update_layout(
        height=max(300, group_max * 30 + 90) * n_rows + 50,
        margin=dict(l=145, r=25, t=45, b=55), showlegend=False,
    )
    return apply_theme_to_fig(fig, theme=th)


def _add_preference_lollipop(
    fig: go.Figure,
    gene_rows: pd.DataFrame,
    gene: str,
    row: int,
    highlight_set=None,
    theme: str | dict = "light",
    lang: str = "en",
):
    th = get_plot_theme(theme)
    is_dark = _is_dark_theme(theme)
    highlight_set = set(highlight_set or [])
    part = gene_rows.sort_values(["difference", "lineage"], ascending=[False, True])
    order = part["lineage"].tolist()
    diffs = part["difference"].to_numpy(dtype=float)
    limit = float(np.nanmax(np.abs(diffs))) or 0.1
    dep_col = "#67A9CF" if is_dark else "#2166AC"
    res_col = "#E06B64" if is_dark else "#B2182B"
    stems = {dep_col: ([], []), res_col: ([], [])}
    for lineage, diff in zip(order, diffs):
        xs, ys = stems[dep_col if diff < 0 else res_col]
        xs.extend((0, diff, None))
        ys.extend((lineage, lineage, None))
    for color, (xs, ys) in stems.items():
        if not xs:
            continue
        fig.add_trace(go.Scatter(
            x=xs, y=ys, mode="lines",
            line=dict(color=color, width=1.8),
            hoverinfo="skip", showlegend=False,
        ), row=row, col=1)
    has_extra = "median_difference" in part.columns and "dependent_fraction" in part.columns
    if has_extra:
        customdata = np.column_stack((
            part["mean"], part["panel_mean"], part["n_valid"],
            part["median_difference"], part["dependent_fraction"] * 100.0,
        ))
        hovertemplate = (
            f"<b>%{{y}}</b> · {gene}<br>"
            + _ui(lang, "Mean difference", "均值差") + ": %{x:+.3f}<br>"
            + _ui(lang, "Median difference", "中位数差") + ": %{customdata[3]:+.3f}<br>"
            + _ui(lang, "Cancer-type mean", "该癌种均值") + ": %{customdata[0]:.3f}<br>"
            + _ui(lang, "All-lineage mean", "全队列均值") + ": %{customdata[1]:.3f}<br>"
            + _ui(lang, f"Cell lines < {ESSENTIALITY_THRESHOLD}", f"分数 < {ESSENTIALITY_THRESHOLD} 占比")
            + ": %{customdata[4]:.1f}%<br>"
            + _ui(lang, "Valid cell lines", "有效细胞系数") + ": %{customdata[2]:.0f}<extra></extra>"
        )
    else:
        customdata = np.column_stack((part["mean"], part["panel_mean"], part["n_valid"]))
        hovertemplate = (
            f"<b>%{{y}}</b> · {gene}<br>"
            + _ui(lang, "Difference", "差值") + ": %{x:+.3f}<br>"
            + _ui(lang, "Cancer-type mean", "该癌种均值") + ": %{customdata[0]:.3f}<br>"
            + _ui(lang, "All-lineage mean", "全队列均值") + ": %{customdata[1]:.3f}<br>"
            + _ui(lang, "Valid cell lines", "有效细胞系数") + ": %{customdata[2]:.0f}<extra></extra>"
        )
    ring_widths = [2.4 if lin in highlight_set else 0 for lin in order]
    ring_colors = [th["text"] if lin in highlight_set else "rgba(0,0,0,0)" for lin in order]
    fig.add_trace(go.Scatter(
        x=diffs, y=order, mode="markers+text",
        text=[f"{value:+.2f}" for value in diffs],
        textposition=["middle left" if value < 0 else "middle right" for value in diffs],
        textfont=dict(size=12, color=th["plot_text"]),
        cliponaxis=False,
        marker=dict(
            size=[13 if lin in highlight_set else 11 for lin in order],
            color=diffs, colorscale=preference_colorscale(theme),
            cmid=0, cmin=-limit, cmax=limit, showscale=False,
            line=dict(width=ring_widths, color=ring_colors),
        ),
        customdata=customdata,
        hovertemplate=hovertemplate,
        showlegend=False,
    ), row=row, col=1)
    fig.add_vline(x=0, line_width=1.4, line_color=th["plot_axis"], row=row, col=1)
    quiet_grid = "rgba(255,255,255,0.08)" if is_dark else "rgba(29,29,31,0.06)"
    fig.update_yaxes(
        automargin=True, categoryorder="array", categoryarray=order,
        showgrid=True, gridcolor=quiet_grid, gridwidth=1, row=row, col=1,
    )
    fig.update_xaxes(
        range=[-limit * 1.48, limit * 1.58], zeroline=False, showgrid=False,
        title_text=_ui(lang, "← more dependent  ·  less dependent →", "← 更依赖  ·  更不依赖 →"),
        row=row, col=1,
    )


def create_preference_plot(
    table: pd.DataFrame,
    highlight_lineages=None,
    top_k: int | None = None,
    theme: str | dict = "light",
    lang: str = "en",
) -> go.Figure | None:
    """One labeled lollipop per gene. Left of zero is more dependent.

    Each gene keeps its own scale, so a small shift is still readable next to
    a gene whose means move much further. The printed number is the difference.
    """
    shown = table[table["n_valid"] >= PREFERENCE_MIN_N]
    if shown.empty:
        return None
    genes = list(dict.fromkeys(shown["gene"]))
    highlight_set = set(highlight_lineages or [])
    slices = []
    for gene in genes:
        part = shown[shown["gene"] == gene]
        if top_k is not None and len(part) > top_k * 2:
            ordered = part.sort_values(["difference", "lineage"], ascending=[True, True])
            keep_idx = set(ordered.head(top_k).index) | set(ordered.tail(top_k).index)
            if highlight_set:
                keep_idx |= set(part[part["lineage"].isin(highlight_set)].index)
            part = part.loc[part.index.isin(keep_idx)]
        slices.append((gene, part))
    longest = max(len(part) for _, part in slices)
    titles = []
    for gene, part in slices:
        mean = float(part["panel_mean"].iloc[0])
        titles.append(_ui(lang, f"{gene} · all-type mean {mean:.2f}", f"{gene} · 总平均 {mean:.2f}"))
    spacing = 0.12 if len(genes) < 3 else min(0.08, 0.85 / (len(genes) - 1))
    fig = make_subplots(
        rows=len(genes), cols=1, shared_xaxes=False, subplot_titles=titles,
        vertical_spacing=spacing,
    )
    th = get_plot_theme(theme)
    for index, (gene, part) in enumerate(slices, start=1):
        _add_preference_lollipop(fig, part, gene, index, highlight_set=highlight_set, theme=th, lang=lang)
    fig.update_layout(
        height=max(420, len(genes) * (26 * longest + 108)),
        showlegend=False,
        margin=dict(l=8, r=36, t=42, b=56),
    )
    themed = apply_theme_to_fig(fig, theme=th)
    themed.update_annotations(font=dict(size=15, family=FONT_FAMILY, color=th["plot_text"]))
    return themed


def create_preference_heatmap(
    table: pd.DataFrame,
    highlight_lineages=None,
    theme: str | dict = "light",
    lang: str = "en",
) -> go.Figure | None:
    """Compact diverging gene × cancer-type heatmap for comparing many genes at once."""
    shown = table[table["n_valid"] >= PREFERENCE_MIN_N]
    if shown.empty:
        return None
    genes = list(dict.fromkeys(shown["gene"]))
    lineage_order = (
        shown.groupby("lineage")["difference"].mean()
        .sort_values(ascending=True).index.tolist()
    )
    highlight_set = set(highlight_lineages or [])
    y_labels = [f"★ {lin}" if lin in highlight_set else lin for lin in lineage_order]
    pivot = shown.pivot(index="lineage", columns="gene", values="difference").reindex(index=lineage_order, columns=genes)
    z = pivot.to_numpy(dtype=float)
    limit = float(np.nanmax(np.abs(z))) if np.isfinite(z).any() else 0.1
    limit = max(limit, 0.1)
    mean_p = shown.pivot(index="lineage", columns="gene", values="mean").reindex(index=lineage_order, columns=genes)
    panel_p = shown.pivot(index="lineage", columns="gene", values="panel_mean").reindex(index=lineage_order, columns=genes)
    n_p = shown.pivot(index="lineage", columns="gene", values="n_valid").reindex(index=lineage_order, columns=genes)
    dep_p = (
        shown.pivot(index="lineage", columns="gene", values="dependent_fraction").reindex(index=lineage_order, columns=genes) * 100.0
        if "dependent_fraction" in shown.columns else None
    )
    if dep_p is not None:
        customdata = np.stack(
            (mean_p.to_numpy(dtype=float), panel_p.to_numpy(dtype=float),
             n_p.to_numpy(dtype=float), dep_p.to_numpy(dtype=float)),
            axis=-1,
        )
        hover = (
            "<b>%{y}</b> · %{x}<br>"
            + _ui(lang, "Mean difference", "均值差") + ": %{z:+.3f}<br>"
            + _ui(lang, "Cancer-type mean", "该癌种均值") + ": %{customdata[0]:.3f}<br>"
            + _ui(lang, "All-lineage mean", "全队列均值") + ": %{customdata[1]:.3f}<br>"
            + _ui(lang, f"Cell lines < {ESSENTIALITY_THRESHOLD}", f"分数 < {ESSENTIALITY_THRESHOLD} 占比")
            + ": %{customdata[3]:.1f}%<br>"
            + _ui(lang, "Valid cell lines", "有效细胞系数") + ": %{customdata[2]:.0f}<extra></extra>"
        )
    else:
        customdata = np.stack(
            (mean_p.to_numpy(dtype=float), panel_p.to_numpy(dtype=float), n_p.to_numpy(dtype=float)),
            axis=-1,
        )
        hover = (
            "<b>%{y}</b> · %{x}<br>"
            + _ui(lang, "Difference", "差值") + ": %{z:+.3f}<br>"
            + _ui(lang, "Cancer-type mean", "该癌种均值") + ": %{customdata[0]:.3f}<br>"
            + _ui(lang, "All-lineage mean", "全队列均值") + ": %{customdata[1]:.3f}<br>"
            + _ui(lang, "Valid cell lines", "有效细胞系数") + ": %{customdata[2]:.0f}<extra></extra>"
        )
    show_text = len(genes) <= 12
    text_matrix = [[f"{val:+.2f}" if np.isfinite(val) else "" for val in row] for row in z] if show_text else None
    fig = go.Figure(go.Heatmap(
        z=z, x=genes, y=y_labels,
        colorscale=preference_colorscale(theme), zmid=0, zmin=-limit, zmax=limit,
        xgap=2, ygap=2,
        text=text_matrix, texttemplate="%{text}" if show_text else None,
        textfont=dict(size=11),
        customdata=customdata, hovertemplate=hover,
        colorbar=dict(title=_ui(lang, "Mean diff", "均值差"), thickness=12, outlinewidth=0),
    ))
    fig.update_layout(
        height=max(420, len(lineage_order) * 24 + 120),
        margin=dict(l=20, r=24, t=54, b=24),
    )
    themed = apply_theme_to_fig(fig, theme=theme)
    themed.update_yaxes(autorange="reversed", automargin=True, showline=False, ticks="")
    themed.update_xaxes(side="top", automargin=True, showline=False, ticks="")
    return themed


def create_correlation_heatmap(
    corr: pd.DataFrame,
    counts: pd.DataFrame,
    method: str = "pearson",
    theme: str | dict = "light",
    lang: str = "en",
) -> go.Figure:
    """Pairwise gene co-dependency correlation heatmap with gallery tile spacing."""
    genes = list(corr.columns)
    z = corr.to_numpy(dtype=float)
    n_mat = counts.to_numpy(dtype=int)
    show_text = len(genes) <= 12
    text_matrix = [[f"{val:.2f}" if np.isfinite(val) else "" for val in row] for row in z] if show_text else None
    metric_label = "Pearson r" if method == "pearson" else "Spearman ρ"
    colorscale = preference_colorscale(theme)
    cell = 56 if len(genes) <= 8 else (44 if len(genes) <= 14 else 30)
    height = max(360, len(genes) * cell + 115)
    width = min(940, max(480, len(genes) * cell + 215))
    fig = go.Figure(go.Heatmap(
        z=z, x=genes, y=genes,
        colorscale=colorscale,
        zmin=-1.0, zmax=1.0, zmid=0.0,
        xgap=2.5, ygap=2.5,
        text=text_matrix, texttemplate="%{text}" if show_text else None,
        textfont=dict(size=11.5),
        customdata=n_mat,
        hovertemplate=(
            "<b>%{x} × %{y}</b><br>"
            + metric_label + ": %{z:+.3f}<br>"
            + _ui(lang, "Shared cell lines", "共有有效细胞系") + ": %{customdata}<extra></extra>"
        ),
        colorbar=dict(title=metric_label, thickness=12, outlinewidth=0),
    ))
    fig.update_layout(
        width=width, height=height,
        margin=dict(l=24, r=48, t=58, b=24),
    )
    themed = apply_theme_to_fig(fig, theme=theme)
    themed.update_yaxes(autorange="reversed", automargin=True, showline=False, ticks="", scaleanchor="x", scaleratio=1)
    themed.update_xaxes(side="top", automargin=True, showline=False, ticks="", constrain="domain")
    return themed
