#!/usr/bin/env python3
"""Reproduce preprint figures for CRISPR Score Analyzer.

Run from the repository root:
    python preprint/generate_figures.py

Composes Figure 1 (live Light and Charcoal Darkroom browser interface panels when
ui_light.png and ui_dark.png are present) and regenerates all quantitative panels
for Figure 2 and Figure 3 directly from data/depmap_chronos.npz using the exact
analysis and CNS (NPG + ColorBrewer RdBu) plotting helpers from analysis_core.py,
codependency.py, app.py, and figure_export.py.
"""

from __future__ import annotations

import ast
import sys
from html import escape
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from PIL import Image, ImageDraw, ImageFont
from plotly.subplots import make_subplots

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis_core import (  # noqa: E402
    analyze_dataset,
    gene_correlation_matrix,
    get_lineage_frame,
    lineage_preference,
    read_builtin_matrix,
)
from codependency import (  # noqa: E402
    create_codependency_scan_plot,
    create_gene_pair_scatter,
    gene_pair_scatter_frame,
    genomewide_codependency,
)
from figure_export import render_figure_bytes  # noqa: E402

OUT = ROOT / "preprint" / "figures"
OUT.mkdir(parents=True, exist_ok=True)


def load_app_plot_helpers() -> dict:
    """Load plotting helpers directly from app.py without launching Streamlit."""
    names = {
        "get_theme", "ui", "apply_theme_to_fig", "_rank_label_annotations",
        "_ranking_fingerprint", "_cached_background_trace",
        "build_rank_figure", "create_rank_plot", "create_multilayer_rank_plot",
        "_preference_scale", "_add_preference_lollipop", "create_preference_plot",
        "create_correlation_heatmap",
    }
    constants = {"THEMES", "PLOT_COLORS", "FONT_FAMILY", "ESSENTIALITY_THRESHOLD", "PREFERENCE_MIN_N"}
    tree = ast.parse((ROOT / "app.py").read_text(encoding="utf-8"))
    nodes = [
        node for node in tree.body
        if (isinstance(node, ast.FunctionDef) and node.name in names)
        or (
            isinstance(node, ast.Assign)
            and any(isinstance(t, ast.Name) and t.id in constants for t in node.targets)
        )
    ]
    ns = {
        "go": go,
        "np": np,
        "pd": pd,
        "escape": escape,
        "make_subplots": make_subplots,
        "st": SimpleNamespace(session_state={"lang": "en", "theme": "light"}),
    }
    exec(compile(ast.Module(body=nodes, type_ignores=[]), "app.py", "exec"), ns)
    return ns


def save_panel(fig: go.Figure, stem: str, width: int, height: int) -> Path:
    png_bytes = render_figure_bytes(fig, "png", width=width, height=height, scale=3)
    png_path = OUT / f"{stem}.png"
    png_path.write_bytes(png_bytes)
    print(f"Saved {png_path.name} ({len(png_bytes) // 1024} KB)")
    return png_path


def load_label_font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    for path in (
        "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    ):
        if Path(path).exists():
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def compose_2x2(
    panels: list[tuple[str, Path]],
    out_path: Path,
    canvas_w: int = 2400,
    pad_x: int = 44,
    pad_y: int = 40,
    gap_x: int = 36,
    gap_y: int = 36,
) -> None:
    """Arrange four panel PNGs into a crisp 2x2 publication figure with bold A/B/C/D labels."""
    col_w = (canvas_w - 2 * pad_x - gap_x) // 2
    imgs = []
    for label, path in panels:
        im = Image.open(path).convert("RGB")
        scale = col_w / im.width
        resized = im.resize((col_w, int(round(im.height * scale))), Image.Resampling.LANCZOS)
        imgs.append((label, resized))

    row1_h = max(imgs[0][1].height, imgs[1][1].height)
    row2_h = max(imgs[2][1].height, imgs[3][1].height)
    label_band = 48
    canvas_h = pad_y + label_band + row1_h + gap_y + label_band + row2_h + pad_y

    canvas = Image.new("RGB", (canvas_w, canvas_h), (255, 255, 255))
    draw = ImageDraw.Draw(canvas)
    font = load_label_font(48)

    coords = [
        (pad_x, pad_y),
        (pad_x + col_w + gap_x, pad_y),
        (pad_x, pad_y + label_band + row1_h + gap_y),
        (pad_x + col_w + gap_x, pad_y + label_band + row1_h + gap_y),
    ]
    row_heights = [row1_h, row1_h, row2_h, row2_h]

    for (label, im), (x, y), rh in zip(imgs, coords, row_heights):
        draw.text((x + 6, y), label, fill=(29, 29, 31), font=font)
        # Center shorter panel vertically within its row
        y_offset = y + label_band + (rh - im.height) // 2
        canvas.paste(im, (x, y_offset))

    canvas.save(out_path, format="PNG", dpi=(300, 300), optimize=True)
    print(f"Composed {out_path.name} ({canvas.size[0]}x{canvas.size[1]} px, {out_path.stat().st_size // 1024} KB)")


def compose_ui_overview(
    panels: list[tuple[str, str, Path]],
    out_path: Path,
    canvas_w: int = 2400,
    pad_x: int = 38,
    pad_y: int = 34,
    gap_x: int = 36,
) -> None:
    """Arrange two live browser UI captures side-by-side with crisp borders and A/B titles."""
    col_w = (canvas_w - 2 * pad_x - gap_x) // 2
    imgs = []
    for label, subtitle, path in panels:
        im = Image.open(path).convert("RGB")
        scale = col_w / im.width
        resized = im.resize((col_w, int(round(im.height * scale))), Image.Resampling.LANCZOS)
        imgs.append((label, subtitle, resized))

    row_h = max(im.height for _, _, im in imgs)
    label_band = 54
    canvas_h = pad_y + label_band + row_h + pad_y
    canvas = Image.new("RGB", (canvas_w, canvas_h), (255, 255, 255))
    draw = ImageDraw.Draw(canvas)
    font_label = load_label_font(44)
    font_sub = load_label_font(26)

    for idx, (label, subtitle, im) in enumerate(imgs):
        x = pad_x + idx * (col_w + gap_x)
        y = pad_y
        draw.text((x + 4, y), label, fill=(29, 29, 31), font=font_label)
        draw.text((x + 48, y + 12), subtitle, fill=(72, 72, 78), font=font_sub)
        img_y = y + label_band
        canvas.paste(im, (x, img_y))
        # Crisp 2px neutral frame around each viewport capture
        draw.rectangle(
            [x, img_y, x + im.width - 1, img_y + im.height - 1],
            outline=(210, 212, 218),
            width=2,
        )

    canvas.save(out_path, format="PNG", dpi=(300, 300), optimize=True)
    print(f"Composed {out_path.name} ({canvas.size[0]}x{canvas.size[1]} px, {out_path.stat().st_size // 1024} KB)")


def annotate_partial_n(scan_df: pd.DataFrame, cohort_n: int) -> pd.DataFrame:
    """Append (n=...) to gene labels whenever n_shared < cohort_n for transparency."""
    out = scan_df.copy()
    mask = out["n_shared"] < cohort_n
    out.loc[mask, "gene"] = out.loc[mask].apply(
        lambda r: f"{r['gene']} (n={int(r['n_shared'])})", axis=1
    )
    return out


def boost_panel_fonts(
    fig: go.Figure,
    *,
    title_size: int = 15,
    axis_title_size: int = 15,
    x_tick_size: int = 13,
    y_tick_size: int = 13,
    legend_size: int = 12,
    annot_size: float | None = 13,
    trace_text_size: float | None = None,
) -> go.Figure:
    """Apply print-legible typography sizes for multi-panel manuscript figures."""
    fig.update_layout(
        title_font=dict(size=title_size),
        legend=dict(font=dict(size=legend_size)),
    )
    fig.update_xaxes(
        title_font=dict(size=axis_title_size),
        tickfont=dict(size=x_tick_size),
    )
    fig.update_yaxes(
        title_font=dict(size=axis_title_size),
        tickfont=dict(size=y_tick_size),
    )
    if annot_size is not None and fig.layout.annotations:
        for ann in fig.layout.annotations:
            ann.font.size = annot_size
    if trace_text_size is not None:
        for tr in fig.data:
            if hasattr(tr, "textfont") and tr.text is not None:
                tr.textfont = dict(size=trace_text_size, color=getattr(tr.textfont, "color", None))
    return fig


def main() -> None:
    ns = load_app_plot_helpers()
    th = ns["get_theme"]()

    # ------------------------------------------------------------------
    # Figure 1: Live browser interface overview (Light & Charcoal Darkroom)
    # ------------------------------------------------------------------
    ui_light = OUT / "ui_light.png"
    ui_dark = OUT / "ui_dark.png"
    if ui_light.is_file() and ui_dark.is_file():
        compose_ui_overview(
            [
                ("A", "Light mode (#f5f5f7) · four-module dock & NPG ranking curve", ui_light),
                ("B", "Charcoal dark mode (#0f0f11) · identical view in darkroom palette", ui_dark),
            ],
            OUT / "Figure1.png",
        )

    df = read_builtin_matrix(ROOT / "data" / "depmap_chronos.npz")
    rankings, _ = analyze_dataset(df)
    n_lines = len(df)

    # ------------------------------------------------------------------
    # Figure 2A: Pan-cancer default five-gene ranking (NPG palette)
    # ------------------------------------------------------------------
    fig2a = ns["create_rank_plot"](
        rankings,
        ["KRAS", "SOX10", "BRAF", "ESR1", "TP53"],
        essential_gene="MYC",
        nonessential_gene="PTEN",
        n_cell_lines=n_lines,
        show_labels=True,
        point_size=4,
        width=620,
        height=450,
    )
    fig2a.update_yaxes(title_text="Mean Chronos gene effect")
    boost_panel_fonts(fig2a, title_size=14, axis_title_size=15, x_tick_size=13, y_tick_size=13, legend_size=11.5, annot_size=13)
    fig2a.update_layout(margin=dict(l=64, r=22, t=42, b=106))
    p2a = save_panel(fig2a, "fig2A_rank", 620, 450)

    # ------------------------------------------------------------------
    # Figure 2B: Two-layer gene-set annotation (NPG slate + deep indigo)
    # ------------------------------------------------------------------
    fig2b = ns["create_multilayer_rank_plot"](
        rankings,
        ["CDK1", "CDK2", "CDK4", "CDK6", "CCNB1", "CCND1", "CCNE1"],
        ["PLK1", "AURKA"],
        bg_color="#8E99A2",
        hl_color="#3C5488",
        essential_gene="MYC",
        nonessential_gene="PTEN",
        n_cell_lines=n_lines,
        show_labels=True,
        point_size=4,
        width=620,
        height=450,
    )
    fig2b.update_yaxes(title_text="Mean Chronos gene effect")
    for tr in fig2b.data:
        if tr.name == "Background gene set":
            tr.name = "Background: CDK1/2/4/6, CCNB1/D1/E1"
        elif tr.name == "Highlight genes":
            tr.name = "Highlight: PLK1, AURKA"
    boost_panel_fonts(fig2b, title_size=14, axis_title_size=15, x_tick_size=13, y_tick_size=13, legend_size=11.5, annot_size=13)
    fig2b.update_layout(margin=dict(l=64, r=22, t=42, b=106))
    p2b = save_panel(fig2b, "fig2B_multilayer", 620, 450)

    # ------------------------------------------------------------------
    # Figure 2C & 2D: Lineage preference lollipops WITH sample size (n=...)
    # ------------------------------------------------------------------
    lin_df = get_lineage_frame(df, ["SOX10", "ESR1"], rankings, "lineage")
    pref_df = lineage_preference(lin_df, threshold=-0.5)
    pref_df = pref_df.copy()
    pref_df["lineage"] = pref_df.apply(
        lambda r: f"{r['lineage']} (n={int(r['n_valid'])})", axis=1
    )

    sox_pref = pref_df[pref_df["gene"] == "SOX10"].copy()
    fig2c = ns["create_preference_plot"](sox_pref, top_k=4)
    fig2c.update_xaxes(title_text="← more dependent · Δ mean Chronos · less dependent →")
    boost_panel_fonts(fig2c, title_size=15, axis_title_size=13.5, x_tick_size=13, y_tick_size=12.5, annot_size=14.5, trace_text_size=12.5)
    fig2c.update_layout(width=620, height=365, margin=dict(l=205, r=36, t=46, b=58))
    p2c = save_panel(fig2c, "fig2C_sox10_pref", 620, 365)

    esr_pref = pref_df[pref_df["gene"] == "ESR1"].copy()
    fig2d = ns["create_preference_plot"](esr_pref, top_k=4)
    fig2d.update_xaxes(title_text="← more dependent · Δ mean Chronos · less dependent →")
    boost_panel_fonts(fig2d, title_size=15, axis_title_size=13.5, x_tick_size=13, y_tick_size=12.5, annot_size=14.5, trace_text_size=12.5)
    fig2d.update_layout(width=620, height=365, margin=dict(l=222, r=36, t=46, b=58))
    p2d = save_panel(fig2d, "fig2D_esr1_pref", 620, 365)

    compose_2x2(
        [("A", p2a), ("B", p2b), ("C", p2c), ("D", p2d)],
        OUT / "Figure2.png",
    )

    # ------------------------------------------------------------------
    # Figure 3A: Pan-cancer MDM2 scan (explicitly flagging partial n_shared)
    # ------------------------------------------------------------------
    mdm2_scan = genomewide_codependency(df, "MDM2", rankings, method="pearson")
    mdm2_labeled = annotate_partial_n(mdm2_scan, cohort_n=n_lines)
    fig3a = create_codependency_scan_plot(mdm2_labeled, "MDM2", method="pearson", top_n=7, th=th, lang="en")
    boost_panel_fonts(fig3a, title_size=15, axis_title_size=13.5, x_tick_size=13, y_tick_size=13, annot_size=None, trace_text_size=12.5)
    fig3a.update_layout(width=620, height=470, margin=dict(l=140, r=38, t=50, b=58))
    p3a = save_panel(fig3a, "fig3A_mdm2_scan", 620, 470)

    # ------------------------------------------------------------------
    # Figure 3B: Multi-gene correlation matrix (TP53–MDM2 axis, ColorBrewer RdBu)
    # ------------------------------------------------------------------
    tp53_genes = ["TP53", "CHEK2", "MDM2", "MDM4", "PPM1D", "USP7"]
    corr, counts = gene_correlation_matrix(df, tp53_genes, rankings, method="pearson")
    fig3b = ns["create_correlation_heatmap"](corr, counts, method="pearson")
    boost_panel_fonts(fig3b, title_size=15, axis_title_size=14, x_tick_size=13.5, y_tick_size=13.5, annot_size=None, trace_text_size=13)
    for tr in fig3b.data:
        if getattr(tr, "colorbar", None) is not None:
            tr.colorbar.title.font = dict(size=13.5)
            tr.colorbar.tickfont = dict(size=12.5)
    fig3b.update_layout(width=540, height=470, margin=dict(l=30, r=42, t=58, b=24))
    p3b = save_panel(fig3b, "fig3B_p53_matrix", 540, 470)

    # ------------------------------------------------------------------
    # Figure 3C: Cell-line scatter (MDM2 vs TP53, 12-lineage NPG palette)
    # ------------------------------------------------------------------
    scatter_df, stats = gene_pair_scatter_frame(df, "MDM2", "TP53", rankings, "lineage")
    fig3c = create_gene_pair_scatter(scatter_df, stats, "MDM2", "TP53", th=th, lang="en")
    fig3c.update_xaxes(title_text="MDM2 · Chronos gene effect")
    fig3c.update_yaxes(title_text="TP53 · Chronos gene effect")
    boost_panel_fonts(fig3c, title_size=14, axis_title_size=15, x_tick_size=13, y_tick_size=13, legend_size=11.5, annot_size=None)
    fig3c.update_layout(width=660, height=470, margin=dict(l=68, r=20, t=50, b=58))
    p3c = save_panel(fig3c, "fig3C_mdm2_tp53_scatter", 660, 470)

    # ------------------------------------------------------------------
    # Figure 3D: Cohort-scoped BRAF scan within Skin (n = 75, ColorBrewer RdBu)
    # ------------------------------------------------------------------
    skin_df = df[df["lineage"] == "Skin"]
    skin_ranks, _ = analyze_dataset(skin_df)
    braf_skin = genomewide_codependency(skin_df, "BRAF", skin_ranks, method="pearson")
    braf_labeled = annotate_partial_n(braf_skin, cohort_n=len(skin_df))
    fig3d = create_codependency_scan_plot(braf_labeled, "BRAF", method="pearson", top_n=7, th=th, lang="en")
    fig3d.update_layout(
        title=dict(
            text=f"<b>BRAF</b> · Skin lineage only (n = {len(skin_df)} lines, {len(braf_skin):,} genes)",
            font=dict(size=14.5),
            x=0.02,
            xanchor="left",
        ),
        width=620,
        height=470,
        margin=dict(l=105, r=38, t=50, b=58),
    )
    boost_panel_fonts(fig3d, title_size=14.5, axis_title_size=13.5, x_tick_size=13, y_tick_size=13, annot_size=None, trace_text_size=12.5)
    p3d = save_panel(fig3d, "fig3D_braf_skin_scan", 620, 470)

    compose_2x2(
        [("A", p3a), ("B", p3b), ("C", p3c), ("D", p3d)],
        OUT / "Figure3.png",
    )


if __name__ == "__main__":
    main()
