#!/usr/bin/env python3
"""Reproduce Figure 1 and Figure 2 for the CRISPR Score Analyzer bioRxiv preprint.

Run from the repository root:
    python preprint/generate_figures.py

Reads data/depmap_chronos.npz using the exact analysis and plotting helpers from
analysis_core.py, codependency.py, app.py, and figure_export.py, then composes
300-DPI composite PNGs (preprint/figures/Figure1.png and Figure2.png).
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
    font = load_label_font(44)

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


def annotate_partial_n(scan_df: pd.DataFrame, cohort_n: int) -> pd.DataFrame:
    """Append (n=...) to gene labels whenever n_shared < cohort_n for transparency."""
    out = scan_df.copy()
    mask = out["n_shared"] < cohort_n
    out.loc[mask, "gene"] = out.loc[mask].apply(
        lambda r: f"{r['gene']} (n={int(r['n_shared'])})", axis=1
    )
    return out


def main() -> None:
    ns = load_app_plot_helpers()
    th = ns["get_theme"]()

    df = read_builtin_matrix(ROOT / "data" / "depmap_chronos.npz")
    rankings, _ = analyze_dataset(df)
    n_lines = len(df)

    # ------------------------------------------------------------------
    # Figure 1A: Pan-cancer default five-gene ranking
    # ------------------------------------------------------------------
    fig1a = ns["create_rank_plot"](
        rankings,
        ["KRAS", "SOX10", "BRAF", "ESR1", "TP53"],
        essential_gene="MYC",
        nonessential_gene="PTEN",
        n_cell_lines=n_lines,
        show_labels=True,
        point_size=4,
        width=620,
        height=440,
    )
    fig1a.update_yaxes(title_text="Mean Chronos gene effect")
    p1a = save_panel(fig1a, "fig1A_rank", 620, 440)

    # ------------------------------------------------------------------
    # Figure 1B: Two-layer gene-set annotation with explicit background gene list
    # ------------------------------------------------------------------
    fig1b = ns["create_multilayer_rank_plot"](
        rankings,
        ["CDK1", "CDK2", "CDK4", "CDK6", "CCNB1", "CCND1", "CCNE1"],
        ["PLK1", "AURKA"],
        bg_color="#8E8E93",
        hl_color="#0071E3",
        essential_gene="MYC",
        nonessential_gene="PTEN",
        n_cell_lines=n_lines,
        show_labels=True,
        point_size=4,
        width=620,
        height=440,
    )
    fig1b.update_yaxes(title_text="Mean Chronos gene effect")
    for tr in fig1b.data:
        if tr.name == "Background gene set":
            tr.name = "Background: CDK1/2/4/6, CCNB1/D1/E1"
        elif tr.name == "Highlight genes":
            tr.name = "Highlight: PLK1, AURKA"
    p1b = save_panel(fig1b, "fig1B_multilayer", 620, 440)

    # ------------------------------------------------------------------
    # Figure 1C & 1D: Lineage preference lollipops WITH sample size (n=...)
    # ------------------------------------------------------------------
    lin_df = get_lineage_frame(df, ["SOX10", "ESR1"], rankings, "lineage")
    pref_df = lineage_preference(lin_df, threshold=-0.5)
    pref_df = pref_df.copy()
    pref_df["lineage"] = pref_df.apply(
        lambda r: f"{r['lineage']} (n={int(r['n_valid'])})", axis=1
    )

    sox_pref = pref_df[pref_df["gene"] == "SOX10"].copy()
    fig1c = ns["create_preference_plot"](sox_pref, top_k=4)
    fig1c.update_xaxes(title_text="← more dependent · Δ mean Chronos · less dependent →")
    fig1c.update_layout(width=620, height=340, margin=dict(l=170, r=32, t=42, b=52))
    p1c = save_panel(fig1c, "fig1C_sox10_pref", 620, 340)

    esr_pref = pref_df[pref_df["gene"] == "ESR1"].copy()
    fig1d = ns["create_preference_plot"](esr_pref, top_k=4)
    fig1d.update_xaxes(title_text="← more dependent · Δ mean Chronos · less dependent →")
    fig1d.update_layout(width=620, height=340, margin=dict(l=185, r=32, t=42, b=52))
    p1d = save_panel(fig1d, "fig1D_esr1_pref", 620, 340)

    compose_2x2(
        [("A", p1a), ("B", p1b), ("C", p1c), ("D", p1d)],
        OUT / "Figure1.png",
    )

    # ------------------------------------------------------------------
    # Figure 2A: Pan-cancer MDM2 scan (explicitly flagging partial n_shared)
    # ------------------------------------------------------------------
    mdm2_scan = genomewide_codependency(df, "MDM2", rankings, method="pearson")
    mdm2_labeled = annotate_partial_n(mdm2_scan, cohort_n=n_lines)
    fig2a = create_codependency_scan_plot(mdm2_labeled, "MDM2", method="pearson", top_n=7, th=th, lang="en")
    fig2a.update_layout(width=620, height=460, margin=dict(l=118, r=32, t=48, b=54))
    p2a = save_panel(fig2a, "fig2A_mdm2_scan", 620, 460)

    # ------------------------------------------------------------------
    # Figure 2B: Multi-gene correlation matrix (TP53–MDM2 axis)
    # ------------------------------------------------------------------
    tp53_genes = ["TP53", "CHEK2", "MDM2", "MDM4", "PPM1D", "USP7"]
    corr, counts = gene_correlation_matrix(df, tp53_genes, rankings, method="pearson")
    fig2b = ns["create_correlation_heatmap"](corr, counts, method="pearson")
    fig2b.update_layout(width=540, height=460, margin=dict(l=20, r=36, t=54, b=20))
    p2b = save_panel(fig2b, "fig2B_p53_matrix", 540, 460)

    # ------------------------------------------------------------------
    # Figure 2C: Cell-line scatter (MDM2 vs TP53) with Chronos axis labels
    # ------------------------------------------------------------------
    scatter_df, stats = gene_pair_scatter_frame(df, "MDM2", "TP53", rankings, "lineage")
    fig2c = create_gene_pair_scatter(scatter_df, stats, "MDM2", "TP53", th=th, lang="en")
    fig2c.update_xaxes(title_text="MDM2 · Chronos gene effect")
    fig2c.update_yaxes(title_text="TP53 · Chronos gene effect")
    fig2c.update_layout(width=660, height=460, margin=dict(l=64, r=20, t=46, b=52))
    p2c = save_panel(fig2c, "fig2C_mdm2_tp53_scatter", 660, 460)

    # ------------------------------------------------------------------
    # Figure 2D: Cohort-scoped BRAF scan within Skin (n = 75)
    # ------------------------------------------------------------------
    skin_df = df[df["lineage"] == "Skin"]
    skin_ranks, _ = analyze_dataset(skin_df)
    braf_skin = genomewide_codependency(skin_df, "BRAF", skin_ranks, method="pearson")
    braf_labeled = annotate_partial_n(braf_skin, cohort_n=len(skin_df))
    fig2d = create_codependency_scan_plot(braf_labeled, "BRAF", method="pearson", top_n=7, th=th, lang="en")
    fig2d.update_layout(
        title=dict(
            text=f"<b>BRAF</b> · Skin lineage only (n = {len(skin_df)} lines, {len(braf_skin):,} genes)",
            font=dict(size=13),
            x=0.02,
            xanchor="left",
        ),
        width=620,
        height=460,
        margin=dict(l=90, r=32, t=48, b=54),
    )
    p2d = save_panel(fig2d, "fig2D_braf_skin_scan", 620, 460)

    compose_2x2(
        [("A", p2a), ("B", p2b), ("C", p2c), ("D", p2d)],
        OUT / "Figure2.png",
    )


if __name__ == "__main__":
    main()
