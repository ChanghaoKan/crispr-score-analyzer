#!/usr/bin/env python3
"""Reproduce preprint figures for CRISPR Score Analyzer.

Run from the repository root:
    python preprint/generate_figures.py

Generates Figure 1 (vector software architecture and analytical workflow
schematic, exported to Figure1.svg and 300-DPI Figure1.png) and regenerates all
quantitative panels for Figure 2 and Figure 3 directly from
data/depmap_chronos.npz using the exact analysis and plotting routines from
crispr_score_analyzer.
"""

from __future__ import annotations

import ast
import shutil
import subprocess
import sys
import tempfile
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


def build_architecture_figure1(svg_path: Path, png_path: Path) -> None:
    """Generate Figure 1: eLife-style vector architecture and workflow schematic."""
    svg = """<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="1180" height="510" viewBox="0 0 1180 510">
  <defs>
    <marker id="arrow-indigo" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="6.5" markerHeight="6.5" orient="auto-start-reverse">
      <path d="M 0 1.2 L 9 5 L 0 8.8 z" fill="#3C5488"/>
    </marker>
    <marker id="arrow-slate" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="6.5" markerHeight="6.5" orient="auto-start-reverse">
      <path d="M 0 1.2 L 9 5 L 0 8.8 z" fill="#4A5568"/>
    </marker>
    <style>
      .stage-letter { font-family: "Liberation Sans", Arial, Helvetica, sans-serif; font-size: 16px; font-weight: 700; fill: #1A1D21; }
      .stage-title { font-family: "Liberation Sans", Arial, Helvetica, sans-serif; font-size: 13.5px; font-weight: 700; fill: #1A1D21; }
      .stage-title-core { font-family: "Liberation Sans", Arial, Helvetica, sans-serif; font-size: 13.5px; font-weight: 700; fill: #3C5488; }
      .stage-sub { font-family: "Liberation Sans", Arial, Helvetica, sans-serif; font-size: 11px; fill: #5A6577; }
      .box-title { font-family: "Liberation Sans", Arial, Helvetica, sans-serif; font-size: 12.5px; font-weight: 700; fill: #1A1D21; }
      .box-body { font-family: "Liberation Sans", Arial, Helvetica, sans-serif; font-size: 11.2px; fill: #374151; }
      .box-note { font-family: "Liberation Sans", Arial, Helvetica, sans-serif; font-size: 10.2px; fill: #5A6577; }
      .mono-tag { font-family: "Liberation Mono", "Courier New", monospace; font-size: 11.2px; font-weight: 700; fill: #3C5488; }
      .sub-divider { font-family: "Liberation Sans", Arial, Helvetica, sans-serif; font-size: 10px; font-weight: 700; fill: #5A6577; letter-spacing: 0.7px; }
    </style>
  </defs>

  <!-- Pure white background -->
  <rect x="0" y="0" width="1180" height="510" fill="#FFFFFF"/>

  <!-- Stage Headers -->
  <text x="20" y="28" class="stage-letter">A</text>
  <text x="41" y="27" class="stage-title">Inputs and harmonization</text>
  <line x1="20" y1="38" x2="305" y2="38" stroke="#CBD5E1" stroke-width="1.1"/>

  <text x="365" y="28" class="stage-letter" style="fill:#3C5488;">B</text>
  <text x="386" y="27" class="stage-title-core">Shared analytical core</text>
  <text x="820" y="27" text-anchor="end" class="stage-sub">NumPy / pandas · no SciPy</text>
  <line x1="365" y1="38" x2="820" y2="38" stroke="#3C5488" stroke-width="1.4"/>

  <text x="880" y="28" class="stage-letter">C</text>
  <text x="901" y="27" class="stage-title">Interfaces and outputs</text>
  <line x1="880" y1="38" x2="1160" y2="38" stroke="#CBD5E1" stroke-width="1.1"/>

  <!-- ==================== COLUMN A: INPUTS ==================== -->
  <!-- Box A1 -->
  <rect x="20" y="54" width="285" height="132" rx="4" fill="#F8F9FB" stroke="#CBD5E1" stroke-width="1.1"/>
  <rect x="20" y="54" width="4" height="132" rx="2" fill="#3C5488"/>
  <text x="36" y="80" class="box-title">Built-in DepMap Public 26Q3</text>
  <text x="36" y="103" class="box-body">Chronos gene effect matrix</text>
  <text x="36" y="123" class="box-body">1,361 cancer cell lines · 30 lineages</text>
  <text x="36" y="143" class="box-body">18,434 genes · float32 binary storage</text>
  <text x="36" y="166" class="box-note">Non-commercial DepMap Terms</text>

  <!-- Box A2 -->
  <rect x="20" y="204" width="285" height="132" rx="4" fill="#F8F9FB" stroke="#CBD5E1" stroke-width="1.1"/>
  <rect x="20" y="204" width="4" height="132" rx="2" fill="#3C5488"/>
  <text x="36" y="230" class="box-title">Custom CRISPR screen matrices</text>
  <text x="36" y="253" class="box-body">User-supplied UTF-8 CSV tables</text>
  <text x="36" y="273" class="box-body">Strict gene-header schema guard</text>
  <text x="36" y="293" class="box-body">Duplicate &amp; missingness diagnostics</text>
  <text x="36" y="316" class="box-note">Optional disease-subtype filtering</text>

  <!-- Box A3 -->
  <rect x="20" y="354" width="285" height="132" rx="4" fill="#F8F9FB" stroke="#CBD5E1" stroke-width="1.1"/>
  <rect x="20" y="354" width="4" height="132" rx="2" fill="#3C5488"/>
  <text x="36" y="380" class="box-title">Offline HGNC harmonization</text>
  <text x="36" y="403" class="box-body">40,402 unambiguous alias mappings</text>
  <text x="36" y="423" class="box-body">Batch lists &amp; DepMap-style headers</text>
  <text x="36" y="443" class="box-body">Explicit alias &amp; unmatched audit log</text>
  <text x="36" y="466" class="box-note">Zero external network fetches</text>

  <!-- Converging Bus Arrows A -> B -->
  <line x1="305" y1="120" x2="335" y2="120" stroke="#3C5488" stroke-width="1.35"/>
  <line x1="305" y1="270" x2="335" y2="270" stroke="#3C5488" stroke-width="1.35"/>
  <line x1="305" y1="420" x2="335" y2="420" stroke="#3C5488" stroke-width="1.35"/>
  <line x1="335" y1="120" x2="335" y2="420" stroke="#3C5488" stroke-width="1.35"/>
  <line x1="335" y1="270" x2="361" y2="270" stroke="#3C5488" stroke-width="1.45" marker-end="url(#arrow-indigo)"/>

  <!-- ==================== COLUMN B: SHARED CORE ==================== -->
  <rect x="365" y="54" width="455" height="432" rx="5" fill="#F3F6FA" stroke="#3C5488" stroke-width="1.25"/>
  <text x="383" y="77" class="mono-tag">crispr_score_analyzer</text>
  <text x="537" y="77" class="stage-sub">·  float64 descriptive estimators</text>

  <!-- Module B1 -->
  <rect x="381" y="92" width="423" height="118" rx="4" fill="#FFFFFF" stroke="#B8C5DC" stroke-width="1.05"/>
  <text x="397" y="116" class="box-title">1. Cohort-scoped essentiality ranking</text>
  <text x="397" y="140" class="box-body">Active cohort mean <tspan font-style="italic">x&#772;</tspan><tspan dy="3" font-size="8.5">g,C</tspan><tspan dy="-3"> and empirical percentile </tspan><tspan font-style="italic">P</tspan><tspan dy="3" font-size="8.5">g,C</tspan></text>
  <text x="397" y="163" class="box-body">Reference anchors (<tspan font-style="italic">MYC</tspan>, <tspan font-style="italic">PTEN</tspan>) and −0.5 visual guide</text>
  <text x="397" y="186" class="box-body">Multi-gene batch curves and two-layer gene-set overlays</text>

  <!-- Module B2 -->
  <rect x="381" y="225" width="423" height="118" rx="4" fill="#FFFFFF" stroke="#B8C5DC" stroke-width="1.05"/>
  <text x="397" y="249" class="box-title">2. Lineage selectivity and subtype skew</text>
  <text x="397" y="273" class="box-body">Joint mean shift &#916;<tspan font-style="italic">&#956;</tspan>, median shift &#916;med, and fraction <tspan font-style="italic">f</tspan> (&lt; −0.5)</text>
  <text x="397" y="296" class="box-body">Distinguishes uniform lineage effects from bimodal subtypes</text>
  <text x="397" y="319" class="box-body">Paired most- and least-dependent cell-line model tables</text>

  <!-- Module B3 -->
  <rect x="381" y="358" width="423" height="118" rx="4" fill="#FFFFFF" stroke="#B8C5DC" stroke-width="1.05"/>
  <text x="397" y="382" class="box-title">3. Pairwise-complete co-dependency engine</text>
  <text x="397" y="406" class="box-body">2,048-column float64 block scans with zero imputation</text>
  <text x="397" y="429" class="box-body">Exact BLAS Pearson <tspan font-style="italic">r</tspan> and average-rank Spearman <tspan font-style="italic">&#961;</tspan> with <tspan font-style="italic">n</tspan><tspan dy="3" font-size="8.5">shared</tspan></text>
  <text x="397" y="452" class="box-body">Multi-gene correlation matrix (2–25 genes) and lineage scatter</text>

  <!-- Branching Bus Arrows B -> C -->
  <line x1="820" y1="270" x2="848" y2="270" stroke="#4A5568" stroke-width="1.35"/>
  <line x1="848" y1="100" x2="848" y2="440" stroke="#4A5568" stroke-width="1.35"/>
  <line x1="848" y1="100" x2="876" y2="100" stroke="#4A5568" stroke-width="1.35" marker-end="url(#arrow-slate)"/>
  <line x1="848" y1="205" x2="876" y2="205" stroke="#4A5568" stroke-width="1.35" marker-end="url(#arrow-slate)"/>
  <line x1="848" y1="335" x2="876" y2="335" stroke="#4A5568" stroke-width="1.35" marker-end="url(#arrow-slate)"/>
  <line x1="848" y1="440" x2="876" y2="440" stroke="#4A5568" stroke-width="1.35" marker-end="url(#arrow-slate)"/>

  <!-- ==================== COLUMN C: INTERFACES & OUTPUTS ==================== -->
  <!-- Box C1 -->
  <rect x="880" y="54" width="280" height="92" rx="4" fill="#FFFFFF" stroke="#CBD5E1" stroke-width="1.1"/>
  <rect x="880" y="54" width="4" height="92" rx="2" fill="#00796B"/>
  <text x="896" y="79" class="box-title">Interactive browser platform</text>
  <text x="896" y="102" class="box-body">Streamlit web application (EN / ZH)</text>
  <text x="896" y="123" class="box-body">Cohort filters, drill-downs &amp; scatter</text>

  <!-- Box C2 -->
  <rect x="880" y="160" width="280" height="92" rx="4" fill="#FFFFFF" stroke="#CBD5E1" stroke-width="1.1"/>
  <rect x="880" y="160" width="4" height="92" rx="2" fill="#00796B"/>
  <text x="896" y="185" class="box-title">Python API and CLI</text>
  <text x="896" y="208" class="box-body">Scriptable multi-omic analysis pipelines</text>
  <text x="896" y="229" class="box-body">Subcommands: info, rank, codep, app</text>

  <!-- Sub-group Divider -->
  <text x="880" y="276" class="sub-divider">REPRODUCIBLE OUTPUTS</text>
  <line x1="1032" y1="272" x2="1160" y2="272" stroke="#E2E8F0" stroke-width="1"/>

  <!-- Box C3 -->
  <rect x="880" y="290" width="280" height="92" rx="4" fill="#F8F9FB" stroke="#CBD5E1" stroke-width="1.1"/>
  <rect x="880" y="290" width="4" height="92" rx="2" fill="#4A5568"/>
  <text x="896" y="315" class="box-title">Publication vector graphics</text>
  <text x="896" y="338" class="box-body">White-background PDF, SVG &amp; PNG</text>
  <text x="896" y="359" class="box-body">Colorblind-safe NPG &amp; RdBu palettes</text>

  <!-- Box C4 -->
  <rect x="880" y="394" width="280" height="92" rx="4" fill="#F8F9FB" stroke="#CBD5E1" stroke-width="1.1"/>
  <rect x="880" y="394" width="4" height="92" rx="2" fill="#4A5568"/>
  <text x="896" y="419" class="box-title">Self-contained provenance ZIP</text>
  <text x="896" y="442" class="box-body">Ranked tables &amp; schema diagnostics</text>
  <text x="896" y="463" class="box-body">Cohort scope, versions &amp; SHA-256 log</text>
</svg>
"""
    svg_path.write_text(svg, encoding="utf-8")
    print(f"Wrote {svg_path.name} ({svg_path.stat().st_size // 1024} KB)")

    chrome = shutil.which("google-chrome") or shutil.which("chromium") or shutil.which("chromium-browser")
    if chrome:
        html_doc = (
            "<!DOCTYPE html><html><body style='margin:0;padding:0;background:#fff;overflow:hidden;'>"
            f"{svg}</body></html>"
        )
        with tempfile.NamedTemporaryFile("w", suffix=".html", encoding="utf-8", delete=False) as tf:
            tf.write(html_doc)
            html_path = Path(tf.name)
        try:
            raw_png = png_path.with_suffix(".raw.png")
            subprocess.run(
                [
                    chrome,
                    "--headless=new",
                    "--no-sandbox",
                    "--disable-gpu",
                    "--hide-scrollbars",
                    "--window-size=1180,510",
                    "--force-device-scale-factor=2",
                    f"--screenshot={raw_png}",
                    html_path.as_uri(),
                ],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            im = Image.open(raw_png).convert("RGB")
            im = im.crop((0, 0, min(im.width, 2360), min(im.height, 1020)))
            im.save(png_path, format="PNG", dpi=(300, 300), optimize=True)
            raw_png.unlink(missing_ok=True)
            print(f"Rendered {png_path.name} ({im.size[0]}x{im.size[1]} px, {png_path.stat().st_size // 1024} KB)")
        finally:
            html_path.unlink(missing_ok=True)
    else:
        print(f"Chrome not found; kept vector {svg_path.name} without PNG rasterization.")


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
    # Figure 1: Vector software architecture and analytical workflow
    # ------------------------------------------------------------------
    build_architecture_figure1(OUT / "Figure1.svg", OUT / "Figure1.png")

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
