#!/usr/bin/env python3
"""Build preprint/manuscript.pdf from print-typeset HTML via headless Chromium / Playwright.

Run from the repository root:
    python3 preprint/build_pdf.py
"""

from __future__ import annotations

import base64
import tempfile
from pathlib import Path
from playwright.sync_api import sync_playwright

PREPRINT_DIR = Path(__file__).resolve().parent
FIGURES_DIR = PREPRINT_DIR / "figures"
OUTPUT_PDF = PREPRINT_DIR / "manuscript.pdf"


def img_to_data_uri(path: Path) -> str:
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{data}"


def build_html() -> str:
    fig1_uri = img_to_data_uri(FIGURES_DIR / "Figure1.png")
    fig2_uri = img_to_data_uri(FIGURES_DIR / "Figure2.png")
    fig3_uri = img_to_data_uri(FIGURES_DIR / "Figure3.png")

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<title>CRISPR Score Analyzer — Preprint</title>
<style>
  @page {{
    size: A4 portrait;
    margin: 14.2mm 17mm 14.2mm 17mm;
  }}

  * {{
    box-sizing: border-box;
  }}

  html, body {{
    margin: 0;
    padding: 0;
    background: #ffffff;
    color: #1d1d1f;
    font-family: "Liberation Sans", -apple-system, BlinkMacSystemFont, "Helvetica Neue", Arial, sans-serif;
    font-size: 8.55pt;
    line-height: 1.40;
    -webkit-font-smoothing: antialiased;
    text-rendering: optimizeLegibility;
  }}

  /* Header / Title Block (NPG Deep Indigo #3C5488 editorial accent) */
  .doc-kicker {{
    display: flex;
    justify-content: space-between;
    align-items: baseline;
    font-size: 7.2pt;
    font-weight: 700;
    letter-spacing: 0.085em;
    text-transform: uppercase;
    color: #3C5488;
    padding-bottom: 3pt;
    border-bottom: 1.1pt solid #1d1d1f;
    margin-bottom: 6.5pt;
  }}

  .doc-kicker .kicker-right {{
    color: #6e6e73;
    font-weight: 600;
    letter-spacing: 0.055em;
  }}

  h1.title {{
    font-size: 15.2pt;
    line-height: 1.20;
    font-weight: 700;
    letter-spacing: -0.018em;
    color: #1d1d1f;
    margin: 0 0 4.8pt 0;
  }}

  .author-block {{
    margin: 0 0 6.5pt 0;
  }}

  .author-name {{
    font-size: 9.4pt;
    font-weight: 700;
    color: #1d1d1f;
    margin: 0 0 1.5pt 0;
  }}

  .affiliations {{
    font-size: 7.8pt;
    line-height: 1.38;
    color: #6e6e73;
    margin: 0;
    text-align: left;
  }}

  /* Structured Abstract Card */
  .abstract-card {{
    background: #f5f5f7;
    border-left: 2.4pt solid #3C5488;
    border-radius: 0 4.5pt 4.5pt 0;
    padding: 6.5pt 9.5pt 6.2pt 9.5pt;
    margin: 0 0 7.5pt 0;
  }}

  .abstract-label {{
    font-size: 7.3pt;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 0.085em;
    color: #3C5488;
    margin: 0 0 3pt 0;
  }}

  .abstract-card p {{
    font-size: 7.95pt;
    line-height: 1.36;
    color: #1d1d1f;
    margin: 0 0 3.6pt 0;
    text-align: justify;
    hyphens: auto;
  }}

  .abstract-card p.avail {{
    text-align: left;
    hyphens: none;
  }}

  .abstract-card p:last-child {{
    margin-bottom: 0;
  }}

  .keywords {{
    font-size: 7.6pt !important;
    color: #424245 !important;
    padding-top: 3pt;
    border-top: 0.5pt solid #d2d2d7;
    margin-top: 4pt !important;
    text-align: left !important;
  }}

  /* Sections */
  h2 {{
    font-size: 10.5pt;
    font-weight: 700;
    letter-spacing: -0.01em;
    color: #1d1d1f;
    margin: 9.5pt 0 4pt 0;
    padding-bottom: 2.2pt;
    border-bottom: 0.6pt solid #d2d2d7;
    break-after: avoid;
    page-break-after: avoid;
  }}

  h2 .sec-num {{
    color: #3C5488;
    margin-right: 4pt;
  }}

  h3 {{
    font-size: 9.15pt;
    font-weight: 700;
    color: #1d1d1f;
    margin: 7.5pt 0 3pt 0;
    break-after: avoid;
    page-break-after: avoid;
  }}

  p {{
    margin: 0 0 4.6pt 0;
    color: #2c2c2e;
    text-align: justify;
    hyphens: auto;
    orphans: 3;
    widows: 3;
  }}

  ol, ul {{
    margin: 2.2pt 0 5pt 0;
    padding-left: 14pt;
  }}

  li {{
    margin-bottom: 2.6pt;
    color: #2c2c2e;
    text-align: justify;
    hyphens: auto;
  }}

  li:last-child {{
    margin-bottom: 0;
  }}

  ul.sublist {{
    margin: 2.2pt 0 1.5pt 0;
    padding-left: 12.5pt;
    list-style-type: disc;
  }}

  ul.sublist li {{
    margin-bottom: 2pt;
  }}

  /* Inline formatting */
  .gene {{
    font-style: italic;
  }}

  .nowrap {{
    white-space: nowrap;
  }}

  .keep-together {{
    break-inside: avoid;
    page-break-inside: avoid;
  }}

  code {{
    font-family: "Liberation Mono", "Courier New", monospace;
    font-size: 7.55pt;
    background: #f5f5f7;
    color: #1d1d1f;
    padding: 0.3pt 2.2pt;
    border-radius: 2pt;
    white-space: nowrap;
  }}

  code.wrap-code {{
    white-space: normal;
    overflow-wrap: anywhere;
  }}

  .abstract-card code {{
    background: #e8e8ed;
  }}

  a {{
    color: #2166AC;
    text-decoration: none;
  }}

  strong {{
    color: #1d1d1f;
  }}

  sup {{
    font-size: 68%;
    line-height: 0;
    vertical-align: super;
  }}

  /* Table 1 */
  .table-block {{
    margin: 7pt 0 8pt 0;
    break-inside: avoid;
    page-break-inside: avoid;
  }}

  .table-caption {{
    font-size: 8pt;
    line-height: 1.38;
    color: #1d1d1f;
    margin: 0 0 4pt 0;
    text-align: left;
  }}

  table.cmp-table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 7.35pt;
    line-height: 1.30;
    border-top: 1.1pt solid #1d1d1f;
    border-bottom: 1.1pt solid #1d1d1f;
  }}

  table.cmp-table thead th {{
    background: #f5f5f7;
    color: #1d1d1f;
    font-weight: 700;
    text-align: left;
    padding: 4pt 5pt;
    border-bottom: 0.7pt solid #86868b;
    vertical-align: bottom;
  }}

  table.cmp-table thead th.col-csa {{
    color: #3C5488;
  }}

  table.cmp-table tbody td {{
    padding: 3.3pt 5pt;
    border-bottom: 0.4pt solid #e5e5ea;
    vertical-align: top;
    color: #2c2c2e;
  }}

  table.cmp-table tbody tr:last-child td {{
    border-bottom: none;
  }}

  table.cmp-table tbody tr:nth-child(even) td {{
    background: #fafafa;
  }}

  table.cmp-table td.cap-col {{
    font-weight: 600;
    color: #1d1d1f;
    width: 26%;
  }}

  table.cmp-table td.csa-col {{
    color: #1d1d1f;
    width: 30%;
  }}

  /* Figures */
  .figure-block {{
    margin: 6.5pt 0 7.5pt 0;
    padding: 5.5pt 6.5pt 5.5pt 6.5pt;
    border: 0.5pt solid #e5e5ea;
    border-radius: 4.5pt;
    background: #ffffff;
    break-inside: avoid;
    page-break-inside: avoid;
  }}

  .figure-block img {{
    display: block;
    width: 99%;
    height: auto;
    margin: 0 auto 5pt auto;
  }}

  .figure-block.fig2 img {{
    width: 93%;
  }}

  .figure-block.fig3 img {{
    width: 94%;
  }}

  .figure-caption {{
    font-size: 7.7pt;
    line-height: 1.38;
    color: #424245;
    margin: 0;
    text-align: justify;
    hyphens: auto;
  }}

  .figure-caption strong {{
    color: #1d1d1f;
  }}

  /* Scope callout list */
  .scope-list {{
    background: #f5f5f7;
    border-radius: 4.5pt;
    padding: 5pt 10pt 5pt 19pt;
    margin: 3pt 0 5.5pt 0;
    break-inside: avoid;
    page-break-inside: avoid;
  }}

  .scope-list li {{
    font-size: 7.9pt;
    line-height: 1.35;
    margin-bottom: 2.2pt;
  }}

  .scope-list li:last-child {{
    margin-bottom: 0;
  }}

  /* Disclosures */
  .disclosure-block p {{
    font-size: 7.85pt;
    line-height: 1.34;
    margin-bottom: 2.6pt;
    text-align: left;
    hyphens: none;
  }}

  /* Author-year hanging-indent References */
  ul.references {{
    list-style: none;
    padding-left: 0;
    margin: 3pt 0 0 0;
  }}

  ul.references li {{
    font-size: 7.5pt;
    line-height: 1.32;
    color: #424245;
    padding-left: 13pt;
    text-indent: -13pt;
    margin-bottom: 2.2pt;
    text-align: left;
    hyphens: none;
  }}
</style>
</head>
<body>

  <!-- Header & Title -->
  <div class="doc-kicker">
    <span>Bioinformatics &middot; Application Note</span>
    <span class="kicker-right">OSF Preprints</span>
  </div>

  <h1 class="title">CRISPR Score Analyzer: an interactive, cohort-scoped browser platform for DepMap gene essentiality and co-dependency analysis</h1>

  <div class="author-block">
    <p class="author-name">Changhao Kan<sup>1,*</sup></p>
    <p class="affiliations">
      <sup>1</sup>Institute of Molecular Physiology, Shenzhen Bay Laboratory, Shenzhen, China<br />
      <sup>*</sup>Correspondence: <a href="mailto:kanchanghao@szbl.ac.cn"><code>kanchanghao@szbl.ac.cn</code></a> &nbsp;&middot;&nbsp;
      ORCID: <a href="https://orcid.org/0009-0000-4014-5923"><code>0009-0000-4014-5923</code></a>
    </p>
  </div>

  <!-- Structured Abstract -->
  <div class="abstract-card">
    <div class="abstract-label">Abstract</div>
    <p>
      <strong>Summary:</strong> Genome-scale CRISPR-Cas9 loss-of-function screens in the Cancer Dependency Map (DepMap) are widely used to prioritize candidate cancer targets, evaluate lineage-restricted vulnerabilities, and examine functional gene relationships. Tasks beyond single-gene portal lookups (placing multi-gene candidate panels on a cohort-wide essentiality curve, separating pan-cancer essentiality from lineage-specific dependency, pairing sensitive and resistant cell lines for wet-lab validation, running single-lineage genome-wide co-dependency scans, and exporting vector figures with analysis provenance) still typically require ad-hoc scripting. We present <strong>CRISPR Score Analyzer (v1.2.0)</strong>, an open-source Python/Streamlit browser platform bundling the DepMap Public 26Q3 Chronos gene effect matrix (1,361 cancer cell lines &times; 18,434 genes across 30 Oncotree lineages) and 40,402 unambiguous HGNC alias mappings, with support for validated custom CSV uploads. Across four modules (<em>Gene ranking</em>, <em>Cancer-type comparison</em>, <em>Gene-set annotation</em>, and <em>Co-dependency</em>), users can analyze pan-cancer or lineage-restricted cohorts, run chunked vectorized Pearson or Spearman scans in seconds using NumPy and pandas without SciPy, view colorblind-safe Nature Publishing Group (NPG) and 7-stop ColorBrewer <code>RdBu</code> figures across light and deep-charcoal dark themes, and export white-paper PDF/SVG/300-DPI PNG figures alongside SHA-256 ZIP provenance bundles. Benchmark examples recover lineage-selective dependencies (<span class="gene">SOX10</span> in skin, <span class="gene">ESR1</span> in ER-positive breast models), the inverse correlation between <span class="gene">MDM2</span> and <span class="gene">TP53</span> Chronos gene-effect profiles across cancer cell lines (<span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.702</span>), and higher descriptive <span class="gene">BRAF</span>–MAPK correlations within skin cancer cell lines (<span class="gene">MAPK1</span> <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.798</span>).
    </p>
    <p class="avail">
      <strong>Availability and Implementation:</strong> Source code for v1.2.0 is freely available under the MIT license at
      <a href="https://github.com/ChanghaoKan/crispr-score-analyzer">https://github.com/ChanghaoKan/crispr-score-analyzer</a>
      (software release tag <code>v1.2.0</code>: <a href="https://github.com/ChanghaoKan/crispr-score-analyzer/releases/tag/v1.2.0">https://github.com/ChanghaoKan/crispr-score-analyzer/releases/tag/v1.2.0</a>; commit <code>5c18f9b</code>; Zenodo version 1.2.0 DOI: <a href="https://doi.org/10.5281/zenodo.24458535"><code>10.5281/zenodo.24458535</code></a>; concept DOI: <code>10.5281/zenodo.19607602</code>).
      A hosted web instance runs at
      <a href="https://crispr-score-analyzer-szbl-denglab.streamlit.app/">https://crispr-score-analyzer-szbl-denglab.streamlit.app/</a>.
      MIT source archives exclude the bundled DepMap Public 26Q3 matrix, which is governed separately by Broad Institute DepMap Terms and Conditions (non-commercial research use only; reposted in full in <code>data/DEPMAP_TERMS.md</code>). The updated browser UI, CNS plotting palettes, manuscript figures, and statistical reproduction scripts (<code>preprint/generate_figures.py</code> and <code>preprint/reproduce_stats.py</code>) are pinned at repository tag <code>preprint-v2</code> (with analytical kernels unchanged from <code>v1.2.0</code>). Exact reproduction from a clean Python 3.12 environment requires checking out <code>preprint-v2</code>, verifying <code>data/depmap_chronos.npz</code> against <code>data/depmap_chronos.provenance.json</code> (SHA-256 <code class="wrap-code">ddf274da2aa5f89cf630bc8150dcd3a2c722466bd18ff6fb6d18fe161f8a0c67</code>, or rebuilding via <code>scripts/build_builtin_matrix.py</code>), installing <code>requirements.txt</code>, and running both scripts.
    </p>
    <p class="keywords">
      <strong>Keywords:</strong> DepMap; CRISPR-Cas9 screen; Chronos gene effect; gene essentiality; co-dependency; cancer cell lines; reproducibility
    </p>
  </div>

  <!-- Section 1 -->
  <h2><span class="sec-num">1.</span>Introduction</h2>
  <p>
    Pooled genome-scale CRISPR-Cas9 loss-of-function screens across human cancer cell lines, led by the Broad Institute Cancer Dependency Map (DepMap) (Tsherniak <em>et al.</em>, 2017; Meyers <em>et al.</em>, 2017; Dempster <em>et al.</em>, 2021), serve as a standard reference for oncology target prioritization, pathway mapping, and experimental model selection. Using copy-number- and kinetics-corrected Chronos gene effect scores (Dempster <em>et al.</em>, 2021), investigators routinely check whether candidate genes from transcriptomic, proteomic, or genetic screens act as pan-cancer essentials or context-restricted vulnerabilities.
  </p>
  <p>
    Although the DepMap Portal (<a href="https://depmap.org/portal/">https://depmap.org/portal/</a>) and companion resources such as shinyDepMap (Shimada <em>et al.</em>, 2021) and co-essentiality atlases (Wainberg <em>et al.</em>, 2021) handle routine single-gene lookups well, four practical gaps frequently appear during day-to-day panel prioritization:
  </p>
  <ol>
    <li><strong>Multi-gene batch context:</strong> Positioning a candidate panel of 5–50 genes together along a cohort-wide essentiality curve, or overlaying a prioritized subset onto a pathway background, usually requires repetitive single-gene lookups or ad-hoc R/Python scripting.</li>
    <li><strong>Linking lineage summaries to wet-lab model selection:</strong> Pan-cancer averages often dilute lineage- or subtype-restricted dependencies, whereas selecting concrete sensitive and resistant cell lines inside a chosen lineage scope requires cross-referencing multiple tables.</li>
    <li><strong>Direct cohort-scoped genome-wide co-dependency scanning:</strong> Default portal co-dependency tables report pan-cancer correlations in which inter-lineage variance can confound pathway co-essentiality; running a live single-lineage or subtype-restricted genome-wide scan typically requires multi-step Custom Analysis setup.</li>
    <li><strong>Reproducible publication export:</strong> Generating vector-editable figures (PDF/SVG) alongside exact numeric tables, data checksums, input-matching diagnostics, and environment metadata remains fragmented across web portals.</li>
  </ol>
  <p>
    We developed <strong>CRISPR Score Analyzer</strong>, an open-source browser application built around explicit input-state isolation, transparent descriptive statistics, colorblind-safe vector figure export, and self-contained provenance archiving.
  </p>

  <!-- Section 2 -->
  <h2><span class="sec-num">2.</span>Implementation and Core Capabilities</h2>

  <h3>2.1 Bundled Reference Data, Construction, Licensing, and Symbol Resolution</h3>
  <p style="text-align:left; hyphens:none;">
    CRISPR Score Analyzer is written in Python (Streamlit, pandas, NumPy, Plotly, Kaleido, and Pillow) with English and Chinese interface support. For deterministic startup without remote downloads, the repository bundles a compact <code>float32</code> binary matrix (<code>data/depmap_chronos.npz</code>, 88 MB; maximum absolute <code>float32</code> rounding deviation <span class="nowrap">1.19&nbsp;&times;&nbsp;10<sup>&minus;7</sup></span> across the first 100 genes) built from the <strong>DepMap Public 26Q3</strong> release (released October 5, 2026; <a href="https://forum.depmap.org/t/announcing-the-26q3-release/4692">https://forum.depmap.org/t/announcing-the-26q3-release/4692</a>; Data Portal: <a href="https://depmap.org/portal/data_page/">https://depmap.org/portal/data_page/</a>). Using <code>scripts/build_builtin_matrix.py</code>, official <code>CRISPRGeneEffect.csv</code> (497,139,964 bytes; MD5 <code>a0a87b48137996f620437604b8c2eb79</code>) and <code>Model.csv</code> (789,805 bytes; MD5 <code>a7423dbf345fcf5a822ab5fc7e854e5c</code>), retrieved on 2026-11-21, were validated for unique <code>ModelID</code> rows and <code>SYMBOL (Entrez ID)</code> headers. All <strong>1,361 cancer cell lines</strong> and <strong>18,434 genes</strong> in <code>CRISPRGeneEffect.csv</code> were retained verbatim without row filtering and left-joined on <code>ModelID</code> against deduplicated <code>Model.csv</code> to attach <code>OncotreeLineage</code> (1,361/1,361 matched across <strong>30 Oncotree lineages</strong>) and <code>CellLineName</code> (falling back to <code>ModelID</code> when blank), with SHA-256 and MD5 checksums recorded in <code>data/depmap_chronos.provenance.json</code>.
  </p>
  <p>
    Software and data licensing are kept strictly separate. The MIT license covers only application source code, whereas rehosting <code>data/depmap_chronos.npz</code> in the GitHub repository and hosted Streamlit instance is governed by the Broad Institute DepMap Terms and Conditions (non-commercial research use only; commercial use requires a separate Broad license; full terms are reposted in <code>data/DEPMAP_TERMS.md</code> and inside the hosted web app per DepMap’s <em>Continuity of Terms of Use</em> clause). MIT GitHub and Zenodo source release archives exclude the matrix via <code>.gitattributes</code>. Optional subtype filtering (<code>OncotreePrimaryDisease</code> / <code>primary_disease</code>) activates automatically when present in custom uploaded CSVs.
  </p>
  <p>
    Researchers can also upload custom UTF-8 CRISPR score CSV matrices. To prevent numeric metadata columns from being parsed as genes, automatic gene detection strictly requires DepMap-style <code>SYMBOL (Entrez ID)</code> headers; plain-symbol tables require explicit manual score-column selection in the sidebar. Ambiguous duplicate gene symbols and repeated metadata columns are excluded automatically and logged in a downloadable <code>diagnostics.csv</code> report.
  </p>
  <p>
    Gene-list inputs accept arbitrary ASCII or full-width delimiters (newlines, spaces, tabs, commas, and semicolons), <code>.txt/.csv/.tsv</code> file uploads, and pasted DepMap headers (e.g., <code>MYC (4609)</code>). An offline HGNC symbol map (<code>data/gene_aliases.json</code>, <strong>40,402</strong> unambiguous previous-symbol and alias mappings; Seal <em>et al.</em>, 2023) resolves legacy names (such as <span class="nowrap"><span class="gene">HER2</span>&nbsp;&rarr;&nbsp;<span class="gene">ERBB2</span></span>) automatically. Input matching reports resolved aliases, deduplicated entries, and unmatched queries explicitly before and after execution while leaving exact matches unflagged.
  </p>

  <h3>2.2 Interface Design and Modular Analytical Workflows</h3>
  <p>
    The browser interface (<strong>Figure 1</strong>) separates translucent navigation chrome from solid, zero-color-cast analytical surfaces across a neutral light theme (<code>#f5f5f7</code>) and a deep-charcoal dark theme (<code>#0f0f11</code>, avoiding pure-black halation on desktop monitors). A sticky four-module segmented dock remains accessible while scrolling through figures and tables, and all four modules share a single sidebar cohort scope (<strong>Cancer types</strong> and optional <strong>Disease / subtype</strong> for custom datasets; <strong>Table 1</strong>):
  </p>
  <ol>
    <li><strong>Gene ranking:</strong> Ranks all valid genes within the active cohort by ascending mean finite Chronos gene effect (<code>rank</code>, percentile&nbsp;=&nbsp;100&nbsp;&times;&nbsp;rank&nbsp;/&nbsp;<em>N</em>, <code>n_valid</code>, and missingness counts). Interactive scatter curves render background genes via WebGL (<code>Scattergl</code>) for responsiveness and convert all traces to vector paths during PDF/SVG export. Two configurable reference anchors (default <span class="gene">MYC</span> and <span class="gene">PTEN</span>, with alias resolution) and a dashed &minus;0.5 screening reference line orient the curve alongside collision-avoided transparent labels with leader lines.</li>
    <li><strong>Cancer-type comparison:</strong> Generates horizontal per-lineage boxplots showing valid sample counts (<em>n</em>) per group, sorted by median score or alphabetically, with optional 2-column grids, single-axis multi-gene comparison, and cell-line point overlays annotated with display name and <code>ModelID</code>. Submissions exceeding eight matched genes prompt explicit subset selection while preserving the full submitted list in metadata.</li>
    <li><strong>Gene-set annotation:</strong> Plots a background gene set beneath a prioritized highlight layer on the cohort ranking curve with configurable colors, enabling pathway-level positioning without conflating visual inspection with formal enrichment testing.</li>
    <li><strong>Co-dependency:</strong> Provides two complementary modes with instant switching between <strong>Pearson <em>r</em></strong> and <strong>Spearman <em>&rho;</em></strong>:
      <ul class="sublist">
        <li><em>Genome-wide scan (1 gene):</em> Correlates a single query gene (via a searchable dropdown supporting approved symbols and HGNC aliases) against every other ranked gene in the active sidebar scope using exact pairwise-complete finite observations in <code>float64</code> precision across 2,048-gene column chunks. Pairs must satisfy <code>min_shared</code>&nbsp;=&nbsp;min(10, max(3, <code>n_valid</code>)) (<span class="nowrap"><code>n_shared</code>&nbsp;&ge;&nbsp;10</span> for cohorts of at least 10 lines) and non-zero variance; Spearman <em>&rho;</em> resolves ties via average ranks using NumPy and pandas without SciPy. Outputs include a diverging top-partner lollipop chart, top-50 positive and negative tables showing <code>n_shared</code>, one-click transfer of top partners to a correlation matrix, and full-genome CSV export.</li>
        <li><em>Multi-gene matrix (2–25 genes):</em> Computes an interactive pairwise correlation heatmap paired with a <strong>lineage-colored cell-line scatter plot</strong>, ordinary least-squares trendline, and both Pearson <em>r</em> and Spearman <em>&rho;</em> statistics for any selected gene pair.</li>
      </ul>
    </li>
  </ol>

  <!-- Figure 1: UI Overview -->
  <div class="figure-block">
    <img src="{fig1_uri}" alt="Figure 1" />
    <p class="figure-caption">
      <strong>Figure 1. Interactive browser interface of CRISPR Score Analyzer across neutral light and deep-charcoal dark themes.</strong>
      <strong>(A)</strong> Light mode (<code>#f5f5f7</code> canvas) showing the sidebar cohort controls (<code>DepMap Public 26Q3 · Chronos Gene Effect</code>), sticky four-module navigation dock (<em>Gene ranking</em>, <em>Cancer-type comparison</em>, <em>Gene-set annotation</em>, <em>Co-dependency</em>), gene-list input with <em>Examples</em>/<em>Search</em> popovers, status strip, and solid white-paper NPG ranking curve (labeled <em>Mean CRISPR score</em> in the generic UI to accommodate custom uploads using alternative scoring metrics).
      <strong>(B)</strong> Scrolled deep-charcoal dark mode (<code>#0f0f11</code> canvas) showing the floating frosted navigation dock over the darkroom plot, colorblind-safe tabular summary (<span class="gene">KRAS</span>, <span class="gene">SOX10</span>, <span class="gene">BRAF</span>, <span class="gene">ESR1</span>, <span class="gene">TP53</span>), export controls, and contextual deep-dive tabs (<em>Lineage preference</em> and <em>Top cell lines</em>).
    </p>
  </div>

  <!-- Table 1 -->
  <div class="table-block">
    <p class="table-caption"><strong>Table 1.</strong> Workflow comparison across CRISPR Score Analyzer and existing DepMap web interfaces, verified November 2026 on equivalent gene-ranking, lineage-comparison, and co-dependency tasks.</p>
    <table class="cmp-table">
      <thead>
        <tr>
          <th>Workflow dimension</th>
          <th>DepMap Portal (26Q3; depmap.org/portal)</th>
          <th>shinyDepMap (eLife 2021; shinyapps.io/depmap)</th>
          <th class="col-csa">CRISPR Score Analyzer (v1.2.0 / <code>preprint-v2</code>)</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td class="cap-col">Multi-gene placement on cohort ranking curve</td>
          <td>Single-gene default (<code>Gene Page</code>)</td>
          <td>Efficacy vs selectivity scatter (<code>Efficacy &amp; Selectivity</code>)</td>
          <td class="csa-col">Batch ranking + two-layer gene-set annotation</td>
        </tr>
        <tr>
          <td class="cap-col">Cohort-scoped genome-wide co-dependency scan</td>
          <td>Pan-cancer default; subsets via multi-step <code>Custom Analysis</code></td>
          <td>Pan-cancer clusters (<code>Connectivity</code>)</td>
          <td class="csa-col">Single-step live scoped scan (Pearson <em>r</em> &amp; Spearman <em>&rho;</em>)</td>
        </tr>
        <tr>
          <td class="cap-col">Pairwise multi-gene matrix &amp; lineage scatter</td>
          <td><code>Data Explorer 2</code> / <code>Custom Analysis</code> small multiples</td>
          <td>Cluster network view</td>
          <td class="csa-col">2–25 gene heatmap + integrated lineage scatter</td>
        </tr>
        <tr>
          <td class="cap-col">Lineage preference vs all-lineage mean</td>
          <td>Per-gene lineage &amp; subtype distributions</td>
          <td>Not provided</td>
          <td class="csa-col">Diverging Top/Bottom-6 lollipop &amp; multi-gene heatmap</td>
        </tr>
        <tr>
          <td class="cap-col">Paired sensitive / resistant model lookup</td>
          <td>Table filtering on <code>Gene Page</code></td>
          <td>Top sensitive cell lines</td>
          <td class="csa-col">Paired most/least dependent tables + CSV export</td>
        </tr>
        <tr>
          <td class="cap-col">Custom CRISPR score matrix upload &amp; diagnostics</td>
          <td>Single-feature custom vector (not full matrix)</td>
          <td>Not supported</td>
          <td class="csa-col">Full UTF-8 CSV matrix + strict schema + <code>diagnostics.csv</code></td>
        </tr>
        <tr>
          <td class="cap-col">Offline batch HGNC alias resolution</td>
          <td>Online single-gene search</td>
          <td>Not provided</td>
          <td class="csa-col">40,402 unambiguous offline mappings (batch &amp; single)</td>
        </tr>
        <tr>
          <td class="cap-col">Colorblind-safe figures, themes &amp; provenance export</td>
          <td>Separate chart and CSV downloads</td>
          <td>Separate PNG and CSV</td>
          <td class="csa-col">NPG/<code>RdBu</code> palettes, light/dark UI, white-paper PDF/SVG/PNG + SHA-256 ZIP</td>
        </tr>
      </tbody>
    </table>
  </div>

  <h3>2.3 Deep-Dive Views, Publication Palettes, State Integrity, and Performance</h3>
  <p>
    Beneath ranking, comparison, and annotation results, two contextual deep-dive tabs support immediate follow-up without leaving the active module: (i) <strong>Lineage preference</strong>, comparing every cancer lineage (<span class="nowrap"><code>n_valid</code>&nbsp;&ge;&nbsp;5</span>) against each gene’s full-dataset mean (independent of the sidebar filter, outlining active lineages) in Top/Bottom-6 lollipop or diverging heatmap views alongside median difference and the fraction of lines scoring <span class="nowrap">&lt;&nbsp;&minus;0.5</span>; and (ii) <strong>Top cell lines</strong>, listing the most dependent and least dependent models in the active scope for experimental selection.
  </p>
  <p>
    All interactive plots and exports use colorblind-safe publication palettes: Nature Publishing Group (<code>NPG</code>) categorical hues (<code>#3C5488</code> deep indigo for query genes, <code>#E64B35</code> vermilion for <span class="gene">MYC</span>, <code>#00A087</code> teal for <span class="gene">PTEN</span>, and a 12-lineage NPG scatter palette, eliminating red–green conflicts) paired with a 7-stop ColorBrewer <code>RdBu</code> diverging scale (<code>#B2182B</code> deep crimson for positive correlation or less-dependent shifts; <code>#2166AC</code> deep steel blue for negative correlation or more-dependent shifts). All figure exports enforce white paper backgrounds, marker halos, transparent label boxes, and automatic dark-to-light <code>RdBu</code> conversion when exporting from dark mode, caching rendered bytes by exact figure fingerprint and dimensions (PNG output capped at 24 million pixels).
  </p>
  <p>
    To prevent mismatch between displayed outputs and edited inputs, gene-list edits remain uncommitted drafts until explicitly run, and changes to the dataset, column mapping, or sidebar scope invalidate prior results and exports. An on-demand <strong>Complete analysis ZIP</strong> bundles <code>selected_results.csv</code>, <code>full_rankings.csv</code>, <code>plot_data.csv</code>, <code>lineage_preference.csv</code>, <code>diagnostics.csv</code>, and <code>metadata.json</code> (recording data SHA-256, release version, Zenodo DOI, cohort parameters, input matching, and Python package versions).
  </p>
  <p>
    On a 2-vCPU Linux container (Python 3.12.16, NumPy 2.5.5, pandas 2.3.3, Streamlit 1.67.0, Plotly 6.9.0, Kaleido 0.2.1), median wall-clock times across three runs on the bundled 88 MB DepMap Public 26Q3 matrix are <strong>1.01 s</strong> for matrix loading (0.97–1.33 s), <strong>0.41 s</strong> for ranking all 18,434 genes (0.38–0.73 s), and <strong>0.82 s</strong> (Pearson <em>r</em>, 0.66–0.85 s) or <strong>6.21 s</strong> (Spearman <em>&rho;</em>, 6.12–6.77 s) for a pan-cancer genome-wide co-dependency scan across 18,433 partner genes. Analytical and export behavior is verified by <strong>90 automated tests</strong> covering CSV parsing, missing/infinite handling, exact pairwise correlation parity, UI state isolation, citation integrity, and vector/DPI figure exports via GitHub Actions continuous integration.
  </p>

  <!-- Section 3 -->
  <h2><span class="sec-num">3.</span>Benchmark Use Cases</h2>

  <h3>3.1 Separating Pan-Cancer Essentiality from Lineage-Selective Vulnerability</h3>
  <p>
    On initial launch with DepMap Public 26Q3 (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;1,361</span> cell lines; 1,712 genes with pan-cancer mean <span class="nowrap">&lt;&nbsp;&minus;0.5</span>, range &minus;5.250 to +0.398), CRISPR Score Analyzer displays a default five-gene panel (<span class="gene">KRAS</span>, <span class="gene">SOX10</span>, <span class="gene">BRAF</span>, <span class="gene">ESR1</span>, <span class="gene">TP53</span>) anchored by <span class="gene">MYC</span> (mean &minus;1.971, rank 187, 1.01%) and <span class="gene">PTEN</span> (mean +0.377, rank 18,428, 99.97%) (<strong>Figure 2A</strong>). Across all lineages, <span class="gene">KRAS</span> falls below the &minus;0.5 screening reference (mean &minus;0.729, rank 1,263, 6.85%), whereas tumor suppressor <span class="gene">TP53</span> sits at the positive extreme (mean +0.381, rank 18,430, 99.98%), consistent with prior DepMap reports of knockout-induced proliferation advantage in <span class="gene">TP53</span>-wild-type models (Tsherniak <em>et al.</em>, 2017; Dempster <em>et al.</em>, 2021), although somatic mutation calls are not bundled in the default matrix. By contrast, two established lineage-restricted dependencies appear modest or near-neutral in pan-cancer means: <span class="gene">SOX10</span> (mean &minus;0.347, rank 2,225, 12.07%) and <span class="gene">ESR1</span> (mean &minus;0.065, rank 7,347, 39.86%).
  </p>

  <!-- Figure 2: Ranking & Lineage Preference -->
  <div class="figure-block fig2">
    <img src="{fig2_uri}" alt="Figure 2" />
    <p class="figure-caption">
      <strong>Figure 2. Multi-gene Chronos gene effect positioning and lineage preference across DepMap Public 26Q3 (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;1,361</span> cell lines, 18,434 genes; Nature NPG and ColorBrewer <code>RdBu</code> palettes).</strong>
      <strong>(A)</strong> Pan-cancer mean Chronos gene effect ranking of five benchmark cancer genes (<span class="gene">KRAS</span>, <span class="gene">SOX10</span>, <span class="gene">BRAF</span>, <span class="gene">ESR1</span>, <span class="gene">TP53</span>) relative to <span class="gene">MYC</span> and <span class="gene">PTEN</span> reference anchors and the &minus;0.5 screening reference line.
      <strong>(B)</strong> Two-layer gene-set annotation highlighting <span class="gene">PLK1</span> and <span class="gene">AURKA</span> over seven cell-cycle CDK/cyclin background genes (<span class="gene">CDK1</span>, <span class="gene">CDK2</span>, <span class="gene">CDK4</span>, <span class="gene">CDK6</span>, <span class="gene">CCNB1</span>, <span class="gene">CCND1</span>, <span class="gene">CCNE1</span>; listed in panel legend).
      <strong>(C, D)</strong> Diverging lineage preference lollipops (top and bottom 4 lineages with <span class="nowrap"><em>n</em>&nbsp;&ge;&nbsp;5</span>, sample sizes shown on y-axis; note gene-specific x-axis scales) showing mean Chronos difference relative to each gene’s all-lineage mean for <span class="gene">SOX10</span> (<strong>C</strong>, Skin and Eye shifts) and <span class="gene">ESR1</span> (<strong>D</strong>, Breast alongside small-sample Ampulla of Vater <span class="nowrap"><em>n</em>&nbsp;=&nbsp;9</span>).
    </p>
  </div>

  <p>
    The <strong>Lineage preference</strong> tab separates mean shifts from sample-size effects and distributional skew (<strong>Figure 2C, D</strong>, noting gene-specific x-axis scales). Relative to its all-lineage mean (&minus;0.347), <span class="gene">SOX10</span> shows a consistent skin-lineage shift (<span class="nowrap">mean&nbsp;=&nbsp;&minus;1.487</span>, <span class="nowrap">&Delta;&nbsp;=&nbsp;&minus;1.140</span>, <span class="nowrap">median&nbsp;&Delta;&nbsp;=&nbsp;&minus;1.368</span>, <span class="nowrap">78.7%&nbsp;&lt;&nbsp;&minus;0.5</span>, <span class="nowrap"><em>n</em>&nbsp;=&nbsp;75</span>). In ocular models (Eye, <span class="nowrap"><em>n</em>&nbsp;=&nbsp;15</span>), the mean shift (<span class="nowrap">&Delta;&nbsp;=&nbsp;&minus;0.535</span>, <span class="nowrap">mean&nbsp;=&nbsp;&minus;0.882</span>) exceeds the median shift (<span class="nowrap">median&nbsp;&Delta;&nbsp;=&nbsp;&minus;0.174</span>), reflecting strong dependency restricted to a uveal melanoma subset (e.g., <code>Omm2.5</code> &minus;2.752, <code>UPMM3</code> &minus;2.113) alongside near-neutral retinoblastoma and retinal pigment epithelial lines. For <span class="gene">ESR1</span>, Breast (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;67</span>) yields a modest lineage-mean shift (<span class="nowrap">mean&nbsp;=&nbsp;&minus;0.285</span>, <span class="nowrap">&Delta;&nbsp;=&nbsp;&minus;0.220</span>, <span class="nowrap">median&nbsp;&Delta;&nbsp;=&nbsp;&minus;0.061</span>, <span class="nowrap">17.9%&nbsp;&lt;&nbsp;&minus;0.5</span>), followed closely in Figure 2D by a small Ampulla of Vater group (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;9</span>, <span class="nowrap">&Delta;&nbsp;=&nbsp;&minus;0.174</span>, <span class="nowrap">median&nbsp;&Delta;&nbsp;=&nbsp;&minus;0.059</span>). Because the bundled matrix stores only <code>OncotreeLineage</code> without breast receptor or PAM50 subtype annotations, lineage-level means compress ER-positive luminal dependency into a modest average. Checking the <strong>Top cell lines</strong> tab resolves this masking by surfacing established ER-positive models (<code>HCC1428</code> [<code>ACH-000352</code>, &minus;1.912], <code>T-47D</code> [<code>ACH-000147</code>, &minus;1.845], <code>MDA-MB-415</code> [<code>ACH-000876</code>, &minus;1.636]) alongside resistant negative controls, while custom CSV uploads containing <code>OncotreePrimaryDisease</code> enable subtype-level filtering. In the <strong>Gene-set annotation</strong> module, mitotic kinases <span class="gene">PLK1</span> (mean &minus;2.691, rank 38, 0.21%) and <span class="gene">AURKA</span> (mean &minus;1.236, rank 626, 3.40%) can be overlaid onto a seven-gene cell-cycle CDK/cyclin background (<span class="gene">CDK1</span> &minus;2.588 through <span class="gene">CCNE1</span> &minus;0.194; <strong>Figure 2B</strong>).
  </p>

  <h3>3.2 Benchmarking Pan-Cancer and Single-Lineage Co-Dependency Exploration</h3>
  <p>
    To evaluate co-dependency workflows in the <strong>Co-dependency</strong> module, we first ran a pan-cancer genome-wide Pearson scan of <span class="gene">MDM2</span> across all 18,433 other ranked genes (0.82 s; <strong>Figure 3A</strong>). The diverging lollipop recovers the <span class="gene">TP53</span>–<span class="gene">MDM2</span> regulatory axis: top positive co-dependencies include <span class="gene">PPM1D</span> (<span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.664</span>), <span class="gene">MDM4</span> (<span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.614</span>), <span class="gene">USP7</span> (<span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.522</span>), <span class="gene">DDX31</span> (<span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.466</span>), and <span class="gene">UBE2D3</span> (<span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.398</span>, all <span class="nowrap"><code>n_shared</code>&nbsp;=&nbsp;1,361</span>), whereas top anti-correlated partners include <span class="gene">TP53</span> (<span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.702</span>), <span class="gene">TP53BP1</span> (<span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.589</span>), <span class="gene">USP28</span> (<span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.531</span>), <span class="gene">CDKN1A</span> (<span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.530</span>), <span class="gene">CHEK2</span> (<span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.491</span>), and <span class="gene">ATM</span> (<span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.461</span>, all <span class="nowrap"><code>n_shared</code>&nbsp;=&nbsp;1,361</span>). Displaying <code>n_shared</code> directly in scan outputs flags <span class="gene">TERF1</span> (<span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.395</span>, rank 6 positive) as partially covered (<span class="nowrap"><code>n_shared</code>&nbsp;=&nbsp;228</span> vs 1,361), preventing sparse-coverage correlations from being weighted equally with full-panel hits. Transferring the six-gene regulatory panel (<span class="gene">TP53</span>, <span class="gene">CHEK2</span>, <span class="gene">MDM2</span>, <span class="gene">MDM4</span>, <span class="gene">PPM1D</span>, <span class="gene">USP7</span>) to multi-gene matrix mode resolves the two-block inverse correlation structure (<strong>Figure 3B</strong>), and the integrated cell-line scatter plot of <span class="gene">MDM2</span> vs <span class="gene">TP53</span> (<span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.702</span>, <span class="nowrap"><em>&rho;</em>&nbsp;=&nbsp;&minus;0.704</span>, <span class="nowrap"><em>n</em>&nbsp;=&nbsp;1,361</span>; <strong>Figure 3C</strong>) shows cell lines across multiple lineages following a similar negative trend. Offline verification (<code>preprint/reproduce_stats.py</code>) confirms that demeaning both genes by <code>OncotreeLineage</code> retains a strong within-lineage pooled correlation (<span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.648</span>, <span class="nowrap"><em>n</em>&nbsp;=&nbsp;1,361</span> across 30 lineages), and all 23 lineages with <span class="nowrap"><em>n</em>&nbsp;&ge;&nbsp;15</span> show negative correlations (median <span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.680</span>; e.g., Lung <span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.802</span> [<span class="nowrap"><em>n</em>&nbsp;=&nbsp;126</span>], Skin <span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.793</span> [<span class="nowrap"><em>n</em>&nbsp;=&nbsp;75</span>], Breast <span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.719</span> [<span class="nowrap"><em>n</em>&nbsp;=&nbsp;67</span>], Bowel <span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.674</span> [<span class="nowrap"><em>n</em>&nbsp;=&nbsp;85</span>], CNS/Brain <span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.663</span> [<span class="nowrap"><em>n</em>&nbsp;=&nbsp;133</span>], Lymphoid <span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.609</span> [<span class="nowrap"><em>n</em>&nbsp;=&nbsp;97</span>]).
  </p>
  <p>
    Pan-cancer correlations of tissue-enriched oncogenes, by contrast, mix pathway co-essentiality with inter-lineage variance. In a pan-cancer scan of <span class="gene">BRAF</span>, core MAPK components show moderate-to-strong descriptive correlations (<span class="gene">MAPK1</span> <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.644</span>, <span class="gene">MAP2K1</span> <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.391</span>, <span class="gene">MAP2K2</span> <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.271</span>) alongside melanoma lineage factors (<span class="gene">SOX10</span> <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.544</span>, <span class="gene">DUSP4</span> <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.474</span>; <span class="gene">SOX5</span> <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.208</span>). Filtering the sidebar scope to <strong>Skin</strong> (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;75</span> cell lines; 17,785 partner genes meeting <span class="nowrap"><code>n_shared</code>&nbsp;&ge;&nbsp;10</span>, with 648 sparser genes excluded) and re-running the scan produces top seven positive correlations for <span class="gene">MAPK1</span> (<span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.798</span>), <span class="gene">SOX10</span> (+0.633), <span class="gene">SOX5</span> (+0.562), <span class="gene">MAP2K1</span> (+0.552), <span class="gene">MAP2K2</span> (+0.550), <span class="gene">UBQLN1</span> (+0.526), and <span class="gene">AKIRIN1</span> (+0.496), with a negative tail led by <span class="gene">STYXL1</span> (&minus;0.534; all <span class="nowrap"><code>n_shared</code>&nbsp;=&nbsp;75</span>; <strong>Figure 3D</strong>). Restricting to Skin therefore increases descriptive RAF–MEK–ERK correlations relative to the pan-cancer scan. An offline OLS partial-correlation check regressing both <span class="gene">BRAF</span> and each partner’s Chronos gene effect on <span class="gene">SOX10</span> across all <span class="nowrap"><em>n</em>&nbsp;=&nbsp;75</span> complete Skin models (zero missing values; <code>preprint/reproduce_stats.py</code>) retains <span class="gene">MAPK1</span> (<span class="nowrap">partial <em>r</em>&nbsp;=&nbsp;+0.665</span>), <span class="gene">MAP2K1</span> (+0.467), and <span class="gene">MAP2K2</span> (+0.398), yet the top list is not a pure MAPK module. The persistence of <span class="gene">SOX10</span> and <span class="gene">SOX5</span> reflects Oncotree Skin heterogeneity spanning <span class="gene">SOX10</span>-high cutaneous melanoma and <span class="gene">SOX10</span>-low non-melanoma or dedifferentiated skin models (e.g., <code>HSC-5</code>, <code>MCC26</code>, <code>RVH421_DAB_R</code>). Post-hoc two-sided Pearson <em>t</em>-tests across all 17,785 Skin partner genes (<span class="nowrap"><code>n_shared</code>&nbsp;=&nbsp;75</span> each; <code>preprint/reproduce_stats.py</code>, computed outside the descriptive browser UI) further show how candidates near <span class="nowrap">|<em>r</em>|&nbsp;&approx;&nbsp;0.50</span> differ sharply from <span class="gene">MAPK1</span> (<span class="nowrap"><em>p</em>&nbsp;=&nbsp;9.72&nbsp;&times;&nbsp;10<sup>&minus;18</sup></span>, Bonferroni <span class="nowrap">1.73&nbsp;&times;&nbsp;10<sup>&minus;13</sup></span>): <span class="gene">AKIRIN1</span> (<span class="nowrap"><em>p</em>&nbsp;=&nbsp;5.95&nbsp;&times;&nbsp;10<sup>&minus;6</sup></span>, Bonferroni 0.106), <span class="gene">UBQLN1</span> (<span class="nowrap"><em>p</em>&nbsp;=&nbsp;1.29&nbsp;&times;&nbsp;10<sup>&minus;6</sup></span>, Bonferroni 0.023), and <span class="gene">STYXL1</span> (<span class="nowrap"><em>p</em>&nbsp;=&nbsp;8.28&nbsp;&times;&nbsp;10<sup>&minus;7</sup></span>, Bonferroni 0.015) remain exploratory leads requiring independent validation.
  </p>

  <!-- Figure 3: Co-dependency -->
  <div class="figure-block fig3">
    <img src="{fig3_uri}" alt="Figure 3" />
    <p class="figure-caption">
      <strong>Figure 3. Pan-cancer and single-lineage genome-wide co-dependency benchmarks (ColorBrewer <code>RdBu</code> correlation convention: positive crimson <code>#B2182B</code>, negative steel blue <code>#2166AC</code>; 12-lineage NPG scatter palette).</strong>
      <strong>(A)</strong> Pan-cancer genome-wide Pearson <em>r</em> scan of <span class="gene">MDM2</span> across 18,433 genes (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;1,361</span> cell lines), displaying the top 7 positive (<span class="gene">PPM1D</span>, <span class="gene">MDM4</span>, <span class="gene">USP7</span>, <span class="gene">DDX31</span>, <span class="gene">UBE2D3</span>, <span class="gene">TERF1</span> [<span class="nowrap"><em>n</em>&nbsp;=&nbsp;228</span>], <span class="gene">PPM1G</span>) and top 7 negative (<span class="gene">TP53</span>, <span class="gene">TP53BP1</span>, <span class="gene">USP28</span>, <span class="gene">CDKN1A</span>, <span class="gene">CHEK2</span>, <span class="gene">ATM</span>, <span class="gene">UBE2K</span>) partners.
      <strong>(B)</strong> Pairwise Pearson correlation matrix for the six-gene <span class="gene">TP53</span>–<span class="gene">MDM2</span> regulatory axis.
      <strong>(C)</strong> Lineage-colored cell-line scatter plot of <span class="gene">MDM2</span> vs <span class="gene">TP53</span> Chronos gene effect (<span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.702</span>, Spearman <span class="nowrap"><em>&rho;</em>&nbsp;=&nbsp;&minus;0.704</span>, <span class="nowrap"><em>n</em>&nbsp;=&nbsp;1,361</span>) with linear trendline.
      <strong>(D)</strong> Cohort-scoped genome-wide <span class="gene">BRAF</span> scan restricted to Skin cell lines (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;75</span>, 17,785 genes with <span class="nowrap"><code>n_shared</code>&nbsp;&ge;&nbsp;10</span>), showing higher descriptive MAPK cascade correlations (<span class="gene">MAPK1</span> +0.80, <span class="gene">MAP2K1</span> +0.55, <span class="gene">MAP2K2</span> +0.55) co-ranking with lineage factors (<span class="gene">SOX10</span> +0.63, <span class="gene">SOX5</span> +0.56) and exploratory hits near <span class="nowrap">|<em>r</em>|&nbsp;&approx;&nbsp;0.50</span> (all <span class="nowrap"><code>n_shared</code>&nbsp;=&nbsp;75</span>).
    </p>
  </div>

  <!-- Section 4 -->
  <div class="keep-together">
    <h2><span class="sec-num">4.</span>Scope, Limitations, and Scientific Interpretation</h2>
    <p>
      CRISPR Score Analyzer is built for exploratory visualization rather than formal statistical hypothesis testing. Four interpretive limits apply to all outputs:
    </p>
    <ul class="scope-list">
      <li><strong>Descriptive metrics only:</strong> Cohort gene ranks, lineage preference mean/median differences, and co-dependency correlation coefficients are descriptive summaries. They do not represent formal tests of selectivity, differential dependency, multiple-testing-corrected significance, or synthetic lethality, and should be followed by statistical modeling or experimental validation.</li>
      <li><strong>Reference threshold:</strong> The &minus;0.5 dashed line serves solely as a visual screening reference on Chronos gene effect scores, not as a formal classifier of common-essential genes.</li>
      <li><strong>Confounding and multiplicity in co-dependency scans:</strong> High correlation between two genes can reflect physical complex membership, shared signaling pathways, residual subtype or cell-state heterogeneity, common-essentiality variance, or small-sample multiplicity when scanning ~18,000 genes within one lineage. Checking shared sample counts (<code>n_shared</code>), Spearman <em>&rho;</em> alongside Pearson <em>r</em>, and lineage-colored cell-line scatter plots is necessary before inferring direct genetic interaction.</li>
      <li><strong>Bundled metadata and custom upload assumptions:</strong> The bundled DepMap Public 26Q3 matrix includes <code>OncotreeLineage</code> without mutation, expression, or PAM50/receptor subtype columns; users requiring subtype stratification should upload custom CSVs containing <code>OncotreePrimaryDisease</code> or equivalent metadata. Uploaded score matrices are assumed to contain one row per unique cell line; duplicate rows are not silently deduplicated.</li>
    </ul>
  </div>

  <div class="keep-together">
    <h2>Acknowledgements and Disclosures</h2>
    <div class="disclosure-block">
      <p><strong>Acknowledgements:</strong> The author thanks Dr. Lin Deng and Dr. Zixi Wang (Institute of Molecular Physiology, Shenzhen Bay Laboratory) for their constructive suggestions and feedback.</p>
      <p><strong>Data Source and Terms:</strong> Reference dependency data derive from DepMap Public 26Q3 Chronos Gene Effect scores (released October 5, 2026; <a href="https://forum.depmap.org/t/announcing-the-26q3-release/4692">https://forum.depmap.org/t/announcing-the-26q3-release/4692</a>; Data Portal: <a href="https://depmap.org/portal/data_page/">https://depmap.org/portal/data_page/</a>). The bundled matrix is governed by Broad Institute DepMap Terms and Conditions (non-commercial research use only; reposted in full in <code>data/DEPMAP_TERMS.md</code> and inside the hosted web app), separate from the MIT-licensed application code.</p>
      <p><strong>AI Assistance Disclosure:</strong> Large language model assistants (Claude, Codex, and related coding agents) were used to assist with software engineering, automated testing, and manuscript editing. All scientific specifications, biological case studies, verification, and final manuscript decisions were reviewed and approved by the author.</p>
      <p><strong>Funding:</strong> This work received no specific external funding.</p>
      <p><strong>Competing Interests:</strong> The author declares no competing financial or non-financial interests.</p>
      <p><strong>Data and Code Availability:</strong> The MIT-licensed software release (<code>v1.2.0</code>, commit <code>5c18f9b</code>) is hosted at <a href="https://github.com/ChanghaoKan/crispr-score-analyzer/releases/tag/v1.2.0">https://github.com/ChanghaoKan/crispr-score-analyzer/releases/tag/v1.2.0</a> and archived on Zenodo (version 1.2.0 DOI: <a href="https://doi.org/10.5281/zenodo.24458535"><code>10.5281/zenodo.24458535</code></a>; concept DOI: <code>10.5281/zenodo.19607602</code>); the hosted web tool is accessible at <a href="https://crispr-score-analyzer-szbl-denglab.streamlit.app/">https://crispr-score-analyzer-szbl-denglab.streamlit.app/</a>. Because MIT release source archives exclude <code>data/depmap_chronos.npz</code> under DepMap terms, the updated UI/CNS plotting palettes together with the manuscript figure and statistical reproduction scripts (<code>preprint/generate_figures.py</code> and <code>preprint/reproduce_stats.py</code>) are pinned at repository tag <code>preprint-v2</code> (analytical computation kernels are unchanged from <code>v1.2.0</code>). From a clean Python 3.12 environment, all quantitative figures (Figure 2 and Figure 3, plus re-compositing Figure 1 from bundled UI captures <code>ui_light.png</code> and <code>ui_dark.png</code>) and statistics can be reproduced via <code class="wrap-code">git clone https://github.com/ChanghaoKan/crispr-score-analyzer.git &amp;&amp; cd crispr-score-analyzer &amp;&amp; git checkout preprint-v2 &amp;&amp; python -m venv .venv &amp;&amp; .venv/bin/pip install -r requirements.txt &amp;&amp; .venv/bin/python preprint/generate_figures.py &amp;&amp; .venv/bin/python preprint/reproduce_stats.py</code>, verifying <code>data/depmap_chronos.npz</code> against <code>data/depmap_chronos.provenance.json</code> (<code>matrix_sha256</code>: <code class="wrap-code">ddf274da2aa5f89cf630bc8150dcd3a2c722466bd18ff6fb6d18fe161f8a0c67</code>, or rebuilding from official 26Q3 files via <code>scripts/build_builtin_matrix.py</code>).</p>
    </div>
  </div>

  <h2>References</h2>
  <ul class="references">
    <li>Dempster, J. M., Boyle, I., Vazquez, F., Root, D. E., Boehm, J. S., Hahn, W. C., Tsherniak, A., &amp; McFarland, J. M. (2021). Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. <em>Genome Biology</em>, 22(1), 343. <a href="https://doi.org/10.1186/s13059-021-02540-7">https://doi.org/10.1186/s13059-021-02540-7</a></li>
    <li>DepMap, Broad. (2026). DepMap Public 26Q3 Chronos Gene Effect dataset (released October 5, 2026). Broad Institute. <a href="https://depmap.org/portal/data_page/">https://depmap.org/portal/data_page/</a> (Release notes: <a href="https://forum.depmap.org/t/announcing-the-26q3-release/4692">https://forum.depmap.org/t/announcing-the-26q3-release/4692</a>)</li>
    <li>Kan, C. (2026). CRISPR Score Analyzer (1.2.0) [Software]. Zenodo. <a href="https://doi.org/10.5281/zenodo.24458535">https://doi.org/10.5281/zenodo.24458535</a> (Concept DOI: <a href="https://doi.org/10.5281/zenodo.19607602">https://doi.org/10.5281/zenodo.19607602</a>)</li>
    <li>Meyers, R. M., Bryan, J. G., McFarland, J. M., Weir, B. A., Sizemore, A. E., Xu, H., <em>et al.</em> (2017). Computational correction of copy number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells. <em>Nature Genetics</em>, 49(12), 1779–1784. <a href="https://doi.org/10.1038/ng.3984">https://doi.org/10.1038/ng.3984</a></li>
    <li>Seal, R. L., Braschi, B., Gray, K., Jones, T. E. M., Tweedie, S., Haim-Vilmovsky, L., &amp; Bruford, E. A. (2023). Genenames.org: the HGNC resources in 2023. <em>Nucleic Acids Research</em>, 51(D1), D1003–D1009. <a href="https://doi.org/10.1093/nar/gkac888">https://doi.org/10.1093/nar/gkac888</a></li>
    <li>Shimada, K., Bachman, J. A., Muhlich, J. L., &amp; Mitchison, T. J. (2021). shinyDepMap, a tool to identify targetable cancer genes and their functional connections from Cancer Dependency Map data. <em>eLife</em>, 10, e57116. <a href="https://doi.org/10.7554/eLife.57116">https://doi.org/10.7554/eLife.57116</a></li>
    <li>Tsherniak, A., Vazquez, F., Montgomery, P. G., Weir, B. A., Kryukov, G., Cowley, G. S., <em>et al.</em> (2017). Defining a Cancer Dependency Map. <em>Cell</em>, 170(3), 564–576.e16. <a href="https://doi.org/10.1016/j.cell.2017.06.010">https://doi.org/10.1016/j.cell.2017.06.010</a></li>
    <li>Wainberg, M., Kamber, R. A., Balsubramani, A., Meyers, R. M., Sinnott-Armstrong, N., Hornburg, D., <em>et al.</em> (2021). A genome-wide atlas of co-essential modules assigns function to uncharacterized genes. <em>Nature Genetics</em>, 53(5), 638–649. <a href="https://doi.org/10.1038/s41588-021-00840-z">https://doi.org/10.1038/s41588-021-00840-z</a></li>
  </ul>

</body>
</html>
"""


def main() -> None:
    html = build_html()
    with tempfile.NamedTemporaryFile("w", suffix=".html", encoding="utf-8", delete=False) as tf:
        tf.write(html)
        html_path = Path(tf.name)

    header_tpl = """
    <div style="width:100%; font-family:'Liberation Sans',Arial,sans-serif; font-size:7pt; color:#86868b; padding:0 17.5mm; display:flex; justify-content:space-between; align-items:center;">
      <span>Kan (2026) &middot; <em>CRISPR Score Analyzer</em> (v1.2.0)</span>
      <span>OSF Preprints</span>
    </div>
    """
    footer_tpl = """
    <div style="width:100%; font-family:'Liberation Sans',Arial,sans-serif; font-size:7pt; color:#86868b; padding:0 17.5mm; display:flex; justify-content:space-between; align-items:center;">
      <span>DepMap Public 26Q3 &middot; Open-Source Application Note</span>
      <span>Page <span class="pageNumber"></span> of <span class="totalPages"></span></span>
    </div>
    """

    try:
        with sync_playwright() as p:
            launch_kwargs = {"args": ["--no-sandbox"]}
            if Path("/usr/bin/google-chrome").exists():
                launch_kwargs["executable_path"] = "/usr/bin/google-chrome"
            browser = p.chromium.launch(**launch_kwargs)
            page = browser.new_page()
            page.goto(html_path.as_uri(), wait_until="networkidle")
            page.pdf(
                path=str(OUTPUT_PDF),
                format="A4",
                print_background=True,
                display_header_footer=True,
                header_template=header_tpl,
                footer_template=footer_tpl,
                margin={
                    "top": "15.5mm",
                    "bottom": "15.5mm",
                    "left": "17.5mm",
                    "right": "17.5mm",
                },
            )
            browser.close()
    finally:
        html_path.unlink(missing_ok=True)

    size_kb = OUTPUT_PDF.stat().st_size / 1024.0
    print(f"Wrote {OUTPUT_PDF} ({size_kb:.1f} KB)")


if __name__ == "__main__":
    main()
