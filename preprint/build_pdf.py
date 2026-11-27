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
    margin: 14mm 17mm 14mm 17mm;
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
    font-size: 8.35pt;
    line-height: 1.38;
    -webkit-font-smoothing: antialiased;
    text-rendering: optimizeLegibility;
  }}

  /* Header / Title Block (editorial NPG Deep Indigo #3C5488 accent) */
  .doc-kicker {{
    display: flex;
    justify-content: space-between;
    align-items: baseline;
    font-size: 7.1pt;
    font-weight: 700;
    letter-spacing: 0.085em;
    text-transform: uppercase;
    color: #3C5488;
    padding-bottom: 2.8pt;
    border-bottom: 1.1pt solid #1d1d1f;
    margin-bottom: 6pt;
  }}

  .doc-kicker .kicker-right {{
    color: #6e6e73;
    font-weight: 600;
    letter-spacing: 0.055em;
  }}

  h1.title {{
    font-size: 14.8pt;
    line-height: 1.19;
    font-weight: 700;
    letter-spacing: -0.018em;
    color: #1d1d1f;
    margin: 0 0 4.5pt 0;
  }}

  .author-block {{
    margin: 0 0 6.2pt 0;
  }}

  .author-name {{
    font-size: 9.2pt;
    font-weight: 700;
    color: #1d1d1f;
    margin: 0 0 1.4pt 0;
  }}

  .affiliations {{
    font-size: 7.75pt;
    line-height: 1.35;
    color: #6e6e73;
    margin: 0;
    text-align: left;
  }}

  /* Structured Abstract Card */
  .abstract-card {{
    background: #f5f5f7;
    border-left: 2.4pt solid #3C5488;
    border-radius: 0 4pt 4pt 0;
    padding: 6.2pt 9.5pt 5.8pt 9.5pt;
    margin: 0 0 7.5pt 0;
  }}

  .abstract-label {{
    font-size: 7.2pt;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 0.085em;
    color: #3C5488;
    margin: 0 0 2.8pt 0;
  }}

  .abstract-card p {{
    font-size: 7.8pt;
    line-height: 1.35;
    color: #1d1d1f;
    margin: 0 0 3.3pt 0;
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
    font-size: 7.5pt !important;
    color: #424245 !important;
    padding-top: 2.8pt;
    border-top: 0.5pt solid #d2d2d7;
    margin-top: 3.5pt !important;
    text-align: left !important;
  }}

  /* Sections */
  h2 {{
    font-size: 10.2pt;
    font-weight: 700;
    letter-spacing: -0.01em;
    color: #1d1d1f;
    margin: 8.5pt 0 3.8pt 0;
    padding-bottom: 2pt;
    border-bottom: 0.6pt solid #d2d2d7;
    break-after: avoid;
    page-break-after: avoid;
  }}

  h2 .sec-num {{
    color: #3C5488;
    margin-right: 4pt;
  }}

  h3 {{
    font-size: 8.95pt;
    font-weight: 700;
    color: #1d1d1f;
    margin: 7pt 0 2.8pt 0;
    break-after: avoid;
    page-break-after: avoid;
  }}

  p {{
    margin: 0 0 4.2pt 0;
    color: #242426;
    text-align: justify;
    hyphens: auto;
    orphans: 3;
    widows: 3;
  }}

  /* Inline & Display Mathematical Typography */
  .math {{
    font-family: "Liberation Serif", "Times New Roman", Georgia, serif;
    font-size: 8.95pt;
    color: #1d1d1f;
    white-space: nowrap;
  }}

  .math i {{
    font-style: italic;
  }}

  .math sub, .math sup {{
    font-size: 68%;
    line-height: 0;
  }}

  .eq-row {{
    display: flex;
    align-items: center;
    justify-content: space-between;
    margin: 3.8pt 0 4.6pt 0;
    padding: 4.5pt 8pt;
    background: #fafafa;
    border-top: 0.45pt solid #e5e5ea;
    border-bottom: 0.45pt solid #e5e5ea;
    break-inside: avoid;
    page-break-inside: avoid;
  }}

  .eq-body {{
    flex: 1;
    text-align: center;
    font-family: "Liberation Serif", "Times New Roman", Georgia, serif;
    font-size: 9.0pt;
    line-height: 1.3;
    color: #1d1d1f;
  }}

  .eq-num {{
    font-family: "Liberation Sans", Arial, sans-serif;
    font-size: 7.8pt;
    color: #6e6e73;
    padding-left: 8pt;
    white-space: nowrap;
  }}

  .eq-gap {{
    display: inline-block;
    width: 16pt;
  }}

  .frac {{
    display: inline-flex;
    flex-direction: column;
    align-items: center;
    vertical-align: middle;
    margin: 0 2pt;
    line-height: 1.15;
  }}

  .frac .num {{
    padding: 1.2pt 3pt 1.8pt 3pt;
    border-bottom: 0.6pt solid #1d1d1f;
    font-size: 8.4pt;
  }}

  .frac .den {{
    padding: 1.8pt 3pt 0.5pt 3pt;
    font-size: 8.4pt;
  }}

  .frac.frac-tall .num {{
    padding: 1.5pt 4pt 5.5pt 4pt;
  }}

  .frac.frac-tall .den {{
    padding: 3.5pt 4pt 4pt 4pt;
  }}

  .sum-op {{
    display: inline-flex;
    flex-direction: column;
    align-items: center;
    vertical-align: -0.32em;
    margin: 0 2pt;
    line-height: 0.95;
  }}

  .sum-op .sym {{
    font-size: 11.8pt;
    font-style: normal;
  }}

  .sum-op .lim {{
    font-size: 6.3pt;
    margin-top: 1.4pt;
  }}

  .overline {{
    display: inline-block;
    border-top: 0.55pt solid #1d1d1f;
    padding-top: 0.4pt;
    line-height: 0.92;
  }}

  /* Inline formatting */
  .gene {{
    font-style: italic;
  }}

  .nowrap {{
    white-space: nowrap;
  }}

  code {{
    font-family: "Liberation Mono", "Courier New", monospace;
    font-size: 7.35pt;
    background: #f5f5f7;
    color: #1d1d1f;
    padding: 0.2pt 1.8pt;
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

  sup, sub {{
    font-size: 68%;
    line-height: 0;
  }}

  /* Academic Tables */
  .table-block {{
    margin: 6pt 0 7pt 0;
    break-inside: avoid;
    page-break-inside: avoid;
  }}

  .table-caption {{
    font-size: 7.8pt;
    line-height: 1.34;
    color: #1d1d1f;
    margin: 0 0 3.5pt 0;
    text-align: left;
  }}

  table.cmp-table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 7.25pt;
    line-height: 1.28;
    border-top: 1.1pt solid #1d1d1f;
    border-bottom: 1.1pt solid #1d1d1f;
  }}

  table.cmp-table thead th {{
    background: #f5f5f7;
    color: #1d1d1f;
    font-weight: 700;
    text-align: left;
    padding: 3.4pt 4.5pt;
    border-bottom: 0.7pt solid #86868b;
    vertical-align: bottom;
  }}

  table.cmp-table thead th.num-col,
  table.cmp-table tbody td.num-col {{
    text-align: right;
    font-variant-numeric: tabular-nums;
    white-space: nowrap;
  }}

  table.cmp-table thead th.col-csa {{
    color: #3C5488;
  }}

  table.cmp-table tbody td {{
    padding: 2.9pt 4.5pt;
    border-bottom: 0.4pt solid #e5e5ea;
    vertical-align: top;
    color: #242426;
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
    width: 25%;
  }}

  table.cmp-table td.csa-col {{
    color: #1d1d1f;
    width: 31%;
  }}

  /* Figures */
  .figure-block {{
    margin: 5.5pt 0 7pt 0;
    padding: 5pt 6pt 5pt 6pt;
    border: 0.5pt solid #e5e5ea;
    border-radius: 4pt;
    background: #ffffff;
    break-inside: avoid;
    page-break-inside: avoid;
  }}

  .figure-block img {{
    display: block;
    width: 98.5%;
    height: auto;
    margin: 0 auto 4.5pt auto;
  }}

  .figure-block.fig2 img {{
    width: 90%;
  }}

  .figure-block.fig3 img {{
    width: 90%;
  }}

  .figure-caption {{
    font-size: 7.55pt;
    line-height: 1.35;
    color: #424245;
    margin: 0;
    text-align: justify;
    hyphens: auto;
  }}

  .figure-caption strong {{
    color: #1d1d1f;
  }}

  /* Disclosures */
  .disclosure-block p {{
    font-size: 7.55pt;
    line-height: 1.31;
    margin-bottom: 2.1pt;
    text-align: left;
    hyphens: none;
  }}

  /* Author-year hanging-indent References */
  ul.references {{
    list-style: none;
    padding-left: 0;
    margin: 2.5pt 0 0 0;
  }}

  ul.references li {{
    font-size: 7.25pt;
    line-height: 1.28;
    color: #424245;
    padding-left: 13pt;
    text-indent: -13pt;
    margin-bottom: 1.7pt;
    text-align: left;
    hyphens: none;
  }}
</style>
</head>
<body>

  <!-- Header & Title -->
  <div class="doc-kicker">
    <span>Bioinformatics &middot; Software &amp; Methods</span>
    <span class="kicker-right">Preprint</span>
  </div>

  <h1 class="title">CRISPR Score Analyzer: cohort-scoped essentiality and co-dependency analysis of DepMap CRISPR-Cas9 screens</h1>

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
      <strong>Background:</strong> Genome-scale CRISPR-Cas9 loss-of-function screens in the Cancer Dependency Map (DepMap) provide a standard reference for prioritizing candidate cancer vulnerabilities and investigating functional gene relationships. When exploratory analyses move beyond single-gene portal lookups, three recurrent methodological issues arise: pan-cancer averages dilute lineage- or subtype-restricted dependencies, multi-library integration introduces uneven gene coverage across cell lines, and global pan-cancer co-essentiality profiles conflate pathway co-functionality with baseline differences between tumor lineages.
    </p>
    <p>
      <strong>Results:</strong> We developed <strong>CRISPR Score Analyzer</strong>, an open-source Python package (<code>crispr_score_analyzer</code>), command-line interface (<code>crispr-score-analyzer</code>), and interactive browser platform built on a shared analytical core over the DepMap Public 26Q3 Chronos dataset (1,361 cancer cell lines across 30 Oncotree lineages and 18,434 genes) and validated user-supplied matrices. The analytical engine implements cohort-scoped essentiality ranking, joint mean and median lineage-selectivity decomposition, exact pairwise-complete <code>float64</code> Pearson and average-rank Spearman genome-wide scans in 2,048-column blocks without SciPy, offline HGNC alias resolution (40,402 unambiguous mappings), vector figure export, and SHA-256 provenance archiving. Across DepMap Public 26Q3, 1,352 genes (7.33%) have incomplete cell-line coverage (<code>n_valid</code> down to 218), motivating explicit pairwise overlap reporting (<code>n_shared</code>) alongside every correlation. Chunked correlation kernels agree with direct reference calculations to machine precision (maximum absolute difference <span class="nowrap">7.77&nbsp;&times;&nbsp;10<sup>&minus;16</sup></span> for Pearson <em>r</em> and <span class="nowrap">5.55&nbsp;&times;&nbsp;10<sup>&minus;17</sup></span> for Spearman <em>&rho;</em>) and execute a pan-cancer scan in 0.74&nbsp;s (Pearson) or 0.10&nbsp;s within a single lineage (Skin, <span class="nowrap"><em>n</em>&nbsp;=&nbsp;75</span>). Benchmarking ten canonical cancer gene pairs distinguishes lineage-invariant epistatic axes (<span class="gene">MDM2</span>–<span class="gene">TP53</span>, pan-cancer <span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.702</span>, lineage-demeaned <span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.648</span>) from context-restricted signaling and transcriptional modules whose pan-cancer rank correlations collapse toward zero (<span class="gene">BRAF</span>–<span class="gene">MAP2K1</span> pan-cancer <span class="nowrap"><em>&rho;</em>&nbsp;=&nbsp;+0.097</span> vs. Skin <span class="nowrap"><em>&rho;</em>&nbsp;=&nbsp;+0.605</span>; <span class="gene">ESR1</span>–<span class="gene">GATA3</span> pan-cancer <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.253</span> vs. Breast <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.596</span>).
    </p>
    <p class="avail">
      <strong>Availability:</strong> Application source code is released under the MIT license at <a href="https://github.com/ChanghaoKan/crispr-score-analyzer">https://github.com/ChanghaoKan/crispr-score-analyzer</a> (archived software release <code>v1.2.0</code>, commit <code>5c18f9b</code>, Zenodo DOI: <a href="https://doi.org/10.5281/zenodo.24458535"><code>10.5281/zenodo.24458535</code></a>, concept DOI: <code>10.5281/zenodo.19607602</code>), with a hosted web instance at <a href="https://crispr-score-analyzer-szbl-denglab.streamlit.app/">https://crispr-score-analyzer-szbl-denglab.streamlit.app/</a>. Python distribution packages and MIT source archives exclude the bundled DepMap Public 26Q3 matrix, which is governed separately by Broad Institute DepMap Terms and Conditions (non-commercial research use only; <code>data/DEPMAP_TERMS.md</code>). This manuscript, the modular Python library and CLI, and all figure and statistical reproduction scripts (<code>preprint/generate_figures.py</code> and <code>preprint/reproduce_stats.py</code>) are pinned at repository tag <code>preprint-v4</code>.
    </p>
    <p class="keywords">
      <strong>Keywords:</strong> DepMap; CRISPR-Cas9 screen; Chronos gene effect; gene essentiality; co-dependency; functional genomics
    </p>
  </div>

  <!-- Section 1 -->
  <h2><span class="sec-num">1.</span>Introduction</h2>
  <p>
    Pooled genome-scale CRISPR-Cas9 loss-of-function screens across human cancer cell lines have become a primary experimental reference for identifying lineage-specific vulnerabilities and mapping functional gene relationships (Tsherniak <em>et al.</em>, 2017; Meyers <em>et al.</em>, 2017; Pacini <em>et al.</em>, 2021). Within the Cancer Dependency Map (DepMap), copy-number- and growth-kinetics-corrected Chronos gene effect scores (Dempster <em>et al.</em>, 2021) allow candidate genes nominated by transcriptomics, proteomics, or focused genetic screens to be positioned rapidly against pan-cancer essentiality distributions.
  </p>
  <p>
    Routine exploration of these matrices, however, exposes three recurring analytical limitations once analysis moves beyond single-gene portal lookups. First, pan-cancer summary statistics frequently obscure lineage-restricted and subtype-restricted dependencies: an oncogene or lineage transcription factor strongly essential in one tumor type can appear near-neutral when averaged across more than a thousand models, while intra-lineage heterogeneity can induce substantial divergence between lineage means and medians. Second, integrating CRISPR screens conducted with different sgRNA libraries introduces uneven gene coverage across cell-line subsets (Pacini <em>et al.</em>, 2021), so genome-wide correlation rankings computed without tracking pairwise sample overlap can place sparsely profiled genes alongside full-cohort associations. Third, pan-cancer co-essentiality profiles conflate direct pathway co-functionality with baseline mean differences across tumor lineages. Global covariance-whitening and dimensionality-reduction methods address correlated sample structure when constructing whole-genome co-essentiality networks (Boyle <em>et al.</em>, 2018; Wainberg <em>et al.</em>, 2021; Gheorghe &amp; Hart, 2022), yet rely on offline matrix transformations rather than direct inspection within a user-defined lineage or disease subset. Existing interactive platforms, including the DepMap Portal (<a href="https://depmap.org/portal/">https://depmap.org/portal/</a>) and shinyDepMap (Shimada <em>et al.</em>, 2021), serve single-gene inspection and global overviews well, whereas multi-gene batch positioning, single-lineage co-dependency scanning, cell-line model selection, and reproducible scripted execution remain fragmented across separate portals and custom notebooks.
  </p>
  <p>
    Here we describe <strong>CRISPR Score Analyzer</strong>, an open-source Python library (<code>crispr_score_analyzer</code>), command-line interface (<code>crispr-score-analyzer</code>), and interactive browser platform sharing a single analytical core over the DepMap Public 26Q3 Chronos dataset (1,361 cell lines across 30 Oncotree lineages and 18,434 genes) as well as validated user-supplied CSV matrices. Rather than introducing a latent-factor normalization model, the software provides transparent, cohort-scoped descriptive estimators: exact pairwise-complete Pearson and Spearman scanning without imputation, joint mean and median lineage-selectivity summaries linked directly to individual cell-line identities, offline HGNC symbol harmonization, and self-contained provenance bundles. Using DepMap Public 26Q3, we audit gene-coverage missingness, verify numerical accuracy and runtime scaling of the chunked correlation kernels, and systematically evaluate pan-cancer versus cohort-scoped co-dependency across ten canonical cancer gene pairs spanning pan-lineage epistatic axes, lineage-restricted signaling cascades, and subtype-dependent transcriptional modules.
  </p>

  <!-- Section 2 -->
  <h2><span class="sec-num">2.</span>Materials and Methods</h2>

  <h3>2.1 Reference Matrix Curation, Schema Validation, and Nomenclature Harmonization</h3>
  <p>
    To enable reproducible local and hosted execution without external network fetches at runtime, the repository includes a compressed <code>float32</code> binary matrix (<code>data/depmap_chronos.npz</code>, 88&nbsp;MB) compiled from <strong>DepMap Public 26Q3</strong> (released October 5, 2026; <a href="https://forum.depmap.org/t/announcing-the-26q3-release/4692">https://forum.depmap.org/t/announcing-the-26q3-release/4692</a>; Data Portal: <a href="https://depmap.org/portal/data_page/">https://depmap.org/portal/data_page/</a>). Using <code>scripts/build_builtin_matrix.py</code>, the upstream <code>CRISPRGeneEffect.csv</code> (497,139,964&nbsp;bytes; MD5 <code>a0a87b48137996f620437604b8c2eb79</code>) and <code>Model.csv</code> (789,805&nbsp;bytes; MD5 <code>a7423dbf345fcf5a822ab5fc7e854e5c</code>), retrieved on 2026-11-21, were validated for unique <code>ModelID</code> rows and DepMap-formatted <code>SYMBOL (Entrez ID)</code> headers. All 1,361 cancer cell lines and 18,434 genes were retained verbatim without row filtering and left-joined on <code>ModelID</code> against deduplicated <code>Model.csv</code> to attach <code>OncotreeLineage</code> (1,361/1,361 matched across 30 Oncotree lineages) and <code>CellLineName</code> (defaulting to <code>ModelID</code> where blank). Storage quantization to <code>float32</code> introduces a maximum absolute rounding error of <span class="nowrap">1.19&nbsp;&times;&nbsp;10<sup>&minus;7</sup></span> across the first 100 genes while reducing file size more than fivefold; upstream and compiled SHA-256/MD5 digests are recorded in <code>data/depmap_chronos.provenance.json</code>.
  </p>
  <p>
    Software and dataset licensing are kept strictly separate. The MIT license governs application source code and the CC0/public-domain HGNC alias map bundled in Python wheels, whereas rehosting <code>data/depmap_chronos.npz</code> in the GitHub repository and hosted web instance falls under the Broad Institute DepMap Terms and Conditions (non-commercial research use only; full text reproduced in <code>data/DEPMAP_TERMS.md</code> and inside the web interface under DepMap’s <em>Continuity of Terms of Use</em> clause). Python wheel/sdist builds (<code>pyproject.toml</code>) and MIT GitHub/Zenodo source archives (<code>.gitattributes</code>) exclude <code>data/depmap_chronos.npz</code>.
  </p>
  <p>
    For custom UTF-8 CSV CRISPR score matrices, automatic gene-column recognition strictly requires <code>SYMBOL (Entrez ID)</code> headers so that numeric metadata fields are not silently parsed as genes; plain-symbol tables require explicit column selection through the sidebar or CLI (<code>--gene-columns</code>). Ambiguous duplicate gene symbols and repeated metadata columns are excluded automatically and itemized in a downloadable <code>diagnostics.csv</code> file. Optional disease-subtype filtering (<code>OncotreePrimaryDisease</code> or <code>primary_disease</code>) activates automatically whenever present in uploaded tables. Gene queries accept arbitrary whitespace or punctuation delimiters, <code>.txt</code>/<code>.csv</code>/<code>.tsv</code> lists, and pasted DepMap headers (e.g., <code>MYC (4609)</code>), and are harmonized offline via <code>data/gene_aliases.json</code>, containing 40,402 unambiguous HGNC previous-symbol and alias mappings (Seal <em>et al.</em>, 2023; e.g., <span class="gene">HER2</span>&nbsp;&rarr;&nbsp;<span class="gene">ERBB2</span>). Resolved aliases, duplicates, and unmatched identifiers are reported explicitly alongside each run.
  </p>

  <h3>2.2 Cohort-Scoped Essentiality and Lineage-Selectivity Estimators</h3>
  <p>
    Let <span class="math"><b>X</b> &isin; (&#8477; &cup; {{NA}})<sup><i>M</i>&times;<i>G</i></sup></span> denote the cell-line-by-gene Chronos matrix, and let <span class="math"><i>C</i> &sube; {{1, &hellip;, <i>M</i>}}</span> be the active cohort defined by the selected cancer lineages (and optional disease subtypes). For each gene <span class="math"><i>g</i></span> with finite observation set <span class="math"><i>C</i><sub><i>g</i></sub> = {{<i>i</i> &isin; <i>C</i> : <i>x</i><sub><i>ig</i></sub> &isin; &#8477;}}</span> and sample count <span class="math"><i>n</i><sub><i>g,C</i></sub> = |<i>C</i><sub><i>g</i></sub>| &gt; 0</span>, cohort-scoped essentiality is summarized by its arithmetic mean in <code>float64</code> precision and ascending empirical percentile rank across all <span class="math"><i>G</i><sub><i>C</i></sub></span> valid genes:
  </p>

  <div class="eq-row">
    <div class="eq-body">
      <span class="overline"><i>x</i></span><sub><i>g,C</i></sub> =
      <span class="frac"><span class="num">1</span><span class="den"><i>n</i><sub><i>g,C</i></sub></span></span>
      <span class="sum-op"><span class="sym">&sum;</span><span class="lim"><i>i</i> &isin; <i>C</i><sub><i>g</i></sub></span></span>
      <i>x</i><sub><i>ig</i></sub>,
      <span class="eq-gap"></span>
      <i>P</i><sub><i>g,C</i></sub> =
      <span class="frac"><span class="num">100 &middot; rank(<span class="overline"><i>x</i></span><sub><i>g,C</i></sub>)</span><span class="den"><i>G</i><sub><i>C</i></sub></span></span>
    </div>
    <div class="eq-num">(1)</div>
  </div>

  <p>
    Non-finite entries (<code>NaN</code>, &plusmn;&infin;) are excluded from sums and tracked separately in <code>n_missing</code> and <code>diagnostics.csv</code>. Because Chronos scores are scaled such that non-essential genes center near 0 and core essential genes center near &minus;1 (Dempster <em>et al.</em>, 2021), a dashed reference line at &minus;0.5 and two configurable reference anchors (default <span class="gene">MYC</span> and <span class="gene">PTEN</span>) orient multi-gene ranking curves and two-layer gene-set overlays without implying a binary essentiality classifier.
  </p>
  <p>
    To quantify lineage selectivity while diagnosing intra-lineage skew, the lineage-preference estimator compares each lineage <span class="math">&ell;</span> meeting <span class="math"><i>n</i><sub><i>g</i>,&ell;</sub> &ge; 5</span> against the full-dataset distribution of gene <span class="math"><i>g</i></span> (independent of active sidebar filters) using three complementary statistics:
  </p>

  <div class="eq-row">
    <div class="eq-body">
      &Delta;<i>&mu;</i><sub><i>g</i>,&ell;</sub> = <span class="overline"><i>x</i></span><sub><i>g</i>,&ell;</sub> &minus; <span class="overline"><i>x</i></span><sub><i>g</i>,all</sub>,
      <span class="eq-gap"></span>
      &Delta;med<sub><i>g</i>,&ell;</sub> = med<sub><i>g</i>,&ell;</sub> &minus; med<sub><i>g</i>,all</sub>,
      <span class="eq-gap"></span>
      <i>f</i><sub><i>g</i>,&ell;</sub> =
      <span class="frac"><span class="num">1</span><span class="den"><i>n</i><sub><i>g</i>,&ell;</sub></span></span>
      <span class="sum-op"><span class="sym">&sum;</span><span class="lim"><i>i</i> &isin; <i>C</i><sub><i>g</i>,&ell;</sub></span></span>
      <b>I</b>(<i>x</i><sub><i>ig</i></sub> &lt; &minus;0.5)
    </div>
    <div class="eq-num">(2)</div>
  </div>

  <p>
    where <span class="math">med</span> denotes the sample median and <span class="math"><b>I</b>(&middot;)</span> is the indicator function. Comparing <span class="math">&Delta;<i>&mu;</i><sub><i>g</i>,&ell;</sub></span> with <span class="math">&Delta;med<sub><i>g</i>,&ell;</sub></span> distinguishes broadly shared lineage dependencies (where mean and median shifts align) from subtype-restricted or bimodal vulnerabilities (where a sensitive subpopulation pulls <span class="math">&Delta;<i>&mu;</i><sub><i>g</i>,&ell;</sub></span> well below <span class="math">&Delta;med<sub><i>g</i>,&ell;</sub></span>). Paired cell-line tables surface the most and least dependent individual models (display name and <code>ModelID</code>) within <span class="math"><i>C</i></span> for direct experimental selection.
  </p>

  <h3>2.3 Pairwise-Complete Vectorized Co-Dependency Estimation</h3>
  <p>
    Given a query gene <span class="math"><i>q</i></span> and candidate partner <span class="math"><i>g</i> &ne; <i>q</i></span> within active cohort <span class="math"><i>C</i></span>, let <span class="math"><i>S</i><sub><i>q,g</i></sub> = {{<i>i</i> &isin; <i>C</i> : <i>x</i><sub><i>iq</i></sub>, <i>x</i><sub><i>ig</i></sub> &isin; &#8477;}}</span> be the shared finite cell-line index set with overlap size <span class="math"><i>n</i><sub><i>q,g</i></sub> = |<i>S</i><sub><i>q,g</i></sub>|</span>. Correlations are computed strictly over <span class="math"><i>S</i><sub><i>q,g</i></sub></span> without missing-value imputation whenever <span class="math"><i>n</i><sub><i>q,g</i></sub> &ge; <i>n</i><sub>min</sub> = min(10, max(3, <i>n</i><sub><i>q,C</i></sub>))</span> and both vectors have strictly positive sample variance over <span class="math"><i>S</i><sub><i>q,g</i></sub></span>:
  </p>

  <div class="eq-row">
    <div class="eq-body">
      <i>r</i><sub><i>q,g</i></sub>(<i>C</i>) =
      <span class="frac frac-tall">
        <span class="num">
          <span class="sum-op"><span class="sym">&sum;</span><span class="lim"><i>i</i> &isin; <i>S</i><sub><i>q,g</i></sub></span></span>
          (<i>x</i><sub><i>iq</i></sub> &minus; <span class="overline"><i>x</i></span><sub><i>q</i></sub><sup>(<i>S</i>)</sup>)(<i>x</i><sub><i>ig</i></sub> &minus; <span class="overline"><i>x</i></span><sub><i>g</i></sub><sup>(<i>S</i>)</sup>)
        </span>
        <span class="den">
          <span style="font-size:11.5pt; vertical-align:-0.08em;">[</span>&thinsp;<span class="sum-op"><span class="sym">&sum;</span><span class="lim"><i>i</i> &isin; <i>S</i><sub><i>q,g</i></sub></span></span>(<i>x</i><sub><i>iq</i></sub> &minus; <span class="overline"><i>x</i></span><sub><i>q</i></sub><sup>(<i>S</i>)</sup>)<sup>2</sup>
          &thinsp;&middot;&thinsp;
          <span class="sum-op"><span class="sym">&sum;</span><span class="lim"><i>i</i> &isin; <i>S</i><sub><i>q,g</i></sub></span></span>(<i>x</i><sub><i>ig</i></sub> &minus; <span class="overline"><i>x</i></span><sub><i>g</i></sub><sup>(<i>S</i>)</sup>)<sup>2</sup>&thinsp;<span style="font-size:11.5pt; vertical-align:-0.08em;">]</span><sup>1/2</sup>
        </span>
      </span>
    </div>
    <div class="eq-num">(3)</div>
  </div>

  <p>
    where <span class="math"><span class="overline"><i>x</i></span><sub><i>q</i></sub><sup>(<i>S</i>)</sup></span> and <span class="math"><span class="overline"><i>x</i></span><sub><i>g</i></sub><sup>(<i>S</i>)</sup></span> denote means computed over <span class="math"><i>S</i><sub><i>q,g</i></sub></span>. Rank correlation <span class="math"><i>&rho;</i><sub><i>q,g</i></sub>(<i>C</i>)</span> applies the same formula after replacing <span class="math"><i>x</i><sub>&middot;<i>q</i></sub></span> and <span class="math"><i>x</i><sub>&middot;<i>g</i></sub></span> with average-tie midranks on <span class="math"><i>S</i><sub><i>q,g</i></sub></span>.
  </p>
  <p>
    To achieve interactive latency while keeping peak memory bounded and avoiding SciPy runtime dependencies, <code>genomewide_codependency</code> processes candidate genes in 2,048-column <code>float64</code> blocks after restricting rows to finite observations of <span class="math"><i>q</i></span>. Within each block, columns complete over <span class="math"><i>C</i><sub><i>q</i></sub></span> are evaluated via a single centered BLAS matrix–vector inner product, whereas columns with missing entries use masked pairwise means and variances for Pearson <em>r</em> and per-column exact pairwise midrank evaluation for Spearman <em>&rho;</em>. Time complexity is <span class="math"><i>O</i>(|<i>C</i>|<i>G</i>)</span> for Pearson scanning with <span class="math"><i>O</i>(2048&nbsp;&middot;&nbsp;|<i>C</i>|)</span> temporary memory. For multi-gene panels (2–25 genes), the same pairwise estimator populates a symmetric correlation matrix linked to lineage-colored cell-line scatter plots reporting both Pearson <em>r</em> and Spearman <em>&rho;</em>.
  </p>

  <h3>2.4 Software Architecture, Provenance, and Benchmark Protocol</h3>
  <p>
    All analytical kernels are packaged in a Streamlit-free PEP 621 library (<code>crispr_score_analyzer</code>, defined in <code>pyproject.toml</code>, requiring only NumPy, pandas, Plotly, Kaleido, and Pillow) organized into four modules: <code>analysis_core</code>, <code>codependency</code>, <code>plotting</code>, and <code>figure_export</code>. This shared core drives both a scriptable CLI (<code>crispr-score-analyzer</code> exposing <code>info</code>, <code>rank</code>, <code>codep</code>, and <code>app</code> subcommands; documented in <code>docs/api.md</code>) and a bilingual (English/Chinese) Streamlit web interface (<strong>Figure 1</strong>; <strong>Table 1</strong>). All interactive figures and white-background PDF, SVG, and 300-DPI PNG exports use colorblind-safe Nature Publishing Group (NPG) categorical hues and a 7-stop ColorBrewer <code>RdBu</code> diverging scale; interactive scatter curves use WebGL (<code>Scattergl</code>) and convert traces to vector paths on export. To prevent stale-state inconsistencies, input edits remain isolated until executed and cohort changes invalidate downstream tables and exports. Each session can export a self-contained provenance ZIP containing <code>selected_results.csv</code>, <code>full_rankings.csv</code>, <code>plot_data.csv</code>, <code>lineage_preference.csv</code>, <code>diagnostics.csv</code>, and <code>metadata.json</code> (data SHA-256, software version, Zenodo DOI, cohort scope, alias resolution log, and package versions). Software correctness is verified offline on deterministic synthetic fixtures by 98 automated unit, API, CLI, UI, and export tests across Python 3.11 and 3.12 via GitHub Actions.
  </p>

  <!-- Figure 1 -->
  <div class="figure-block fig1">
    <img src="{fig1_uri}" alt="Figure 1" />
    <p class="figure-caption">
      <strong>Figure 1. Browser interface of CRISPR Score Analyzer operating on the DepMap Public 26Q3 pan-cancer cohort.</strong>
      <strong>(A)</strong> Light display mode showing sidebar dataset and lineage-scope controls, four-module navigation (<em>Gene ranking</em>, <em>Cancer-type comparison</em>, <em>Gene-set annotation</em>, <em>Co-dependency</em>), batch gene-list input with alias diagnostics, and the pan-cancer ranking curve for five benchmark genes (y-axis labeled <em>Mean CRISPR score</em> in the generic interface to accommodate custom score matrices).
      <strong>(B)</strong> Dark display mode aligned on the identical cohort view; all downloadable PDF, SVG, and PNG exports standardize to a white background.
    </p>
  </div>

  <!-- Table 1 -->
  <div class="table-block">
    <p class="table-caption">
      <strong>Table 1.</strong> Comparison of analytical and reproducibility capabilities across CRISPR Score Analyzer and existing DepMap interfaces (verified November 2026).
    </p>
    <table class="cmp-table">
      <thead>
        <tr>
          <th>Analytical capability</th>
          <th>DepMap Portal (26Q3; depmap.org/portal)</th>
          <th>shinyDepMap (eLife 2021; shinyapps.io/depmap)</th>
          <th class="col-csa">CRISPR Score Analyzer (<code>v1.2.0</code> / <code>preprint-v4</code>)</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td class="cap-col">Multi-gene panel placement on cohort curve</td>
          <td>Single-gene default (<code>Gene Page</code>)</td>
          <td>Efficacy vs. selectivity scatter (<code>Efficacy &amp; Selectivity</code>)</td>
          <td class="csa-col">Batch ranking and two-layer gene-set annotation</td>
        </tr>
        <tr>
          <td class="cap-col">Cohort-scoped genome-wide co-dependency scan</td>
          <td>Pan-cancer default; subsets via multi-step <code>Custom Analysis</code></td>
          <td>Pan-cancer clusters (<code>Connectivity</code>)</td>
          <td class="csa-col">Single-step scoped scan (exact pairwise Pearson <em>r</em> &amp; Spearman <em>&rho;</em>)</td>
        </tr>
        <tr>
          <td class="cap-col">Pairwise multi-gene matrix &amp; lineage scatter</td>
          <td><code>Data Explorer 2</code> / <code>Custom Analysis</code> small multiples</td>
          <td>Cluster network view</td>
          <td class="csa-col">2–25 gene heatmap linked to lineage-colored scatter</td>
        </tr>
        <tr>
          <td class="cap-col">Lineage selectivity &amp; intra-lineage skew</td>
          <td>Per-gene lineage &amp; subtype plots</td>
          <td>Not provided</td>
          <td class="csa-col">Joint &Delta;<em>&mu;</em>, &Delta;med, and fraction &lt; &minus;0.5 lollipop/heatmap</td>
        </tr>
        <tr>
          <td class="cap-col">Sensitive &amp; resistant cell-line selection</td>
          <td>Table filtering on <code>Gene Page</code></td>
          <td>Top sensitive cell lines</td>
          <td class="csa-col">Paired most/least dependent tables with <code>ModelID</code> &amp; CSV export</td>
        </tr>
        <tr>
          <td class="cap-col">Custom matrix ingestion &amp; schema audit</td>
          <td>Single-feature custom vector</td>
          <td>Not supported</td>
          <td class="csa-col">Full UTF-8 CSV matrix, strict schema guard &amp; <code>diagnostics.csv</code></td>
        </tr>
        <tr>
          <td class="cap-col">Offline batch HGNC harmonization</td>
          <td>Online single-gene search</td>
          <td>Not provided</td>
          <td class="csa-col">40,402 unambiguous offline mappings (batch &amp; single)</td>
        </tr>
        <tr>
          <td class="cap-col">Scriptable Python API, CLI &amp; provenance</td>
          <td>Web portal / manual bulk downloads</td>
          <td>R Shiny web application only</td>
          <td class="csa-col"><code>pip</code>-installable Python library, CLI, vector export &amp; SHA-256 ZIP</td>
        </tr>
      </tbody>
    </table>
  </div>

  <p>
    All quantitative metrics reported in Section&nbsp;3 are reproduced offline by <code>preprint/reproduce_stats.py</code> without SciPy. For post-hoc evaluation of single-lineage co-dependency scans, two-sided Pearson <em>p</em>-values are computed from the exact Student-<em>t</em> distribution (<span class="math">df = <i>n</i><sub><i>q,g</i></sub> &minus; 2</span>) via continued-fraction evaluation of the regularized incomplete beta function, accompanied by Bonferroni and Benjamini–Hochberg false discovery rate (FDR) adjustments and ordinary least-squares (OLS) partial correlations on complete observations.
  </p>

  <!-- Section 3 -->
  <h2><span class="sec-num">3.</span>Results</h2>

  <h3>3.1 Coverage Heterogeneity, Numerical Accuracy, and Runtime Scaling in DepMap Public 26Q3</h3>
  <p>
    Auditing the bundled DepMap Public 26Q3 matrix (1,361 cell lines &times; 18,434 genes; 1,712 genes with pan-cancer mean <span class="nowrap">&lt;&nbsp;&minus;0.5</span>, overall mean range &minus;5.250 to +0.398) reveals substantial gene-level missingness resulting from multi-library integration. Whereas 17,082 genes (92.67%) are profiled across all 1,361 cell lines, 1,352 genes (7.33%) have incomplete coverage ranging down to <span class="nowrap"><code>n_valid</code>&nbsp;=&nbsp;218</span>, including 648 genes assayed in fewer than 1,000 cell lines and 519 in fewer than 500. In correlation scans over such heterogeneous matrices, failing to preserve exact pairwise overlap or omitting <code>n_shared</code> from outputs risks equating sparse sub-cohort correlations with full-panel associations.
  </p>
  <p>
    We verified numerical equivalence of the chunked <code>float64</code> correlation engine against direct reference evaluations (<code>pandas.Series.corr</code> for Pearson <em>r</em> and pairwise midrank <code>numpy.corrcoef</code> for Spearman <em>&rho;</em>) across top, bottom, and partially observed partners of <span class="gene">MDM2</span> (193 genes for Pearson, 100 for Spearman). Maximum absolute deviations were <span class="nowrap">7.77&nbsp;&times;&nbsp;10<sup>&minus;16</sup></span> for Pearson <em>r</em> and <span class="nowrap">5.55&nbsp;&times;&nbsp;10<sup>&minus;17</sup></span> for Spearman <em>&rho;</em>, confirming machine-precision fidelity across both complete and incomplete gene columns. On a 2-vCPU Linux environment (Python 3.12.16, NumPy 2.5.5, pandas 2.3.3), median wall-clock runtimes across three replicates were 1.01&nbsp;s (0.97–1.33&nbsp;s) to load the 88&nbsp;MB matrix and 0.41&nbsp;s (0.38–0.73&nbsp;s) to rank all 18,434 genes. A genome-wide pan-cancer scan (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;1,361</span>, 18,433 partner genes) completed in 0.74&nbsp;s (0.73–0.85&nbsp;s) for Pearson <em>r</em> and 6.19&nbsp;s (6.02–6.26&nbsp;s) for Spearman <em>&rho;</em>; restricting to a single lineage (Skin, <span class="nowrap"><em>n</em>&nbsp;=&nbsp;75</span>) reduced scan times to 0.10&nbsp;s (Pearson) and 0.44&nbsp;s (Spearman), supporting interactive exploration across both metrics.
  </p>

  <h3>3.2 Lineage-Selectivity Decomposition and Subtype-Driven Skew</h3>
  <p>
    Placing a five-gene benchmark panel (<span class="gene">KRAS</span>, <span class="gene">SOX10</span>, <span class="gene">BRAF</span>, <span class="gene">ESR1</span>, <span class="gene">TP53</span>) on the pan-cancer Chronos ranking curve (<strong>Figure 2A</strong>) illustrates how global averages compress context-restricted vulnerabilities. Relative to reference anchors <span class="gene">MYC</span> (mean &minus;1.971, rank 187, 1.01%) and <span class="gene">PTEN</span> (mean +0.377, rank 18,428, 99.97%), <span class="gene">KRAS</span> falls below the &minus;0.5 screening reference across all lineages (mean &minus;0.729, rank 1,263, 6.85%), whereas <span class="gene">TP53</span> lies at the positive extreme (mean +0.381, rank 18,430, 99.98%), reflecting proliferative advantage upon <span class="gene">TP53</span> knockout in <span class="gene">TP53</span>-wild-type lines (Tsherniak <em>et al.</em>, 2017; Dempster <em>et al.</em>, 2021). By contrast, two canonical lineage-dependent targets appear modest or near-neutral in pan-cancer means: <span class="gene">SOX10</span> (mean &minus;0.347, rank 2,225, 12.07%) and <span class="gene">ESR1</span> (mean &minus;0.065, rank 7,347, 39.86%). Two-layer gene-set annotation similarly positions mitotic kinases <span class="gene">PLK1</span> (mean &minus;2.691, rank 38, 0.21%) and <span class="gene">AURKA</span> (mean &minus;1.236, rank 626, 3.40%) against seven cell-cycle CDK/cyclin genes (<span class="gene">CDK1</span> &minus;2.588 through <span class="gene">CCNE1</span> &minus;0.194; <strong>Figure 2B</strong>).
  </p>
  <p>
    Evaluating &Delta;<em>&mu;</em> jointly with &Delta;med and <span class="math"><i>f</i><sub><i>g</i>,&ell;</sub></span> (<strong>Figure 2C, D</strong>) separates uniform lineage essentiality from subtype-restricted bimodality. For <span class="gene">SOX10</span>, Skin (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;75</span>) exhibits concordant mean and median shifts (<span class="nowrap">mean&nbsp;=&nbsp;&minus;1.487</span>, <span class="nowrap">&Delta;<em>&mu;</em>&nbsp;=&nbsp;&minus;1.140</span>, <span class="nowrap">&Delta;med&nbsp;=&nbsp;&minus;1.368</span>, <span class="nowrap">78.7%&nbsp;&lt;&nbsp;&minus;0.5</span>), consistent with broad neural-crest dependency in cutaneous melanoma. In Eye (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;15</span>), however, the mean shift (<span class="nowrap">&Delta;<em>&mu;</em>&nbsp;=&nbsp;&minus;0.535</span>, <span class="nowrap">mean&nbsp;=&nbsp;&minus;0.882</span>, <span class="nowrap">46.7%&nbsp;&lt;&nbsp;&minus;0.5</span>) is more than threefold stronger than the median shift (<span class="nowrap">&Delta;med&nbsp;=&nbsp;&minus;0.174</span>), flagging a bimodal distribution driven by uveal melanoma models (<code>Omm2.5</code> [<code>ACH-002018</code>, &minus;2.752], <code>UPMM3</code> [<code>ACH-002925</code>, &minus;2.113], <code>WM3772F</code> [<code>ACH-002687</code>, &minus;1.678]) mixed with near-neutral retinoblastoma lines. A comparable subtype-masking effect appears for <span class="gene">ESR1</span> in Breast (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;67</span>, <span class="nowrap">mean&nbsp;=&nbsp;&minus;0.285</span>, <span class="nowrap">&Delta;<em>&mu;</em>&nbsp;=&nbsp;&minus;0.220</span>, <span class="nowrap">&Delta;med&nbsp;=&nbsp;&minus;0.061</span>, <span class="nowrap">17.9%&nbsp;&lt;&nbsp;&minus;0.5</span>), where Oncotree lineage grouping pools ER-positive luminal models with ER-negative subtypes and produces a lineage shift only slightly larger than a small Ampulla of Vater group (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;9</span>, <span class="nowrap">&Delta;<em>&mu;</em>&nbsp;=&nbsp;&minus;0.174</span>, <span class="nowrap">&Delta;med&nbsp;=&nbsp;&minus;0.059</span>). Inspecting the ranked cell-line table immediately resolves this aggregation by identifying strongly dependent ER-positive luminal models (<code>HCC1428</code> [<code>ACH-000352</code>, &minus;1.912], <code>T-47D</code> [<code>ACH-000147</code>, &minus;1.845], <code>MDA-MB-415</code> [<code>ACH-000876</code>, &minus;1.636]) alongside resistant negative controls.
  </p>

  <!-- Figure 2 -->
  <div class="figure-block fig2">
    <img src="{fig2_uri}" alt="Figure 2" />
    <p class="figure-caption">
      <strong>Figure 2. Multi-gene essentiality positioning and lineage-selectivity decomposition in DepMap Public 26Q3 (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;1,361</span> cell lines, 18,434 genes).</strong>
      <strong>(A)</strong> Pan-cancer mean Chronos gene effect curve positioning five benchmark genes (<span class="gene">KRAS</span>, <span class="gene">SOX10</span>, <span class="gene">BRAF</span>, <span class="gene">ESR1</span>, <span class="gene">TP53</span>) relative to <span class="gene">MYC</span>/<span class="gene">PTEN</span> anchors and the &minus;0.5 reference line.
      <strong>(B)</strong> Two-layer gene-set annotation overlaying <span class="gene">PLK1</span> and <span class="gene">AURKA</span> onto seven cell-cycle CDK/cyclin genes (<span class="gene">CDK1</span>, <span class="gene">CDK2</span>, <span class="gene">CDK4</span>, <span class="gene">CDK6</span>, <span class="gene">CCNB1</span>, <span class="gene">CCND1</span>, <span class="gene">CCNE1</span>).
      <strong>(C, D)</strong> Lineage-preference lollipop plots (top and bottom 4 lineages with <span class="nowrap"><em>n</em>&nbsp;&ge;&nbsp;5</span>; gene-specific x-axis scales) displaying mean Chronos shifts &Delta;<em>&mu;</em> relative to each gene’s all-lineage mean for <span class="gene">SOX10</span> (<strong>C</strong>) and <span class="gene">ESR1</span> (<strong>D</strong>).
    </p>
  </div>

  <h3>3.3 Pan-Cancer Versus Cohort-Scoped Co-Dependency Across Canonical Cancer Axes</h3>
  <p>
    To quantify how cohort scope alters co-dependency inference, we evaluated ten established cancer gene pairs across three estimation regimes in DepMap Public 26Q3: pan-cancer Pearson <em>r</em> and Spearman <em>&rho;</em> (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;1,361</span>), pooled within-lineage Pearson correlation after subtracting lineage means (<span class="math"><i>r</i><sub>within</sub></span> across all 30 lineages), and single-lineage Pearson <em>r</em> and Spearman <em>&rho;</em> within each pair’s primary biological context (<strong>Table 2</strong>).
  </p>

  <!-- Table 2 -->
  <div class="table-block">
    <p class="table-caption">
      <strong>Table 2.</strong> Comparison of pan-cancer, lineage-demeaned, and target-lineage co-dependency across ten canonical cancer gene pairs in DepMap Public 26Q3 (all pan-cancer pairs <span class="nowrap"><em>n</em>&nbsp;=&nbsp;1,361</span>).
    </p>
    <table class="cmp-table">
      <thead>
        <tr>
          <th>Gene pair</th>
          <th>Target lineage (<em>n</em>)</th>
          <th class="num-col">Pan <em>r</em></th>
          <th class="num-col">Pan <em>&rho;</em></th>
          <th class="num-col"><em>r</em><sub>within</sub></th>
          <th class="num-col col-csa">Lineage <em>r</em></th>
          <th class="num-col col-csa">Lineage <em>&rho;</em></th>
          <th>Empirical regime</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td><span class="gene">MDM2</span>–<span class="gene">TP53</span></td>
          <td>Lung (126)</td>
          <td class="num-col">&minus;0.702</td>
          <td class="num-col">&minus;0.704</td>
          <td class="num-col">&minus;0.648</td>
          <td class="num-col"><strong>&minus;0.802</strong></td>
          <td class="num-col">&minus;0.500</td>
          <td>Pan-lineage epistatic</td>
        </tr>
        <tr>
          <td><span class="gene">MDM2</span>–<span class="gene">PPM1D</span></td>
          <td>Lung (126)</td>
          <td class="num-col">+0.664</td>
          <td class="num-col">+0.658</td>
          <td class="num-col">+0.635</td>
          <td class="num-col"><strong>+0.642</strong></td>
          <td class="num-col">+0.584</td>
          <td>Pan-lineage epistatic</td>
        </tr>
        <tr>
          <td><span class="gene">BRAF</span>–<span class="gene">MAPK1</span></td>
          <td>Skin (75)</td>
          <td class="num-col">+0.644</td>
          <td class="num-col">+0.451</td>
          <td class="num-col">+0.538</td>
          <td class="num-col"><strong>+0.798</strong></td>
          <td class="num-col"><strong>+0.787</strong></td>
          <td>Lineage-enriched pathway</td>
        </tr>
        <tr>
          <td><span class="gene">BRAF</span>–<span class="gene">MAP2K1</span></td>
          <td>Skin (75)</td>
          <td class="num-col">+0.391</td>
          <td class="num-col">+0.097</td>
          <td class="num-col">+0.303</td>
          <td class="num-col"><strong>+0.552</strong></td>
          <td class="num-col"><strong>+0.605</strong></td>
          <td>Lineage-enriched pathway</td>
        </tr>
        <tr>
          <td><span class="gene">BRAF</span>–<span class="gene">SOX10</span></td>
          <td>Skin (75)</td>
          <td class="num-col">+0.544</td>
          <td class="num-col">+0.196</td>
          <td class="num-col">+0.294</td>
          <td class="num-col">+0.633</td>
          <td class="num-col">+0.643</td>
          <td>Lineage co-enrichment / subtype</td>
        </tr>
        <tr>
          <td><span class="gene">SOX10</span>–<span class="gene">MITF</span></td>
          <td>Skin (75)</td>
          <td class="num-col">+0.492</td>
          <td class="num-col">+0.153</td>
          <td class="num-col">+0.286</td>
          <td class="num-col">+0.418</td>
          <td class="num-col">+0.384</td>
          <td>Lineage-restricted TF axis</td>
        </tr>
        <tr>
          <td><span class="gene">ESR1</span>–<span class="gene">FOXA1</span></td>
          <td>Breast (67)</td>
          <td class="num-col">+0.324</td>
          <td class="num-col">+0.111</td>
          <td class="num-col">+0.259</td>
          <td class="num-col"><strong>+0.503</strong></td>
          <td class="num-col">+0.346</td>
          <td>Subtype-restricted TF axis</td>
        </tr>
        <tr>
          <td><span class="gene">ESR1</span>–<span class="gene">GATA3</span></td>
          <td>Breast (67)</td>
          <td class="num-col">+0.253</td>
          <td class="num-col">+0.095</td>
          <td class="num-col">+0.228</td>
          <td class="num-col"><strong>+0.596</strong></td>
          <td class="num-col">+0.372</td>
          <td>Subtype-restricted TF axis</td>
        </tr>
        <tr>
          <td><span class="gene">CDK4</span>–<span class="gene">CCND1</span></td>
          <td>Breast (67)</td>
          <td class="num-col">+0.416</td>
          <td class="num-col">+0.380</td>
          <td class="num-col">+0.425</td>
          <td class="num-col"><strong>+0.704</strong></td>
          <td class="num-col"><strong>+0.635</strong></td>
          <td>Context-amplified cell cycle</td>
        </tr>
        <tr>
          <td><span class="gene">KRAS</span>–<span class="gene">RAF1</span></td>
          <td>Bowel (85)</td>
          <td class="num-col">+0.338</td>
          <td class="num-col">+0.255</td>
          <td class="num-col">+0.319</td>
          <td class="num-col"><strong>+0.530</strong></td>
          <td class="num-col">+0.461</td>
          <td>Lineage-amplified effector</td>
        </tr>
      </tbody>
    </table>
  </div>

  <p>
    The benchmark separates three distinct statistical behaviors. First, pan-lineage epistatic pairs (<span class="gene">MDM2</span>–<span class="gene">TP53</span>, <span class="gene">MDM2</span>–<span class="gene">PPM1D</span>) exhibit tight agreement among pan-cancer <em>r</em>, pan-cancer <em>&rho;</em>, and lineage-demeaned <span class="math"><i>r</i><sub>within</sub></span> (<span class="gene">MDM2</span>–<span class="gene">TP53</span>: &minus;0.702, &minus;0.704, &minus;0.648). In a genome-wide pan-cancer scan of <span class="gene">MDM2</span> across 18,433 genes (<strong>Figure 3A</strong>), top positive partners include <span class="gene">PPM1D</span> (+0.664), <span class="gene">MDM4</span> (+0.614), <span class="gene">USP7</span> (+0.522), <span class="gene">DDX31</span> (+0.466), and <span class="gene">UBE2D3</span> (+0.398, all <span class="nowrap"><code>n_shared</code>&nbsp;=&nbsp;1,361</span>), whereas top anti-correlated partners recover the p53 tumor-suppressor axis (<span class="gene">TP53</span> &minus;0.702, <span class="gene">TP53BP1</span> &minus;0.589, <span class="gene">USP28</span> &minus;0.531, <span class="gene">CDKN1A</span> &minus;0.530, <span class="gene">CHEK2</span> &minus;0.491, <span class="gene">ATM</span> &minus;0.461, all <span class="nowrap"><code>n_shared</code>&nbsp;=&nbsp;1,361</span>). Reporting <code>n_shared</code> directly on the figure flags rank-6 positive hit <span class="gene">TERF1</span> (+0.395) as partially covered (<span class="nowrap"><code>n_shared</code>&nbsp;=&nbsp;228</span>), distinguishing it from full-panel associations. Multi-gene matrix inspection of the six-gene regulatory axis resolves its two-block inverse structure (<strong>Figure 3B</strong>), and lineage-stratified analysis confirms negative <span class="gene">MDM2</span>–<span class="gene">TP53</span> correlation across all 23/23 lineages with <span class="nowrap"><em>n</em>&nbsp;&ge;&nbsp;15</span> (median <span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.680</span>; e.g., Ovary/Fallopian Tube &minus;0.818 [<span class="nowrap"><em>n</em>&nbsp;=&nbsp;62</span>], Lung &minus;0.802 [<span class="nowrap"><em>n</em>&nbsp;=&nbsp;126</span>], Skin &minus;0.793 [<span class="nowrap"><em>n</em>&nbsp;=&nbsp;75</span>], Breast &minus;0.719 [<span class="nowrap"><em>n</em>&nbsp;=&nbsp;67</span>], Bowel &minus;0.674 [<span class="nowrap"><em>n</em>&nbsp;=&nbsp;85</span>], CNS/Brain &minus;0.663 [<span class="nowrap"><em>n</em>&nbsp;=&nbsp;133</span>]; <strong>Figure 3C</strong>).
  </p>
  <p>
    Second, context-restricted signaling and transcriptional modules suffer strong attenuation in pan-cancer screens, reflected directly in divergence between pan-cancer Pearson <em>r</em> and Spearman <em>&rho;</em>. When a pathway operates in a small fraction of cell lines, extreme scores in the active lineage generate moderate Pearson covariance while &gt;90% of non-expressing lines contribute rank noise near zero, collapsing pan-cancer <em>&rho;</em> toward zero (<span class="gene">BRAF</span>–<span class="gene">MAP2K1</span>: <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.391</span> vs. <span class="nowrap"><em>&rho;</em>&nbsp;=&nbsp;+0.097</span>; <span class="gene">ESR1</span>–<span class="gene">GATA3</span>: <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.253</span> vs. <span class="nowrap"><em>&rho;</em>&nbsp;=&nbsp;+0.095</span>; <span class="gene">ESR1</span>–<span class="gene">FOXA1</span>: <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.324</span> vs. <span class="nowrap"><em>&rho;</em>&nbsp;=&nbsp;+0.111</span>). Restricting the cohort to the relevant lineage restores parametric and rank concordance simultaneously (<span class="gene">BRAF</span>–<span class="gene">MAP2K1</span> in Skin: <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.552</span>, <span class="nowrap"><em>&rho;</em>&nbsp;=&nbsp;+0.605</span>; <span class="gene">ESR1</span>–<span class="gene">GATA3</span> in Breast: <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.596</span>, <span class="nowrap"><em>&rho;</em>&nbsp;=&nbsp;+0.372</span>; <span class="gene">CDK4</span>–<span class="gene">CCND1</span> in Breast: <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.704</span>, <span class="nowrap"><em>&rho;</em>&nbsp;=&nbsp;+0.635</span>; <span class="gene">KRAS</span>–<span class="gene">RAF1</span> in Bowel: <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.530</span>, <span class="nowrap"><em>&rho;</em>&nbsp;=&nbsp;+0.461</span>).
  </p>

  <!-- Figure 3 -->
  <div class="figure-block fig3">
    <img src="{fig3_uri}" alt="Figure 3" />
    <p class="figure-caption">
      <strong>Figure 3. Pan-cancer and lineage-scoped genome-wide co-dependency benchmarks in DepMap Public 26Q3.</strong>
      <strong>(A)</strong> Pan-cancer Pearson <em>r</em> scan of <span class="gene">MDM2</span> across 18,433 partner genes (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;1,361</span>), showing top 7 positive (<span class="gene">PPM1D</span>, <span class="gene">MDM4</span>, <span class="gene">USP7</span>, <span class="gene">DDX31</span>, <span class="gene">UBE2D3</span>, <span class="gene">TERF1</span> [<span class="nowrap"><em>n</em>&nbsp;=&nbsp;228</span>], <span class="gene">PPM1G</span>) and top 7 negative (<span class="gene">TP53</span>, <span class="gene">TP53BP1</span>, <span class="gene">USP28</span>, <span class="gene">CDKN1A</span>, <span class="gene">CHEK2</span>, <span class="gene">ATM</span>, <span class="gene">UBE2K</span>) partners.
      <strong>(B)</strong> Pairwise Pearson correlation heatmap for the six-gene <span class="gene">TP53</span>–<span class="gene">MDM2</span> regulatory panel.
      <strong>(C)</strong> Lineage-colored cell-line scatter plot of <span class="gene">MDM2</span> vs. <span class="gene">TP53</span> (<span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.702</span>, Spearman <span class="nowrap"><em>&rho;</em>&nbsp;=&nbsp;&minus;0.704</span>, <span class="nowrap"><em>n</em>&nbsp;=&nbsp;1,361</span>) with OLS trendline.
      <strong>(D)</strong> Single-lineage genome-wide <span class="gene">BRAF</span> scan restricted to Skin models (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;75</span>, 17,785 partners with <span class="nowrap"><code>n_shared</code>&nbsp;&ge;&nbsp;10</span>), showing elevated MAPK cascade correlations (<span class="gene">MAPK1</span> +0.80, <span class="gene">MAP2K1</span> +0.55, <span class="gene">MAP2K2</span> +0.55) alongside melanoma lineage factors (<span class="gene">SOX10</span> +0.63, <span class="gene">SOX5</span> +0.56) and exploratory hits near <span class="nowrap">|<em>r</em>|&nbsp;&approx;&nbsp;0.50</span>.
    </p>
  </div>

  <p>
    Third, between-lineage baseline co-enrichment inflates pan-cancer Pearson correlations for lineage-restricted transcription factors (<span class="gene">SOX10</span>–<span class="gene">MITF</span>: pan-cancer <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.492</span> drops to <span class="nowrap"><span class="math"><i>r</i><sub>within</sub></span>&nbsp;=&nbsp;+0.286</span> after lineage demeaning and <span class="nowrap"><em>&rho;</em>&nbsp;=&nbsp;+0.153</span>; <span class="gene">BRAF</span>–<span class="gene">SOX10</span>: +0.544 drops to <span class="nowrap"><span class="math"><i>r</i><sub>within</sub></span>&nbsp;=&nbsp;+0.294</span>). Single-lineage scoping, however, does not automatically remove intra-lineage subtype confounding. In the Skin-restricted <span class="gene">BRAF</span> scan (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;75</span>, 17,785 partner genes meeting <span class="nowrap"><code>n_shared</code>&nbsp;&ge;&nbsp;10</span>, 648 sparser genes excluded; <strong>Figure 3D</strong>), MAPK effectors strengthen relative to pan-cancer values (<span class="gene">MAPK1</span> +0.644&nbsp;&rarr;&nbsp;+0.798, <span class="gene">MAP2K1</span> +0.391&nbsp;&rarr;&nbsp;+0.552, <span class="gene">MAP2K2</span> +0.271&nbsp;&rarr;&nbsp;+0.550), while <span class="gene">SOX10</span> (+0.633) and <span class="gene">SOX5</span> (+0.562) remain among the top hits alongside <span class="gene">UBQLN1</span> (+0.526), <span class="gene">AKIRIN1</span> (+0.496), and negatively correlated <span class="gene">STYXL1</span> (&minus;0.534, all <span class="nowrap"><code>n_shared</code>&nbsp;=&nbsp;75</span>). Regressing both <span class="gene">BRAF</span> and each partner on <span class="gene">SOX10</span> via OLS across all 75 complete Skin models retains MAPK cascade correlations (<span class="gene">MAPK1</span> partial <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.665</span>, <span class="gene">MAP2K1</span> +0.467, <span class="gene">MAP2K2</span> +0.398), indicating that <span class="gene">SOX10</span> co-ranking arises from Oncotree Skin heterogeneity spanning <span class="gene">SOX10</span>-high cutaneous melanoma and <span class="gene">SOX10</span>-low non-melanoma or dedifferentiated lines (e.g., <code>HSC-5</code>, <code>MCC26</code>, <code>RVH421_DAB_R</code>). Post-hoc two-sided Pearson <em>t</em>-tests across all 17,785 Skin tests further separate <span class="gene">MAPK1</span> (<span class="nowrap"><em>p</em>&nbsp;=&nbsp;9.72&nbsp;&times;&nbsp;10<sup>&minus;18</sup></span>, Bonferroni <span class="nowrap">1.73&nbsp;&times;&nbsp;10<sup>&minus;13</sup></span>) from candidates near <span class="nowrap">|<em>r</em>|&nbsp;&approx;&nbsp;0.50</span> such as <span class="gene">STYXL1</span> (<span class="nowrap"><em>p</em>&nbsp;=&nbsp;8.28&nbsp;&times;&nbsp;10<sup>&minus;7</sup></span>, Bonferroni 0.015, BH FDR 0.0025), <span class="gene">UBQLN1</span> (<span class="nowrap"><em>p</em>&nbsp;=&nbsp;1.29&nbsp;&times;&nbsp;10<sup>&minus;6</sup></span>, Bonferroni 0.023, BH FDR 0.0033), and <span class="gene">AKIRIN1</span> (<span class="nowrap"><em>p</em>&nbsp;=&nbsp;5.95&nbsp;&times;&nbsp;10<sup>&minus;6</sup></span>, Bonferroni 0.106, BH FDR 0.0132), showing why moderate single-lineage correlations warrant independent validation rather than immediate mechanistic interpretation.
  </p>

  <!-- Section 4 -->
  <h2><span class="sec-num">4.</span>Discussion</h2>
  <p>
    Co-essentiality analyses of large-scale CRISPR screens balance global statistical normalization against local interpretability. Covariance-whitening methods such as generalized least squares (Wainberg <em>et al.</em>, 2021), olfactory-receptor principal-component removal (Boyle <em>et al.</em>, 2018), and PCA-whitened network pipelines (Gheorghe &amp; Hart, 2022) effectively suppress shared technical variance and dominant mitochondrial or proliferation axes when constructing genome-wide functional networks. Global linear whitening, however, alters the native Chronos scale and can obscure whether a correlation arises across tissues or within a small subset of cell lines. Cohort-scoped descriptive analysis serves a complementary role: operating directly on unwhitened Chronos scores inside an explicit lineage or disease subset, reporting both Pearson <em>r</em> and Spearman <em>&rho;</em> alongside exact <code>n_shared</code>, and linking aggregate correlations to individual cell-line scatter plots lets investigators inspect biological context before proceeding to network modeling or wet-lab validation.
  </p>
  <p>
    The DepMap Public 26Q3 benchmarks also delineate clear boundaries on interpreting descriptive outputs. Cohort ranks, lineage mean and median shifts, and pairwise correlations summarize empirical screen behavior rather than formal tests of differential dependency, selectivity, or synthetic lethality, and the &minus;0.5 line serves solely as a visual screening reference rather than a common-essential classifier. As seen for <span class="gene">ESR1</span> in Breast, <span class="gene">SOX10</span> in Eye, and <span class="gene">BRAF</span>–<span class="gene">SOX10</span> in Skin, Oncotree lineage labels pool distinct molecular subtypes (luminal versus basal breast models; cutaneous melanoma versus non-melanoma skin models), while scanning ~18,000 genes across <span class="nowrap"><em>n</em>&nbsp;&approx;&nbsp;75</span> cell lines introduces substantial multiple-testing variance around <span class="nowrap">|<em>r</em>|&nbsp;&approx;&nbsp;0.50</span>. Moreover, the bundled 26Q3 matrix intentionally stores only <code>OncotreeLineage</code> and cell-line identifiers without somatic mutation, transcriptomic, or receptor/PAM50 covariates, and uploaded custom tables are assumed to contain one row per unique cell line without silent row deduplication.
  </p>
  <p>
    Separating the analytical library (<code>crispr_score_analyzer</code>) and CLI from the Streamlit interface mitigates these constraints in extended workflows. Investigators who begin interactively to screen a candidate panel or choose sensitive and resistant cell lines can move directly to scripted Python pipelines backed by identical ranking and correlation routines, merge external multi-omic covariates (such as <span class="gene">TP53</span> mutation status, <span class="gene">BRAF</span> V600E genotypes, or PAM50 subtypes), or apply downstream inferential methods including partial correlation, mixed-effects regression, and FDR control directly to package outputs. Maintaining a strict boundary between MIT-licensed software distributions and the Broad Institute DepMap non-commercial reference matrix keeps licensing compliance unambiguous across both academic scripts and institutional deployments.
  </p>

  <!-- Acknowledgements & Disclosures -->
  <h2>Acknowledgements and Disclosures</h2>
  <div class="disclosure-block">
    <p><strong>Acknowledgements:</strong> The author thanks Dr. Lin Deng and Dr. Zixi Wang (Institute of Molecular Physiology, Shenzhen Bay Laboratory) for constructive discussions and feedback.</p>
    <p><strong>Ethics Statement:</strong> This study analyzes publicly available, de-identified cancer cell-line CRISPR dependency scores from DepMap Public 26Q3 and synthetic software test fixtures; no new human subjects, clinical records, or animal experiments were involved, and institutional ethics approval is not applicable.</p>
    <p><strong>Data Source and Terms:</strong> Built-in dependency data derive from DepMap Public 26Q3 Chronos Gene Effect scores (released October 5, 2026; <a href="https://forum.depmap.org/t/announcing-the-26q3-release/4692">https://forum.depmap.org/t/announcing-the-26q3-release/4692</a>; Data Portal: <a href="https://depmap.org/portal/data_page/">https://depmap.org/portal/data_page/</a>). The bundled matrix is governed by Broad Institute DepMap Terms and Conditions (non-commercial research use only; reposted in <code>data/DEPMAP_TERMS.md</code> and inside the hosted web application), separate from the MIT-licensed source code.</p>
    <p><strong>AI Assistance Disclosure:</strong> Large language model coding assistants (Claude, Codex, and Gemini/Antigravity agents) assisted with software refactoring, automated testing, and manuscript editing. All methodological specifications, quantitative analyses, biological interpretations, and final editorial decisions were verified and approved by the author.</p>
    <p><strong>Funding:</strong> This work received no specific external funding.</p>
    <p><strong>Competing Interests:</strong> The author declares no competing financial or non-financial interests.</p>
    <p>
      <strong>Data and Code Availability:</strong> Application source code (<code>v1.2.0</code>, commit <code>5c18f9b</code>) is archived at <a href="https://github.com/ChanghaoKan/crispr-score-analyzer/releases/tag/v1.2.0">https://github.com/ChanghaoKan/crispr-score-analyzer/releases/tag/v1.2.0</a> and on Zenodo (version 1.2.0 DOI: <code>10.5281/zenodo.24458535</code>; concept DOI: <code>10.5281/zenodo.19607602</code>); the hosted web instance runs at <a href="https://crispr-score-analyzer-szbl-denglab.streamlit.app/">https://crispr-score-analyzer-szbl-denglab.streamlit.app/</a>. Because Python wheels/sdists and MIT release archives exclude <code>data/depmap_chronos.npz</code> under DepMap terms, the modular <code>crispr_score_analyzer</code> library, CLI, manuscript sources, and statistical reproduction scripts (<code>preprint/generate_figures.py</code> and <code>preprint/reproduce_stats.py</code>) are pinned at repository tag <code>preprint-v4</code> (analytical kernels unchanged from <code>v1.2.0</code>). From a clean Python 3.12 environment, all figures and quantitative statistics in Tables&nbsp;1–2 and Figures&nbsp;1–3 can be reproduced via <code class="wrap-code">git clone https://github.com/ChanghaoKan/crispr-score-analyzer.git &amp;&amp; cd crispr-score-analyzer &amp;&amp; git checkout preprint-v4 &amp;&amp; python -m venv .venv &amp;&amp; .venv/bin/pip install -e ".[app,dev]" &amp;&amp; .venv/bin/python preprint/generate_figures.py &amp;&amp; .venv/bin/python preprint/reproduce_stats.py</code>, verifying <code>data/depmap_chronos.npz</code> against <code>data/depmap_chronos.provenance.json</code> (<code>matrix_sha256</code>: <code class="wrap-code">ddf274da2aa5f89cf630bc8150dcd3a2c722466bd18ff6fb6d18fe161f8a0c67</code>, or rebuilding from upstream 26Q3 files via <code>scripts/build_builtin_matrix.py</code>).
    </p>
  </div>

  <!-- References -->
  <h2>References</h2>
  <ul class="references">
    <li>Boyle, E. A., Pritchard, J. K., &amp; Greenleaf, W. J. (2018). High-resolution mapping of cancer cell networks using co-functional interactions. <em>Molecular Systems Biology</em>, 14(12), e8594. <a href="https://doi.org/10.15252/msb.20188594">https://doi.org/10.15252/msb.20188594</a></li>
    <li>Dempster, J. M., Boyle, I., Vazquez, F., Root, D. E., Boehm, J. S., Hahn, W. C., Tsherniak, A., &amp; McFarland, J. M. (2021). Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. <em>Genome Biology</em>, 22(1), 343. <a href="https://doi.org/10.1186/s13059-021-02540-7">https://doi.org/10.1186/s13059-021-02540-7</a></li>
    <li>DepMap, Broad. (2026). DepMap Public 26Q3 Chronos Gene Effect dataset (released October 5, 2026). Broad Institute. <a href="https://depmap.org/portal/data_page/">https://depmap.org/portal/data_page/</a> (Release notes: <a href="https://forum.depmap.org/t/announcing-the-26q3-release/4692">https://forum.depmap.org/t/announcing-the-26q3-release/4692</a>)</li>
    <li>Gheorghe, V., &amp; Hart, T. (2022). Optimal construction of a functional interaction network from pooled library CRISPR fitness screens. <em>BMC Bioinformatics</em>, 23(1), 510. <a href="https://doi.org/10.1186/s12859-022-05078-y">https://doi.org/10.1186/s12859-022-05078-y</a></li>
    <li>Kan, C. (2026). CRISPR Score Analyzer (1.2.0) [Software]. Zenodo. <a href="https://doi.org/10.5281/zenodo.24458535">https://doi.org/10.5281/zenodo.24458535</a> (Concept DOI: <a href="https://doi.org/10.5281/zenodo.19607602">https://doi.org/10.5281/zenodo.19607602</a>)</li>
    <li>Meyers, R. M., Bryan, J. G., McFarland, J. M., Weir, B. A., Sizemore, A. E., Xu, H., <em>et al.</em> (2017). Computational correction of copy number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells. <em>Nature Genetics</em>, 49(12), 1779–1784. <a href="https://doi.org/10.1038/ng.3984">https://doi.org/10.1038/ng.3984</a></li>
    <li>Pacini, C., Dempster, J. M., Boyle, I., Gon&ccedil;alves, E., Najgebauer, H., Karakoc, E., <em>et al.</em> (2021). Integrated cross-study datasets of genetic dependencies in cancer. <em>Nature Communications</em>, 12(1), 1661. <a href="https://doi.org/10.1038/s41467-021-21898-7">https://doi.org/10.1038/s41467-021-21898-7</a></li>
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
      <span>Kan (2026) &middot; <em>CRISPR Score Analyzer</em></span>
      <span>Preprint</span>
    </div>
    """
    footer_tpl = """
    <div style="width:100%; font-family:'Liberation Sans',Arial,sans-serif; font-size:7pt; color:#86868b; padding:0 17.5mm; display:flex; justify-content:space-between; align-items:center;">
      <span>DepMap Public 26Q3 &middot; Software &amp; Methods Preprint</span>
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
