#!/usr/bin/env python3
"""Build preprint/manuscript.pdf in clean plain-white academic manuscript style.

Typesets standard Times-metric serif (Liberation Serif), 1-inch margins,
1.5-spaced eLife-structured prose, unstructured abstract, vector Figure 1
architecture schematic, unshaded display equations, booktabs-style three-line
tables, frameless figures, and plain centered page numbers via headless
Chromium / Playwright.

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
    mime = "image/svg+xml" if path.suffix.lower() == ".svg" else "image/png"
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{mime};base64,{data}"


def build_html() -> str:
    fig1_path = FIGURES_DIR / "Figure1.svg"
    if not fig1_path.is_file():
        fig1_path = FIGURES_DIR / "Figure1.png"
    fig1_uri = img_to_data_uri(fig1_path)
    fig2_uri = img_to_data_uri(FIGURES_DIR / "Figure2.png")
    fig3_uri = img_to_data_uri(FIGURES_DIR / "Figure3.png")

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<title>CRISPR Score Analyzer: cohort-scoped essentiality and co-dependency analysis of DepMap CRISPR-Cas9 screens</title>
<style>
  @page {{
    size: A4 portrait;
    margin: 25.4mm;
  }}

  * {{
    box-sizing: border-box;
  }}

  html, body {{
    margin: 0;
    padding: 0;
    background: #ffffff;
    color: #000000;
    font-family: "Liberation Serif", "Times New Roman", Times, Georgia, serif;
    font-size: 11pt;
    line-height: 1.5;
    -webkit-font-smoothing: antialiased;
    text-rendering: optimizeLegibility;
  }}

  /* Plain Academic Title Block */
  .title-block {{
    text-align: center;
    margin: 0 0 18pt 0;
  }}

  h1.title {{
    font-size: 15pt;
    line-height: 1.3;
    font-weight: 700;
    color: #000000;
    margin: 0 0 11pt 0;
  }}

  .author-name {{
    font-size: 11.5pt;
    font-weight: 400;
    color: #000000;
    margin: 0 0 4pt 0;
    text-align: center;
  }}

  .affiliations {{
    font-size: 10pt;
    line-height: 1.45;
    color: #000000;
    margin: 0;
    text-align: center;
  }}

  /* Unstructured Plain Abstract */
  .abstract-section {{
    margin: 0 0 14pt 0;
  }}

  .abstract-heading {{
    font-size: 12pt;
    font-weight: 700;
    color: #000000;
    margin: 0 0 6pt 0;
    text-align: left;
    break-after: avoid;
    page-break-after: avoid;
  }}

  .keywords {{
    font-size: 10.5pt;
    line-height: 1.45;
    color: #000000;
    margin: 8pt 0 0 0;
    text-align: left;
    hyphens: none;
  }}

  /* Section Headings (pure black, no colored accents or rules) */
  h2 {{
    font-size: 12pt;
    line-height: 1.3;
    font-weight: 700;
    color: #000000;
    margin: 16pt 0 6pt 0;
    break-after: avoid;
    page-break-after: avoid;
  }}

  h3 {{
    font-size: 11pt;
    line-height: 1.35;
    font-weight: 700;
    color: #000000;
    margin: 12pt 0 5pt 0;
    break-after: avoid;
    page-break-after: avoid;
  }}

  p {{
    margin: 0 0 7.5pt 0;
    color: #000000;
    text-align: justify;
    hyphens: auto;
    orphans: 3;
    widows: 3;
  }}

  /* Inline & Display Mathematical Typography (plain white, zero shading) */
  .math {{
    font-family: "Liberation Serif", "Times New Roman", Times, Georgia, serif;
    font-size: 11pt;
    color: #000000;
    white-space: nowrap;
  }}

  .math i {{
    font-style: italic;
  }}

  .math sub, .math sup {{
    font-size: 70%;
    line-height: 0;
  }}

  .eq-row {{
    display: flex;
    align-items: center;
    justify-content: space-between;
    margin: 9pt 0 10pt 0;
    padding: 2pt 0;
    background: transparent;
    border: none;
    break-inside: avoid;
    page-break-inside: avoid;
  }}

  .eq-body {{
    flex: 1;
    text-align: center;
    font-family: "Liberation Serif", "Times New Roman", Times, Georgia, serif;
    font-size: 11pt;
    line-height: 1.35;
    color: #000000;
  }}

  .eq-num {{
    font-family: "Liberation Serif", "Times New Roman", Times, serif;
    font-size: 10.5pt;
    color: #000000;
    padding-left: 10pt;
    white-space: nowrap;
  }}

  .eq-gap {{
    display: inline-block;
    width: 22pt;
  }}

  .frac {{
    display: inline-flex;
    flex-direction: column;
    align-items: center;
    vertical-align: middle;
    margin: 0 2.5pt;
    line-height: 1.15;
  }}

  .frac .num {{
    padding: 1.5pt 3.5pt 2.2pt 3.5pt;
    border-bottom: 0.65pt solid #000000;
    font-size: 10.2pt;
  }}

  .frac .den {{
    padding: 2.2pt 3.5pt 0.8pt 3.5pt;
    font-size: 10.2pt;
  }}

  .frac.frac-tall .num {{
    padding: 2pt 4.5pt 6.5pt 4.5pt;
  }}

  .frac.frac-tall .den {{
    padding: 4.5pt 4.5pt 4.5pt 4.5pt;
  }}

  .sum-op {{
    display: inline-flex;
    flex-direction: column;
    align-items: center;
    vertical-align: -0.34em;
    margin: 0 2.5pt;
    line-height: 0.95;
  }}

  .sum-op .sym {{
    font-size: 14.5pt;
    font-style: normal;
  }}

  .sum-op .lim {{
    font-size: 7.4pt;
    margin-top: 1.6pt;
  }}

  .overline {{
    display: inline-block;
    border-top: 0.6pt solid #000000;
    padding-top: 0.5pt;
    line-height: 0.92;
  }}

  /* Inline Formatting (pure plain text, no pill backgrounds, black links) */
  .gene {{
    font-style: italic;
  }}

  .nowrap {{
    white-space: nowrap;
  }}

  code {{
    font-family: "Liberation Mono", "Courier New", Courier, monospace;
    font-size: 9.2pt;
    background: transparent;
    color: #000000;
    padding: 0;
    border-radius: 0;
    white-space: nowrap;
  }}

  code.code-brk {{
    white-space: normal;
  }}

  a {{
    color: #000000;
    text-decoration: none;
  }}

  strong {{
    color: #000000;
  }}

  sup, sub {{
    font-size: 70%;
    line-height: 0;
  }}

  /* Classic Three-Line Academic Tables (booktabs style, pure white) */
  .table-block {{
    margin: 12pt 0 13pt 0;
    break-inside: avoid;
    page-break-inside: avoid;
  }}

  .table-caption {{
    font-size: 9.8pt;
    line-height: 1.4;
    color: #000000;
    margin: 0 0 5pt 0;
    text-align: justify;
    hyphens: auto;
  }}

  table.three-line {{
    width: 100%;
    border-collapse: collapse;
    font-size: 8.8pt;
    line-height: 1.35;
    border-top: 1.2pt solid #000000;
    border-bottom: 1.2pt solid #000000;
    background: #ffffff;
  }}

  table.three-line thead th {{
    background: #ffffff;
    color: #000000;
    font-weight: 700;
    text-align: left;
    padding: 4.5pt 4.5pt;
    border-bottom: 0.65pt solid #000000;
    vertical-align: bottom;
  }}

  table.three-line thead th.num-col,
  table.three-line tbody td.num-col {{
    text-align: right;
    font-variant-numeric: tabular-nums;
    white-space: nowrap;
  }}

  table.three-line tbody td {{
    background: #ffffff;
    padding: 3.6pt 4.5pt;
    border: none;
    vertical-align: top;
    color: #000000;
  }}

  table.three-line td.cap-col {{
    width: 24%;
  }}

  table.three-line td.csa-col {{
    width: 32%;
  }}

  table.three-line code {{
    font-size: 7.8pt;
  }}

  table.three-line.compact thead th,
  table.three-line.compact tbody td {{
    padding: 3.6pt 3.2pt;
  }}

  /* Frameless Embedded Figures on White Paper */
  .figure-block {{
    margin: 12pt 0 13pt 0;
    padding: 0;
    border: none;
    border-radius: 0;
    background: #ffffff;
    break-inside: avoid;
    page-break-inside: avoid;
  }}

  .figure-block img {{
    display: block;
    width: 100%;
    height: auto;
    margin: 0 auto 6pt auto;
  }}

  .figure-block.fig2 img,
  .figure-block.fig3 img {{
    width: 96%;
  }}

  .figure-caption {{
    font-size: 9.6pt;
    line-height: 1.4;
    color: #000000;
    margin: 0;
    text-align: justify;
    hyphens: auto;
  }}

  .figure-caption code {{
    font-size: 8.4pt;
  }}

  /* Additional Information */
  .declarations-block p {{
    font-size: 10.2pt;
    line-height: 1.45;
    margin-bottom: 5.5pt;
    text-align: justify;
    hyphens: auto;
  }}

  .declarations-block p.avail-p {{
    text-align: left;
    hyphens: none;
  }}

  .declarations-block code {{
    font-size: 8.6pt;
  }}

  /* Hanging-Indent Author-Year References */
  ul.references {{
    list-style: none;
    padding-left: 0;
    margin: 4pt 0 0 0;
  }}

  ul.references li {{
    font-size: 9.8pt;
    line-height: 1.42;
    color: #000000;
    padding-left: 18pt;
    text-indent: -18pt;
    margin-bottom: 4.5pt;
    text-align: left;
    hyphens: none;
    overflow-wrap: break-word;
  }}
</style>
</head>
<body>

  <!-- Title & Author Block -->
  <div class="title-block">
    <h1 class="title">CRISPR Score Analyzer: cohort-scoped essentiality and <span class="nowrap">co-dependency</span> analysis of DepMap <span class="nowrap">CRISPR-Cas9</span> screens</h1>
    <p class="author-name">Changhao Kan<sup>1,*</sup></p>
    <p class="affiliations">
      <sup>1</sup>Institute of Molecular Physiology, Shenzhen Bay Laboratory, Shenzhen, China<br />
      <sup>*</sup>Correspondence: <a href="mailto:kanchanghao@szbl.ac.cn">kanchanghao@szbl.ac.cn</a> &nbsp;|&nbsp;
      ORCID: <a href="https://orcid.org/0009-0000-4014-5923">0009-0000-4014-5923</a>
    </p>
  </div>

  <!-- Unstructured Abstract -->
  <div class="abstract-section">
    <h2 class="abstract-heading">Abstract</h2>
    <p>
      Genome-scale CRISPR-Cas9 loss-of-function screens in the Cancer Dependency Map, or DepMap, provide a standard reference for prioritizing cancer vulnerabilities and mapping functional gene relationships. Beyond single-gene portal queries, however, exploratory analyses routinely confront three methodological problems: pan-cancer averages dilute lineage- and subtype-restricted dependencies, multi-library integration produces uneven gene coverage across cell lines, and pan-cancer co-essentiality profiles conflate pathway co-functionality with baseline differences between tumor lineages. Here we present CRISPR Score Analyzer, an open-source Python package, command-line interface, and browser platform built on a single analytical core over DepMap Public 26Q3 Chronos scores spanning 1,361 cancer cell lines, 30 Oncotree lineages, and 18,434 genes, alongside validated custom screen matrices. The analytical engine provides cohort-scoped essentiality ranking, joint mean and median lineage-selectivity decomposition, exact pairwise-complete <code>float64</code> Pearson and average-rank Spearman genome-wide scans in 2,048-column blocks without SciPy, offline HGNC symbol harmonization across 40,402 alias mappings, vector graphics export, and provenance archiving. Auditing DepMap Public 26Q3 shows that 1,352 genes, or 7.33% of the matrix, lack full cell-line coverage with <code>n_valid</code> falling as low as 218, motivating explicit pairwise overlap reporting alongside every correlation. Chunked kernels agree with direct reference calculations to machine precision, within <span class="nowrap">7.77&nbsp;&times;&nbsp;10<sup>&minus;16</sup></span> for Pearson <em>r</em> and <span class="nowrap">5.55&nbsp;&times;&nbsp;10<sup>&minus;17</sup></span> for Spearman <em>&rho;</em>, and complete a pan-cancer Pearson scan in 0.74&nbsp;s or a single-lineage Skin scan in 0.10&nbsp;s. Across ten canonical cancer gene pairs, cohort scoping separates pan-lineage epistatic axes such as <span class="gene">MDM2</span>–<span class="gene">TP53</span> from context-restricted modules whose pan-cancer rank correlations collapse toward zero, including <span class="gene">BRAF</span>–<span class="gene">MAP2K1</span> and <span class="gene">ESR1</span>–<span class="gene">GATA3</span>. Source code is available under the MIT license at <a href="https://github.com/ChanghaoKan/crispr-score-analyzer">https://github.com/<wbr />ChanghaoKan/<wbr />crispr-score-analyzer</a> with a hosted application at <a href="https://crispr-score-analyzer-szbl-denglab.streamlit.app/">https://crispr-score-analyzer-szbl-denglab.<wbr />streamlit.app/</a>.
    </p>
    <p class="keywords">
      <strong>Keywords:</strong> DepMap; CRISPR-Cas9 screen; Chronos gene effect; gene essentiality; co-dependency; functional genomics
    </p>
  </div>

  <!-- Section 1: Introduction -->
  <h2>1. Introduction</h2>
  <p>
    Pooled genome-scale CRISPR-Cas9 loss-of-function screens across human cancer cell lines provide a primary experimental reference for identifying lineage-specific vulnerabilities and mapping functional gene relationships (Meyers <em>et al.</em>, 2017; Pacini <em>et al.</em>, 2021). Within the Cancer Dependency Map, or DepMap (Tsherniak <em>et al.</em>, 2017), copy-number- and growth-kinetics-corrected Chronos gene effect scores (Dempster <em>et al.</em>, 2021) allow candidate genes nominated by transcriptomics, proteomics, or focused genetic screens to be evaluated rapidly against pan-cancer essentiality distributions.
  </p>
  <p>
    Moving beyond single-gene portal lookups, however, exposes three recurring analytical limitations. First, pan-cancer averages frequently obscure lineage- and subtype-restricted dependencies: an oncogene or lineage transcription factor strongly essential in one tumor type can appear near-neutral across more than a thousand models, while intra-lineage heterogeneity can separate lineage means from medians. Second, integrating screens conducted with different sgRNA libraries introduces uneven gene coverage across cell-line subsets (Pacini <em>et al.</em>, 2021), so genome-wide correlation rankings computed without tracking pairwise sample overlap can place sparsely profiled genes alongside full-cohort associations. Third, pan-cancer co-essentiality profiles conflate direct pathway co-functionality with baseline mean differences between tumor lineages. Covariance-whitening and dimensionality-reduction approaches address correlated sample structure when building whole-genome co-essentiality networks (Boyle <em>et al.</em>, 2018; Wainberg <em>et al.</em>, 2021; Gheorghe &amp; Hart, 2022), yet rely on offline matrix transformations rather than direct inspection within a user-defined lineage or disease subset. Existing web resources such as the DepMap Portal (<a href="https://depmap.org/portal/">https://depmap.org/portal/</a>) and shinyDepMap (Shimada <em>et al.</em>, 2021) serve single-gene queries and global overviews well, whereas multi-gene batch positioning, single-lineage co-dependency scanning, paired cell-line selection, and reproducible scripted execution remain fragmented across separate interfaces and custom scripts.
  </p>
  <p>
    Here we describe CRISPR Score Analyzer, an open-source Python library (<code class="code-brk">crispr_<wbr />score_<wbr />analyzer</code>), command-line interface (<code>crispr-score-analyzer</code>), and interactive browser platform built on a shared analytical core over DepMap Public 26Q3 Chronos scores and validated custom matrices. Rather than introducing a latent-factor normalization model, the software provides transparent, cohort-scoped descriptive estimators: exact pairwise-complete Pearson and Spearman scanning without imputation, joint mean and median lineage-selectivity summaries linked directly to individual cell-line identities, offline HGNC symbol harmonization, and self-contained provenance archives. Using DepMap Public 26Q3, we audit gene-coverage missingness, verify numerical accuracy and runtime scaling of the chunked correlation kernels, and benchmark pan-cancer against cohort-scoped co-dependency across ten canonical cancer gene pairs spanning pan-lineage epistatic axes, lineage-restricted signaling cascades, and subtype-dependent transcriptional modules.
  </p>

  <!-- Section 2: Results (eLife order) -->
  <h2>2. Results</h2>

  <h3>2.1 Architecture and analytical scope of CRISPR Score Analyzer</h3>
  <p>
    CRISPR Score Analyzer separates analytical computation from user-interface presentation through a modular, Streamlit-free Python library (<code class="code-brk">crispr_<wbr />score_<wbr />analyzer</code>; archived at release <code>v1.2.0</code>, Kan, 2026) that drives both a scriptable CLI and a bilingual browser application (Figure 1). All cohort-scoped essentiality rankings, lineage-selectivity decompositions, and pairwise-complete co-dependency scans run in <code>float64</code> precision on top of NumPy and pandas without requiring SciPy. Offline HGNC alias resolution harmonizes batch gene queries across 40,402 unambiguous symbol mappings, while plotting and export routines produce white-background vector graphics and self-contained provenance archives. Relative to existing interactive DepMap portals, the platform unifies multi-gene cohort positioning, single-step lineage-scoped co-dependency scanning, joint mean–median skew diagnostics, custom screen ingestion, and scripted Python execution on identical estimators (Table 1).
  </p>

  <!-- Figure 1 (Vector SVG architecture schematic) -->
  <div class="figure-block fig1">
    <img src="{fig1_uri}" alt="Figure 1" />
    <p class="figure-caption">
      <strong>Figure 1. Software architecture and analytical workflow of CRISPR Score Analyzer.</strong>
      <strong>A.</strong> Data ingestion and offline nomenclature layer supporting the bundled DepMap Public 26Q3 Chronos matrix, validated user-supplied UTF-8 CSV screens, and 40,402 unambiguous HGNC alias mappings.
      <strong>B.</strong> Shared <code>crispr_score_analyzer</code> Python analytical core executing <code>float64</code> cohort-scoped essentiality ranking, joint mean and median lineage-selectivity decomposition, and 2,048-column pairwise-complete Pearson and Spearman co-dependency scans without SciPy.
      <strong>C.</strong> Dual execution interfaces—a bilingual Streamlit browser platform and a scriptable Python API and CLI—producing white-background vector figures and self-contained SHA-256 provenance archives.
    </p>
  </div>

  <!-- Table 1 -->
  <div class="table-block">
    <p class="table-caption">
      <strong>Table 1.</strong> Comparison of analytical and reproducibility capabilities across CRISPR Score Analyzer and existing DepMap interfaces, verified November 2026.
    </p>
    <table class="three-line">
      <thead>
        <tr>
          <th>Analytical capability</th>
          <th>DepMap Portal (26Q3)</th>
          <th>shinyDepMap (eLife 2021)</th>
          <th>CRISPR Score Analyzer (<code>v1.2.0</code> / <code>preprint-v6</code>)</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td class="cap-col">Multi-gene panel placement on cohort curve</td>
          <td>Single-gene default</td>
          <td>Efficacy vs. selectivity scatter</td>
          <td class="csa-col">Batch ranking and two-layer gene-set annotation</td>
        </tr>
        <tr>
          <td class="cap-col">Cohort-scoped genome-wide co-dependency scan</td>
          <td>Pan-cancer default; multi-step subsetting</td>
          <td>Pan-cancer connectivity clusters</td>
          <td class="csa-col">Single-step scoped scan with exact pairwise Pearson <em>r</em> and Spearman <em>&rho;</em></td>
        </tr>
        <tr>
          <td class="cap-col">Pairwise multi-gene matrix and lineage scatter</td>
          <td>Multi-panel explorer</td>
          <td>Cluster network view</td>
          <td class="csa-col">2–25 gene heatmap linked to lineage-colored scatter</td>
        </tr>
        <tr>
          <td class="cap-col">Lineage selectivity and intra-lineage skew</td>
          <td>Per-gene lineage and subtype plots</td>
          <td>Not provided</td>
          <td class="csa-col">Joint &Delta;<em>&mu;</em>, &Delta;med, and fraction &lt; &minus;0.5 lollipop and heatmap</td>
        </tr>
        <tr>
          <td class="cap-col">Sensitive and resistant cell-line selection</td>
          <td>Table filtering on gene page</td>
          <td>Top sensitive lines</td>
          <td class="csa-col">Paired most- and least-dependent tables with <code>ModelID</code> export</td>
        </tr>
        <tr>
          <td class="cap-col">Custom matrix ingestion and schema audit</td>
          <td>Single-feature vector</td>
          <td>Not supported</td>
          <td class="csa-col">Full UTF-8 CSV matrix with schema guard and diagnostics</td>
        </tr>
        <tr>
          <td class="cap-col">Offline batch HGNC harmonization</td>
          <td>Online single-gene search</td>
          <td>Not provided</td>
          <td class="csa-col">40,402 unambiguous offline mappings for batch and single queries</td>
        </tr>
        <tr>
          <td class="cap-col">Scriptable Python API, CLI, and provenance</td>
          <td>Web portal and bulk downloads</td>
          <td>R Shiny web app only</td>
          <td class="csa-col">Installable Python library, CLI, vector export, and SHA-256 archive</td>
        </tr>
      </tbody>
    </table>
  </div>

  <h3>2.2 Coverage heterogeneity, numerical fidelity, and runtime scaling in DepMap Public 26Q3</h3>
  <p>
    Auditing the DepMap Public 26Q3 Chronos matrix across 1,361 cancer cell lines and 18,434 genes reveals substantial gene-level missingness arising from multi-library integration. Pan-cancer gene means span &minus;5.250 to +0.398, with 1,712 genes falling below &minus;0.5. Although 17,082 genes, or 92.67%, are profiled across all 1,361 cell lines, 1,352 genes, or 7.33%, have incomplete coverage ranging down to <span class="nowrap"><code>n_valid</code>&nbsp;=&nbsp;218</span>. Within this incomplete subset, 648 genes are assayed in fewer than 1,000 cell lines and 519 in fewer than 500. In genome-wide correlation scans over such heterogeneous matrices, failing to preserve exact pairwise overlap or omitting <code>n_shared</code> from ranked outputs risks conflating sparse sub-cohort associations with full-panel correlations.
  </p>
  <p>
    We verified numerical equivalence of the chunked <code>float64</code> correlation engine against direct reference evaluations using <code>pandas.Series.corr</code> for Pearson <em>r</em> and pairwise midrank <code>numpy.corrcoef</code> for Spearman <em>&rho;</em> across top, bottom, and partially observed partners of <span class="gene">MDM2</span>—193 genes for Pearson and 100 for Spearman. Maximum absolute deviations were <span class="nowrap">7.77&nbsp;&times;&nbsp;10<sup>&minus;16</sup></span> for Pearson <em>r</em> and <span class="nowrap">5.55&nbsp;&times;&nbsp;10<sup>&minus;17</sup></span> for Spearman <em>&rho;</em>, confirming machine-precision agreement across both complete and incomplete gene columns. On a 2-vCPU Linux environment, median wall-clock runtimes across three replicates were 1.01&nbsp;s to load the binary matrix and 0.41&nbsp;s to rank all 18,434 genes. A genome-wide pan-cancer scan across all 1,361 models and 18,433 partner genes finished in 0.74&nbsp;s for Pearson <em>r</em> and 6.19&nbsp;s for Spearman <em>&rho;</em>. Restricting the cohort to a single lineage, Skin (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;75</span>), reduced scan runtimes to 0.10&nbsp;s for Pearson and 0.44&nbsp;s for Spearman, enabling interactive exploration across both metrics.
  </p>

  <h3>2.3 Lineage-selectivity decomposition resolves subtype-driven skew</h3>
  <p>
    Positioning a five-gene benchmark panel—<span class="gene">KRAS</span>, <span class="gene">SOX10</span>, <span class="gene">BRAF</span>, <span class="gene">ESR1</span>, and <span class="gene">TP53</span>—on the pan-cancer Chronos ranking curve illustrates how global averages compress context-restricted vulnerabilities (Figure 2A). Relative to reference anchors <span class="gene">MYC</span> at mean &minus;1.971 (rank 187, 1.01%) and <span class="gene">PTEN</span> at mean +0.377 (rank 18,428, 99.97%), <span class="gene">KRAS</span> falls below the &minus;0.5 visual reference across all lineages at mean &minus;0.729 (rank 1,263, 6.85%), whereas <span class="gene">TP53</span> sits at the positive extreme at mean +0.381 (rank 18,430, 99.98%), reflecting positive growth-rate effects upon tumor-suppressor knockout in <span class="gene">TP53</span>-wild-type lines (Dempster <em>et al.</em>, 2021). By contrast, two canonical lineage-dependent targets appear modest or near-neutral in pan-cancer averages: <span class="gene">SOX10</span> averages &minus;0.347 (rank 2,225, 12.07%) and <span class="gene">ESR1</span> averages &minus;0.065 (rank 7,347, 39.86%). Two-layer gene-set annotation similarly positions mitotic kinases <span class="gene">PLK1</span> at mean &minus;2.691 (rank 38, 0.21%) and <span class="gene">AURKA</span> at mean &minus;1.236 (rank 626, 3.40%) against seven cell-cycle CDK and cyclin genes ranging from <span class="gene">CDK1</span> at &minus;2.588 to <span class="gene">CCNE1</span> at &minus;0.194 (Figure 2B).
  </p>

  <!-- Figure 2 -->
  <div class="figure-block fig2">
    <img src="{fig2_uri}" alt="Figure 2" />
    <p class="figure-caption">
      <strong>Figure 2. Multi-gene essentiality positioning and lineage-selectivity decomposition in DepMap Public 26Q3.</strong>
      <strong>A.</strong> Pan-cancer mean Chronos gene effect curve across 1,361 cell lines and 18,434 genes positioning five benchmark genes (<span class="gene">KRAS</span>, <span class="gene">SOX10</span>, <span class="gene">BRAF</span>, <span class="gene">ESR1</span>, <span class="gene">TP53</span>) relative to <span class="gene">MYC</span> and <span class="gene">PTEN</span> anchors and the &minus;0.5 reference guide.
      <strong>B.</strong> Two-layer gene-set annotation overlaying <span class="gene">PLK1</span> and <span class="gene">AURKA</span> onto seven cell-cycle CDK and cyclin genes (<span class="gene">CDK1</span>, <span class="gene">CDK2</span>, <span class="gene">CDK4</span>, <span class="gene">CDK6</span>, <span class="gene">CCNB1</span>, <span class="gene">CCND1</span>, <span class="gene">CCNE1</span>).
      <strong>C, D.</strong> Lineage-preference lollipop plots showing top and bottom four lineages with <span class="nowrap"><em>n</em>&nbsp;&ge;&nbsp;5</span> on gene-specific scales, displaying mean Chronos shifts &Delta;<em>&mu;</em> relative to each gene’s all-lineage mean for <span class="gene">SOX10</span> (<strong>C</strong>) and <span class="gene">ESR1</span> (<strong>D</strong>).
    </p>
  </div>

  <p>
    Evaluating mean shift &Delta;<em>&mu;</em> jointly with median shift &Delta;med and sensitive fraction <span class="math"><i>f</i><sub><i>g</i>,&ell;</sub></span> separates uniform lineage essentiality from subtype-restricted bimodality (Figure 2C, D). For <span class="gene">SOX10</span>, Skin (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;75</span>) shows concordant mean and median shifts—mean &minus;1.487, <span class="nowrap">&Delta;<em>&mu;</em>&nbsp;=&nbsp;&minus;1.140</span>, <span class="nowrap">&Delta;med&nbsp;=&nbsp;&minus;1.368</span>, and 78.7% of lines below &minus;0.5—consistent with broad neural-crest dependency in cutaneous melanoma. In Eye (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;15</span>), however, the mean shift of <span class="nowrap">&Delta;<em>&mu;</em>&nbsp;=&nbsp;&minus;0.535</span> (mean &minus;0.882, 46.7% below &minus;0.5) is more than threefold stronger than the median shift of <span class="nowrap">&Delta;med&nbsp;=&nbsp;&minus;0.174</span>. This divergence flags a bimodal distribution driven by strongly dependent uveal melanoma models such as <code>Omm2.5</code> at &minus;2.752, <code>UPMM3</code> at &minus;2.113, and <code>WM3772F</code> at &minus;1.678 mixed with near-neutral retinoblastoma lines. A comparable subtype-masking effect arises for <span class="gene">ESR1</span> in Breast (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;67</span>), where Oncotree lineage grouping pools ER-positive luminal models with ER-negative subtypes. The resulting Breast shift—mean &minus;0.285, <span class="nowrap">&Delta;<em>&mu;</em>&nbsp;=&nbsp;&minus;0.220</span>, <span class="nowrap">&Delta;med&nbsp;=&nbsp;&minus;0.061</span>, and 17.9% of lines below &minus;0.5—is only slightly larger than a nine-line Ampulla of Vater group (<span class="nowrap">&Delta;<em>&mu;</em>&nbsp;=&nbsp;&minus;0.174</span>, <span class="nowrap">&Delta;med&nbsp;=&nbsp;&minus;0.059</span>). Inspecting the ranked cell-line table resolves this aggregation immediately by surfacing strongly dependent ER-positive luminal models led by <code>HCC1428</code> at &minus;1.912, <code>T-47D</code> at &minus;1.845, and <code>MDA-MB-415</code> at &minus;1.636 alongside resistant controls.
  </p>

  <h3>2.4 Cohort scope reshapes co-dependency inference across canonical cancer axes</h3>
  <p>
    To quantify how cohort scope alters co-dependency inference, we benchmarked ten established cancer gene pairs across three estimation regimes in DepMap Public 26Q3: pan-cancer Pearson <em>r</em> and Spearman <em>&rho;</em> across all 1,361 models, pooled within-lineage Pearson correlation <span class="math"><i>r</i><sub>within</sub></span> after subtracting lineage means across all 30 lineages, and single-lineage Pearson <em>r</em> and Spearman <em>&rho;</em> within each pair’s primary biological context (Table 2).
  </p>

  <!-- Table 2 -->
  <div class="table-block">
    <p class="table-caption">
      <strong>Table 2.</strong> Pan-cancer, lineage-demeaned, and target-lineage co-dependency across ten canonical cancer gene pairs in DepMap Public 26Q3 (all pan-cancer pairs <span class="nowrap"><em>n</em>&nbsp;=&nbsp;1,361</span>).
    </p>
    <table class="three-line compact">
      <thead>
        <tr>
          <th>Gene pair</th>
          <th>Target lineage (<em>n</em>)</th>
          <th class="num-col">Pan <em>r</em></th>
          <th class="num-col">Pan <em>&rho;</em></th>
          <th class="num-col"><em>r</em><sub>within</sub></th>
          <th class="num-col">Lineage <em>r</em></th>
          <th class="num-col">Lineage <em>&rho;</em></th>
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
          <td class="num-col">&minus;0.802</td>
          <td class="num-col">&minus;0.500</td>
          <td>Pan-lineage epistatic</td>
        </tr>
        <tr>
          <td><span class="gene">MDM2</span>–<span class="gene">PPM1D</span></td>
          <td>Lung (126)</td>
          <td class="num-col">+0.664</td>
          <td class="num-col">+0.658</td>
          <td class="num-col">+0.635</td>
          <td class="num-col">+0.642</td>
          <td class="num-col">+0.584</td>
          <td>Pan-lineage epistatic</td>
        </tr>
        <tr>
          <td><span class="gene">BRAF</span>–<span class="gene">MAPK1</span></td>
          <td>Skin (75)</td>
          <td class="num-col">+0.644</td>
          <td class="num-col">+0.451</td>
          <td class="num-col">+0.538</td>
          <td class="num-col">+0.798</td>
          <td class="num-col">+0.787</td>
          <td>Lineage-enriched pathway</td>
        </tr>
        <tr>
          <td><span class="gene">BRAF</span>–<span class="gene">MAP2K1</span></td>
          <td>Skin (75)</td>
          <td class="num-col">+0.391</td>
          <td class="num-col">+0.097</td>
          <td class="num-col">+0.303</td>
          <td class="num-col">+0.552</td>
          <td class="num-col">+0.605</td>
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
          <td class="num-col">+0.503</td>
          <td class="num-col">+0.346</td>
          <td>Subtype-restricted TF axis</td>
        </tr>
        <tr>
          <td><span class="gene">ESR1</span>–<span class="gene">GATA3</span></td>
          <td>Breast (67)</td>
          <td class="num-col">+0.253</td>
          <td class="num-col">+0.095</td>
          <td class="num-col">+0.228</td>
          <td class="num-col">+0.596</td>
          <td class="num-col">+0.372</td>
          <td>Subtype-restricted TF axis</td>
        </tr>
        <tr>
          <td><span class="gene">CDK4</span>–<span class="gene">CCND1</span></td>
          <td>Breast (67)</td>
          <td class="num-col">+0.416</td>
          <td class="num-col">+0.380</td>
          <td class="num-col">+0.425</td>
          <td class="num-col">+0.704</td>
          <td class="num-col">+0.635</td>
          <td>Context-amplified cell cycle</td>
        </tr>
        <tr>
          <td><span class="gene">KRAS</span>–<span class="gene">RAF1</span></td>
          <td>Bowel (85)</td>
          <td class="num-col">+0.338</td>
          <td class="num-col">+0.255</td>
          <td class="num-col">+0.319</td>
          <td class="num-col">+0.530</td>
          <td class="num-col">+0.461</td>
          <td>Lineage-amplified effector</td>
        </tr>
      </tbody>
    </table>
  </div>

  <p>
    The benchmark separates three distinct statistical regimes. First, pan-lineage epistatic pairs such as <span class="gene">MDM2</span>–<span class="gene">TP53</span> and <span class="gene">MDM2</span>–<span class="gene">PPM1D</span> show close agreement across pan-cancer <em>r</em>, pan-cancer <em>&rho;</em>, and lineage-demeaned <span class="math"><i>r</i><sub>within</sub></span>: for <span class="gene">MDM2</span>–<span class="gene">TP53</span>, the three estimators yield &minus;0.702, &minus;0.704, and &minus;0.648. In a genome-wide pan-cancer scan of <span class="gene">MDM2</span> across 18,433 genes (Figure 3A), top positive partners include <span class="gene">PPM1D</span> (+0.664), <span class="gene">MDM4</span> (+0.614), <span class="gene">USP7</span> (+0.522), <span class="gene">DDX31</span> (+0.466), and <span class="gene">UBE2D3</span> (+0.398), whereas top anti-correlated partners recover the p53 tumor-suppressor axis—<span class="gene">TP53</span> (&minus;0.702), <span class="gene">TP53BP1</span> (&minus;0.589), <span class="gene">USP28</span> (&minus;0.531), <span class="gene">CDKN1A</span> (&minus;0.530), <span class="gene">CHEK2</span> (&minus;0.491), and <span class="gene">ATM</span> (&minus;0.461)—all profiled across all 1,361 cell lines. Annotating <code>n_shared</code> directly on the figure flags rank-6 positive hit <span class="gene">TERF1</span> (+0.395) as partially covered across only 228 shared lines, immediately separating it from full-cohort associations. Multi-gene matrix inspection of the six-gene regulatory axis resolves its two-block inverse structure (Figure 3B). Across all 23 lineages containing at least 15 models, <span class="gene">MDM2</span> and <span class="gene">TP53</span> remained negatively correlated, with a median Pearson <em>r</em> of &minus;0.680 spanning &minus;0.892 in Liver to &minus;0.246 in Cervix (Figure 3C).
  </p>
  <p>
    Second, context-restricted signaling and transcriptional modules suffer sharp attenuation in pan-cancer screens, visible directly in divergence between pan-cancer Pearson <em>r</em> and Spearman <em>&rho;</em>. When a pathway operates in a small minority of cell lines, extreme scores in the active lineage generate moderate Pearson covariance while more than 90% of non-expressing lines contribute rank noise near zero, collapsing pan-cancer <em>&rho;</em> toward zero. Thus pan-cancer <em>r</em> and <em>&rho;</em> diverge from +0.391 to +0.097 for <span class="gene">BRAF</span>–<span class="gene">MAP2K1</span>, from +0.253 to +0.095 for <span class="gene">ESR1</span>–<span class="gene">GATA3</span>, and from +0.324 to +0.111 for <span class="gene">ESR1</span>–<span class="gene">FOXA1</span>. Restricting the cohort to the relevant lineage restores parametric and rank concordance simultaneously, elevating <span class="gene">BRAF</span>–<span class="gene">MAP2K1</span> in Skin to <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.552</span> and <span class="nowrap"><em>&rho;</em>&nbsp;=&nbsp;+0.605</span>, <span class="gene">ESR1</span>–<span class="gene">GATA3</span> in Breast to <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.596</span> and <span class="nowrap"><em>&rho;</em>&nbsp;=&nbsp;+0.372</span>, <span class="gene">CDK4</span>–<span class="gene">CCND1</span> in Breast to <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.704</span> and <span class="nowrap"><em>&rho;</em>&nbsp;=&nbsp;+0.635</span>, and <span class="gene">KRAS</span>–<span class="gene">RAF1</span> in Bowel to <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.530</span> and <span class="nowrap"><em>&rho;</em>&nbsp;=&nbsp;+0.461</span>.
  </p>

  <!-- Figure 3 -->
  <div class="figure-block fig3">
    <img src="{fig3_uri}" alt="Figure 3" />
    <p class="figure-caption">
      <strong>Figure 3. Pan-cancer and lineage-scoped genome-wide co-dependency benchmarks in DepMap Public 26Q3.</strong>
      <strong>A.</strong> Pan-cancer Pearson <em>r</em> scan of <span class="gene">MDM2</span> across 18,433 partner genes (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;1,361</span>), displaying the top seven positive (<span class="gene">PPM1D</span>, <span class="gene">MDM4</span>, <span class="gene">USP7</span>, <span class="gene">DDX31</span>, <span class="gene">UBE2D3</span>, partially profiled <span class="gene">TERF1</span> at <span class="nowrap"><em>n</em>&nbsp;=&nbsp;228</span>, <span class="gene">PPM1G</span>) and top seven negative (<span class="gene">TP53</span>, <span class="gene">TP53BP1</span>, <span class="gene">USP28</span>, <span class="gene">CDKN1A</span>, <span class="gene">CHEK2</span>, <span class="gene">ATM</span>, <span class="gene">UBE2K</span>) partners.
      <strong>B.</strong> Pairwise Pearson correlation heatmap for the six-gene <span class="gene">TP53</span>–<span class="gene">MDM2</span> regulatory panel.
      <strong>C.</strong> Lineage-colored cell-line scatter plot of <span class="gene">MDM2</span> versus <span class="gene">TP53</span> (<span class="nowrap"><em>r</em>&nbsp;=&nbsp;&minus;0.702</span>, Spearman <span class="nowrap"><em>&rho;</em>&nbsp;=&nbsp;&minus;0.704</span>, <span class="nowrap"><em>n</em>&nbsp;=&nbsp;1,361</span>) with OLS trendline.
      <strong>D.</strong> Single-lineage genome-wide <span class="gene">BRAF</span> scan restricted to Skin models (<span class="nowrap"><em>n</em>&nbsp;=&nbsp;75</span>, 17,785 eligible partners), showing elevated MAPK cascade correlations alongside melanoma lineage factors and exploratory hits near <span class="nowrap">|<em>r</em>|&nbsp;&approx;&nbsp;0.50</span>.
    </p>
  </div>

  <p>
    Third, between-lineage baseline co-enrichment inflates pan-cancer Pearson correlations for lineage-restricted transcription factors. After lineage demeaning, <span class="gene">SOX10</span>–<span class="gene">MITF</span> drops from pan-cancer <span class="nowrap"><em>r</em>&nbsp;=&nbsp;+0.492</span> to <span class="nowrap"><span class="math"><i>r</i><sub>within</sub></span>&nbsp;=&nbsp;+0.286</span> alongside pan-cancer <span class="nowrap"><em>&rho;</em>&nbsp;=&nbsp;+0.153</span>, and <span class="gene">BRAF</span>–<span class="gene">SOX10</span> drops from +0.544 to <span class="nowrap"><span class="math"><i>r</i><sub>within</sub></span>&nbsp;=&nbsp;+0.294</span>. Single-lineage scoping, however, does not automatically remove intra-lineage subtype confounding. In the Skin-restricted <span class="gene">BRAF</span> scan across 75 models and 17,785 eligible partners with <span class="nowrap"><code>n_shared</code>&nbsp;&ge;&nbsp;10</span>—excluding 648 sparser genes (Figure 3D)—MAPK effectors strengthen markedly relative to pan-cancer values: <span class="gene">MAPK1</span> rises from +0.644 to +0.798, <span class="gene">MAP2K1</span> from +0.391 to +0.552, and <span class="gene">MAP2K2</span> from +0.271 to +0.550. Yet melanoma lineage factors <span class="gene">SOX10</span> (+0.633) and <span class="gene">SOX5</span> (+0.562) remain among the top hits alongside exploratory candidates near <span class="nowrap">|<em>r</em>|&nbsp;&approx;&nbsp;0.50</span>, including <span class="gene">STYXL1</span> (&minus;0.534), <span class="gene">UBQLN1</span> (+0.526), and <span class="gene">AKIRIN1</span> (+0.496), all complete across 75 models. Regressing both <span class="gene">BRAF</span> and each partner on <span class="gene">SOX10</span> via ordinary least squares across all 75 Skin lines retains strong MAPK cascade partial correlations—+0.665 for <span class="gene">MAPK1</span>, +0.467 for <span class="gene">MAP2K1</span>, and +0.398 for <span class="gene">MAP2K2</span>—indicating that <span class="gene">SOX10</span> co-ranking stems from Oncotree Skin pooling of <span class="gene">SOX10</span>-high cutaneous melanoma with <span class="gene">SOX10</span>-low non-melanoma or dedifferentiated lines. Post-hoc two-sided Pearson <em>t</em>-tests across all 17,785 Skin pairs further separate <span class="gene">MAPK1</span> (<span class="nowrap"><em>p</em>&nbsp;=&nbsp;9.72&nbsp;&times;&nbsp;10<sup>&minus;18</sup></span>, Bonferroni <span class="nowrap">1.73&nbsp;&times;&nbsp;10<sup>&minus;13</sup></span>) from the three moderate candidates near <span class="nowrap">|<em>r</em>|&nbsp;&approx;&nbsp;0.50</span>, which span nominal <em>p</em> from <span class="nowrap">8.28&nbsp;&times;&nbsp;10<sup>&minus;7</sup></span> to <span class="nowrap">5.95&nbsp;&times;&nbsp;10<sup>&minus;6</sup></span>, Bonferroni-adjusted <em>p</em> from 0.015 to 0.106, and Benjamini–Hochberg FDR from 0.0025 to 0.0132, illustrating why moderate single-lineage correlations warrant independent validation rather than immediate mechanistic interpretation.
  </p>

  <!-- Section 3: Discussion -->
  <h2>3. Discussion</h2>
  <p>
    Co-essentiality analyses of large-scale CRISPR screens trade off global statistical normalization against local biological interpretability. Covariance-whitening methods such as generalized least squares (Wainberg <em>et al.</em>, 2021), olfactory-receptor principal-component removal (Boyle <em>et al.</em>, 2018), and PCA-whitened network pipelines (Gheorghe &amp; Hart, 2022) effectively suppress shared technical variance and dominant mitochondrial or proliferation axes when constructing genome-wide functional networks. Global linear whitening, however, transforms the native Chronos scale and can obscure whether an association holds across tissues or arises from a small subset of models. Cohort-scoped descriptive analysis serves a complementary purpose: operating directly on unwhitened Chronos scores within an explicit lineage or disease subset, reporting both Pearson <em>r</em> and Spearman <em>&rho;</em> alongside exact <code>n_shared</code>, and linking aggregate correlations directly to individual cell-line scatter plots lets investigators inspect biological context before proceeding to network modeling or wet-lab follow-up.
  </p>
  <p>
    The DepMap Public 26Q3 benchmarks also clarify practical boundaries on interpreting descriptive outputs. Cohort ranks, lineage mean and median shifts, and pairwise correlations summarize empirical screen behavior rather than formal tests of differential dependency, selectivity, or synthetic lethality, and the &minus;0.5 line serves solely as a visual screening reference rather than a common-essential classifier. As illustrated by <span class="gene">ESR1</span> in Breast, <span class="gene">SOX10</span> in Eye, and <span class="gene">BRAF</span>–<span class="gene">SOX10</span> in Skin, Oncotree lineage labels pool distinct molecular subtypes, while scanning roughly 18,000 genes across 75 cell lines introduces substantial multiple-testing variance around <span class="nowrap">|<em>r</em>|&nbsp;&approx;&nbsp;0.50</span>. Furthermore, the bundled 26Q3 matrix intentionally stores only Oncotree lineage and cell-line identifiers without somatic mutation, expression, or receptor covariates, and user-uploaded custom tables are assumed to contain one row per unique cell line without silent row deduplication.
  </p>
  <p>
    Decoupling the <code class="code-brk">crispr_<wbr />score_<wbr />analyzer</code> library and CLI from the Streamlit interface addresses these boundaries in extended workflows. Investigators who start interactively to position a candidate gene panel or select sensitive and resistant cell lines can transition directly to scripted Python pipelines backed by identical ranking and correlation routines, merge external multi-omic covariates such as <span class="gene">TP53</span> mutation status, <span class="gene">BRAF</span> V600E genotypes, or PAM50 subtypes, and apply downstream inferential models including partial correlation, mixed-effects regression, and false discovery rate control. Keeping MIT-licensed software distributions strictly separate from the Broad Institute DepMap non-commercial reference matrix also preserves clear licensing compliance across academic pipelines and institutional deployments.
  </p>

  <!-- Section 4: Materials and Methods (eLife order) -->
  <h2>4. Materials and Methods</h2>

  <h3>4.1 Reference dataset curation and nomenclature harmonization</h3>
  <p>
    To support reproducible local and hosted execution without runtime network fetches, the repository bundles a compressed <code>float32</code> binary matrix compiled from DepMap Public 26Q3 Chronos Gene Effect scores and matching cell-line metadata, released October 5, 2026 (DepMap, 2026). All 1,361 cancer cell lines and 18,434 genes were retained verbatim after validating unique <code>ModelID</code> rows and DepMap-style <code>SYMBOL (Entrez ID)</code> column headers, then left-joined on <code>ModelID</code> against deduplicated model metadata to attach <code>OncotreeLineage</code>—matched across all 1,361 models and 30 lineages—and <code>CellLineName</code>, defaulting to <code>ModelID</code> where blank. Quantization to <code>float32</code> introduces a maximum absolute rounding error of <span class="nowrap">1.19&nbsp;&times;&nbsp;10<sup>&minus;7</sup></span> across the first 100 genes while reducing storage more than fivefold; upstream and compiled checksums and build metadata are archived in the repository provenance manifest.
  </p>
  <p>
    Software and dataset licensing are maintained separately. Application source code and the public-domain HGNC alias map distributed in Python packages are released under the MIT license, whereas rehosting the compiled DepMap matrix in the repository and hosted web instance is governed by Broad Institute DepMap Terms and Conditions for non-commercial research use only. Python wheel builds and MIT release archives exclude the bundled DepMap matrix.
  </p>
  <p>
    For user-supplied UTF-8 CSV CRISPR score matrices, automatic gene-column recognition requires <code>SYMBOL (Entrez ID)</code> headers so numeric metadata columns are never silently parsed as genes; plain-symbol tables require explicit column selection in the sidebar or CLI. Ambiguous duplicate gene symbols and repeated metadata fields are excluded automatically and logged in a downloadable diagnostics table. Disease-subtype filtering activates automatically whenever <code class="code-brk">Oncotree<wbr />PrimaryDisease</code> or <code>primary_disease</code> is present in custom uploads. Gene queries accept arbitrary whitespace or punctuation delimiters, text or delimited lists, and pasted DepMap headers, and are harmonized offline via 40,402 unambiguous HGNC previous-symbol and alias mappings such as <span class="gene">HER2</span>&nbsp;&rarr;&nbsp;<span class="gene">ERBB2</span> (Seal <em>et al.</em>, 2023). Resolved aliases, duplicates, and unmatched symbols are reported explicitly with every run.
  </p>

  <h3>4.2 Cohort-scoped essentiality and lineage-selectivity estimators</h3>
  <p>
    Let <span class="math"><b>X</b> &isin; (&#8477; &cup; {{NA}})<sup><i>M</i>&times;<i>G</i></sup></span> denote the cell-line-by-gene Chronos matrix, and let <span class="math"><i>C</i> &sube; {{1, &hellip;, <i>M</i>}}</span> denote the active cohort defined by the selected cancer lineages and optional disease subtypes. For each gene <span class="math"><i>g</i></span> with finite observation set <span class="math"><i>C</i><sub><i>g</i></sub> = {{<i>i</i> &isin; <i>C</i> : <i>x</i><sub><i>ig</i></sub> &isin; &#8477;}}</span> and sample count <span class="math"><i>n</i><sub><i>g,C</i></sub> = |<i>C</i><sub><i>g</i></sub>| &gt; 0</span>, cohort-scoped essentiality is summarized by its <code>float64</code> arithmetic mean and ascending empirical percentile rank across all <span class="math"><i>G</i><sub><i>C</i></sub></span> valid genes:
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
    Non-finite entries are excluded from sums and tracked separately in missingness diagnostics. Because Chronos scores are scaled so that non-essential genes center near 0 and core essential genes center near &minus;1 (Dempster <em>et al.</em>, 2021), a dashed visual guide at &minus;0.5 and two configurable reference anchors—defaulting to <span class="gene">MYC</span> and <span class="gene">PTEN</span>—orient multi-gene ranking curves and two-layer gene-set overlays without implying a binary classifier.
  </p>
  <p>
    To quantify lineage selectivity while diagnosing intra-lineage skew, the lineage-preference estimator compares each lineage <span class="math">&ell;</span> meeting <span class="math"><i>n</i><sub><i>g</i>,&ell;</sub> &ge; 5</span> against the full-dataset distribution of gene <span class="math"><i>g</i></span>, independent of active sidebar filters, using three complementary statistics:
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
    where <span class="math">med</span> denotes the sample median and <span class="math"><b>I</b>(&middot;)</span> is the indicator function. Comparing <span class="math">&Delta;<i>&mu;</i><sub><i>g</i>,&ell;</sub></span> against <span class="math">&Delta;med<sub><i>g</i>,&ell;</sub></span> distinguishes broadly shared lineage dependencies from bimodal or subtype-restricted vulnerabilities where a sensitive subpopulation pulls the mean well below the median. Paired cell-line tables expose the most and least dependent models with display names and <code>ModelID</code> identifiers within <span class="math"><i>C</i></span> for direct experimental selection.
  </p>

  <h3>4.3 Pairwise-complete vectorized co-dependency estimation</h3>
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
          <span style="font-size:13.5pt; vertical-align:-0.08em;">[</span>&thinsp;<span class="sum-op"><span class="sym">&sum;</span><span class="lim"><i>i</i> &isin; <i>S</i><sub><i>q,g</i></sub></span></span>(<i>x</i><sub><i>iq</i></sub> &minus; <span class="overline"><i>x</i></span><sub><i>q</i></sub><sup>(<i>S</i>)</sup>)<sup>2</sup>
          &thinsp;&middot;&thinsp;
          <span class="sum-op"><span class="sym">&sum;</span><span class="lim"><i>i</i> &isin; <i>S</i><sub><i>q,g</i></sub></span></span>(<i>x</i><sub><i>ig</i></sub> &minus; <span class="overline"><i>x</i></span><sub><i>g</i></sub><sup>(<i>S</i>)</sup>)<sup>2</sup>&thinsp;<span style="font-size:13.5pt; vertical-align:-0.08em;">]</span><sup>1/2</sup>
        </span>
      </span>
    </div>
    <div class="eq-num">(3)</div>
  </div>

  <p>
    where <span class="math"><span class="overline"><i>x</i></span><sub><i>q</i></sub><sup>(<i>S</i>)</sup></span> and <span class="math"><span class="overline"><i>x</i></span><sub><i>g</i></sub><sup>(<i>S</i>)</sup></span> denote means computed over <span class="math"><i>S</i><sub><i>q,g</i></sub></span>. Spearman rank correlation <span class="math"><i>&rho;</i><sub><i>q,g</i></sub>(<i>C</i>)</span> applies the identical formula after replacing <span class="math"><i>x</i><sub>&middot;<i>q</i></sub></span> and <span class="math"><i>x</i><sub>&middot;<i>g</i></sub></span> with average-tie midranks on <span class="math"><i>S</i><sub><i>q,g</i></sub></span>.
  </p>
  <p>
    To maintain interactive latency with bounded memory and no SciPy dependency, genome-wide scanning processes candidate genes in 2,048-column <code>float64</code> blocks after restricting rows to finite observations of <span class="math"><i>q</i></span>. Complete columns within each block are evaluated via a centered BLAS matrix–vector inner product, whereas incomplete columns use masked pairwise means and variances for Pearson <em>r</em> and per-column exact pairwise midranks for Spearman <em>&rho;</em>. Time complexity is <span class="math"><i>O</i>(|<i>C</i>|<i>G</i>)</span> for Pearson scans with <span class="math"><i>O</i>(2048&nbsp;&middot;&nbsp;|<i>C</i>|)</span> temporary memory. For multi-gene panels of 2 to 25 genes, the same pairwise estimator populates a symmetric correlation matrix linked to lineage-colored cell-line scatter plots reporting both Pearson <em>r</em> and Spearman <em>&rho;</em>.
  </p>

  <h3>4.4 Software implementation, post-hoc benchmarking, and verification</h3>
  <p>
    All analytical routines are packaged in the PEP 621 library <code class="code-brk">crispr_<wbr />score_<wbr />analyzer</code> requiring only NumPy, pandas, Plotly, Kaleido, and Pillow, exposing four modules for essentiality ranking, co-dependency scanning, plotting, and graphics export. The shared core powers the <code>crispr-score-analyzer</code> CLI and a bilingual English/Chinese Streamlit web interface. All interactive plots and white-background PDF, SVG, and 300-DPI PNG exports apply colorblind-safe Nature Publishing Group categorical hues and a 7-stop ColorBrewer <code>RdBu</code> diverging palette, converting WebGL scatter curves to vector paths on export. Input edits remain isolated until executed, and cohort changes invalidate downstream tables and exports to avoid stale state. Each session can export a self-contained ZIP archive containing ranked tables, lineage-preference summaries, schema and alias diagnostics, and run metadata with data SHA-256 digests and software versions. Software behavior is verified offline on deterministic synthetic fixtures by 98 automated unit, API, CLI, UI, and export tests across Python 3.11 and 3.12 via GitHub Actions.
  </p>
  <p>
    All figures and quantitative metrics reported in Section&nbsp;2 are reproduced offline without SciPy by <code class="code-brk">preprint/<wbr />generate_figures.py</code> and <code class="code-brk">preprint/<wbr />reproduce_stats.py</code>. Post-hoc two-sided Pearson <em>p</em>-values in the Skin-restricted <span class="gene">BRAF</span> scan are evaluated from the exact Student-<em>t</em> distribution with <span class="math">df = <i>n</i><sub><i>q,g</i></sub> &minus; 2</span> via continued-fraction evaluation of the regularized incomplete beta function, alongside Bonferroni and Benjamini–Hochberg FDR adjustments and ordinary least-squares partial correlations on complete observations.
  </p>

  <!-- Additional Information -->
  <h2>Additional Information</h2>
  <div class="declarations-block">
    <p><strong>Acknowledgements:</strong> The author thanks Dr. Lin Deng and Dr. Zixi Wang at the Institute of Molecular Physiology, Shenzhen Bay Laboratory, for constructive discussions and feedback.</p>
    <p><strong>Funding:</strong> This work received no specific external funding.</p>
    <p><strong>Competing interests:</strong> The author declares no competing financial or non-financial interests.</p>
    <p><strong>Ethics:</strong> This study analyzes publicly available, de-identified cancer cell-line CRISPR dependency scores from DepMap Public 26Q3 and synthetic software test fixtures. No new human subjects, clinical records, or animal experiments were involved, and institutional ethics approval is not applicable.</p>
    <p><strong>AI assistance disclosure:</strong> Large language model coding assistants—Claude, Codex, and Gemini/Antigravity agents—assisted with software refactoring, automated testing, and manuscript editing. All methodological specifications, quantitative analyses, biological interpretations, and final editorial decisions were verified and approved by the author.</p>
    <p class="avail-p">
      <strong>Data and code availability:</strong> Application source code is released under the MIT license at <a href="https://github.com/ChanghaoKan/crispr-score-analyzer">https://github.com/ChanghaoKan/crispr-score-analyzer</a>, with archived release <code>v1.2.0</code> deposited on Zenodo (Kan, 2026; version DOI: <a href="https://doi.org/10.5281/zenodo.24458535">10.5281/zenodo.24458535</a>; concept DOI: <a href="https://doi.org/10.5281/zenodo.19607602">10.5281/zenodo.19607602</a>) and a hosted web application at <a href="https://crispr-score-analyzer-szbl-denglab.streamlit.app/">https://crispr-score-analyzer-szbl-denglab.streamlit.app/</a>. Built-in dependency data derive from DepMap Public 26Q3 Chronos Gene Effect scores (DepMap, 2026; <a href="https://depmap.org/portal/data_page/">https://depmap.org/portal/data_page/</a>) and are governed separately by Broad Institute DepMap Terms and Conditions for non-commercial research use only. The modular Python library, CLI, manuscript sources, and reproduction scripts (<code>preprint/generate_figures.py</code> and <code>preprint/reproduce_stats.py</code>) are pinned at repository tag <code>preprint-v6</code>, with upstream and compiled file checksums and setup instructions documented in the repository provenance manifest and README.
    </p>
  </div>

  <!-- References -->
  <h2>References</h2>
  <ul class="references">
    <li>Boyle, E. A., Pritchard, J. K., &amp; Greenleaf, W. J. (2018). High-resolution mapping of cancer cell networks using co-functional interactions. <em>Molecular Systems Biology</em>, 14(12), e8594. <a href="https://doi.org/10.15252/msb.20188594">https://doi.org/10.15252/msb.20188594</a></li>
    <li>Dempster, J. M., Boyle, I., Vazquez, F., Root, D. E., Boehm, J. S., Hahn, W. C., <em>et al.</em> (2021). Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. <em>Genome Biology</em>, 22(1), 343. <a href="https://doi.org/10.1186/s13059-021-02540-7">https://doi.org/10.1186/s13059-021-02540-7</a></li>
    <li>DepMap. (2026). DepMap Public 26Q3 Chronos Gene Effect dataset (released October 5, 2026). Broad Institute of MIT and Harvard. <a href="https://depmap.org/portal/data_page/">https://depmap.org/portal/data_page/</a> (Release notes: <a href="https://forum.depmap.org/t/announcing-the-26q3-release/4692">https://forum.depmap.org/t/announcing-the-26q3-release/4692</a>)</li>
    <li>Gheorghe, V., &amp; Hart, T. (2022). Optimal construction of a functional interaction network from pooled library CRISPR fitness screens. <em>BMC Bioinformatics</em>, 23(1), 510. <a href="https://doi.org/10.1186/s12859-022-05078-y">https://doi.org/10.1186/s12859-022-05078-y</a></li>
    <li>Kan, C. (2026). CRISPR Score Analyzer (v1.2.0) [Software]. Zenodo. <a href="https://doi.org/10.5281/zenodo.24458535">https://doi.org/10.5281/zenodo.24458535</a> (Concept DOI: <a href="https://doi.org/10.5281/zenodo.19607602">https://doi.org/10.5281/zenodo.19607602</a>)</li>
    <li>Meyers, R. M., Bryan, J. G., McFarland, J. M., Weir, B. A., Sizemore, A. E., Xu, H., <em>et al.</em> (2017). Computational correction of copy number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells. <em>Nature Genetics</em>, 49(12), 1779–1784. <a href="https://doi.org/10.1038/ng.3984">https://doi.org/10.1038/ng.3984</a></li>
    <li>Pacini, C., Dempster, J. M., Boyle, I., Gon&ccedil;alves, E., Najgebauer, H., Karakoc, E., <em>et al.</em> (2021). Integrated cross-study datasets of genetic dependencies in cancer. <em>Nature Communications</em>, 12(1), 1661. <a href="https://doi.org/10.1038/s41467-021-21898-7">https://doi.org/10.1038/s41467-021-21898-7</a></li>
    <li>Seal, R. L., Braschi, B., Gray, K., Jones, T. E. M., Tweedie, S., Haim-Vilmovsky, L., <em>et al.</em> (2023). Genenames.org: the HGNC resources in 2023. <em>Nucleic Acids Research</em>, 51(D1), D1003–D1009. <a href="https://doi.org/10.1093/nar/gkac888">https://doi.org/10.1093/nar/gkac888</a></li>
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

    header_tpl = "<div></div>"
    footer_tpl = """
    <div style="width:100%; font-family:'Liberation Serif','Times New Roman',Times,serif; font-size:9.5pt; color:#000000; text-align:center; padding-bottom:6mm;">
      <span class="pageNumber"></span>
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
                    "top": "25.4mm",
                    "bottom": "25.4mm",
                    "left": "25.4mm",
                    "right": "25.4mm",
                },
            )
            browser.close()
    finally:
        html_path.unlink(missing_ok=True)

    size_kb = OUTPUT_PDF.stat().st_size / 1024.0
    print(f"Wrote {OUTPUT_PDF} ({size_kb:.1f} KB)")


if __name__ == "__main__":
    main()
