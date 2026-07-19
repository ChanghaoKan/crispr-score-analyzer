<div align="center">

# 🧬 CRISPR Score Analyzer

### An interactive platform for gene essentiality analysis using DepMap CRISPR screens

[![Streamlit App](https://img.shields.io/badge/Streamlit-Live%20Demo-FF4B4B?logo=streamlit&logoColor=white)](https://crispr-score-analyzer-szbl-denglab.streamlit.app/)
[![HuggingFace Data](https://img.shields.io/badge/🤗-Dataset-yellow)](https://huggingface.co/datasets/ChanghaoKan/crispr-depmap)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9+-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![DOI](https://zenodo.org/badge/1125701784.svg)](https://doi.org/10.5281/zenodo.19607602)

**[🚀 Launch App](https://crispr-score-analyzer-szbl-denglab.streamlit.app/)** ·
**[📊 Dataset](https://huggingface.co/datasets/ChanghaoKan/crispr-depmap)** · 
**[🐛 Report Issue](https://github.com/ChanghaoKan/crispr-score-analyzer/issues)**

</div>

---

## 📖 Overview

**CRISPR Score Analyzer** is a lightweight, interactive web application for exploring gene essentiality from the [DepMap Cancer Dependency Map](https://depmap.org) (Chronos-corrected CRISPR screens across 1,000+ cancer cell lines). It lets researchers quickly locate genes of interest within a genome-wide dependency ranking, visualize lineage-specific dependencies, and export publication-quality figures — all without writing a line of code.

An experimental **CRISPR gene dependency × drug sensitivity** module remains in the source code for continued validation, but it is currently disabled in the public interface.

The tool is built for wet-lab biologists, early-stage computational researchers, and anyone preparing figures for manuscripts involving CRISPR dependency or drug-response data.

## ✨ Features

- **🔍 Gene dependency ranking** — locate any gene(s) on a genome-wide essentiality curve with reference anchors (MYC, PTEN customizable)
- **📦 Lineage-specific boxplots** — compare CRISPR scores across cancer types for selected genes
- **🎯 Multi-layer annotation** — overlay gene sets and highlight candidates with custom colors
- **📄 Publication-ready export** — download figures as **PDF (vector)**, **PNG (300 DPI)**, or **SVG (editable)**
- **🌐 Bilingual interface** — English / 中文 one-click switch
- **🎨 Light / Dark theme** — theme-aware UI and charts
- **📤 Bring your own data** — upload custom CSV if you have local DepMap-format scores
- **⚡ Cloud-hosted** — no installation required; data streamed from Hugging Face CDN

## 🖼️ Screenshots

<div align="center">

| Gene Ranking | Lineage Boxplot |
|:---:|:---:|
|<img width="935" height="634" alt="image" src="https://github.com/user-attachments/assets/89a2b476-d0bf-430a-b6f3-a2c3f84540a6" />|<img width="1375" height="625" alt="image" src="https://github.com/user-attachments/assets/18ff825e-bb35-4597-b145-47f059641242" />

</div>

> The public application exposes three validated gene-essentiality analysis tabs; the screenshots above show the ranking and lineage views.

## 🚀 Quick Start

### Option 1: Use the live app (recommended)

Simply visit **[crispr-score-analyzer-szbl-denglab.streamlit.app](https://crispr-score-analyzer-szbl-denglab.streamlit.app/)** — no installation, no login.

### Option 2: Run locally

```bash
# Clone the repository
git clone https://github.com/ChanghaoKan/crispr-score-analyzer.git
cd crispr-score-analyzer

# Install dependencies (Python 3.9+)
pip install -r requirements.txt

# Launch
streamlit run app.py
```


## 📊 Data Sources

All datasets are hosted on Hugging Face: 🔗 **[huggingface.co/datasets/ChanghaoKan/crispr-depmap](https://huggingface.co/datasets/ChanghaoKan/crispr-depmap)**

**Gene essentiality (Tabs 1–3)** uses the **DepMap Public 25Q3 Chronos CRISPR** dependency scores, subsetted for efficient cloud delivery.

The disabled experimental **Gene × Drug correlation** module additionally uses **DepMap Public 26Q1 CRISPR**, the **GDSC2 AUC** drug-response matrix, and the DepMap compound annotation table. All files are aligned by DepMap `ModelID` (`ACH-xxxxxx`).

| File (HF) | Role | Index / Key |
|---|---|---|
| `CRISPR_(DepMap_Public_25Q3+Score,_Chronos)_subsetted.csv` | Essentiality browser (Tabs 1–3) | `ModelID` × `GENE (entrez)` |
| `CRISPRGeneEffect_26Q1.csv` | Gene dependency for correlation | `ModelID` × `GENE (entrez)` |
| `GDSC2_AUC_Matrix.csv` | Drug sensitivity (AUC) | `ModelID` × `CompoundID (DPC-xxxxxx)` |
| `PortalCompounds.csv` | Compound → name / target lookup | `CompoundID` |

| Property | Value |
|---|---|
| Source | [DepMap Portal](https://depmap.org/portal/) (Broad Institute) · GDSC2 (Sanger) |
| CRISPR release | 25Q3 (browser) · 26Q1 (correlation) |
| Score type | Chronos (corrected CRISPR dependency) |
| CRISPR interpretation | Lower = more essential; ≈0 = neutral; >0 = enriched |
| GDSC2 AUC interpretation | Lower = more sensitive; higher = more resistant |

> **Score-direction note (important).** Both CRISPR Gene Effect and GDSC2 AUC are "lower = more extreme." A **positive** Spearman ρ therefore means *more gene-dependent lines are also more drug-sensitive* — the shared-vulnerability direction. The app states this explicitly with every result so the sign is never ambiguous.

## 📋 Expected CSV Formats

### Essentiality browser (Tabs 1–3)

If uploading your own data, use this schema:

| `cell_line` | `lineage` | `MYC (4609)` | `PTEN (5728)` | ... |
|---|---|---|---|---|
| ACH-000001 | Liver | -1.95 | 0.12 | ... |
| ACH-000002 | Lung | -2.13 | 0.08 | ... |

- **First columns**: metadata (cell line ID, lineage, etc.)
- **Gene columns**: numeric CRISPR scores, headers in `GENE_SYMBOL (ENTREZ_ID)` format
- At least one column containing the word `lineage` is needed for boxplot analysis

### Experimental Gene × Drug correlation (currently hidden)

These three files are read from the Hugging Face dataset repo (filenames are configurable at the top of `app.py`). The **first column** of both matrices must be the DepMap `ModelID`.

**`CRISPRGeneEffect_26Q1.csv`** — `ModelID` × gene

| ModelID | `DUSP6 (1848)` | `A1BG (1)` | ... |
|---|---|---|---|
| ACH-000001 | -0.42 | -0.05 | ... |

**`GDSC2_AUC_Matrix.csv`** — `ModelID` × compound

| ModelID | `DPC-000001` | `DPC-000002` | ... |
|---|---|---|---|
| ACH-000001 | 0.95 | 0.83 | ... |

**`PortalCompounds.csv`** — compound annotation (used for the searchable drug picker)

| CompoundID | CompoundName | GeneSymbolOfTargets | TargetOrMechanism | ... |
|---|---|---|---|---|
| DPC-000001 | (+)-CAMPTOTHECIN | TOP1 | TOPOISOMERASE INHIBITOR | ... |

> The drug picker matches your query against compound **name, ID, and target** — so searching `PARP` returns all PARP inhibitors even if you don't know their names.

## 🧭 Usage Guide

### Tab 1 — Gene Ranking

Paste a gene list (one per line or comma-separated) to locate your genes within the genome-wide dependency ranking. Useful for:
- Quickly assessing whether a pathway's members are broadly essential
- Preparing Figure 1–style panels showing candidate gene position

### Tab 2 — Lineage Boxplot

Visualize how a gene's dependency varies across cancer types. Useful for:
- Identifying lineage-selective dependencies
- Supporting lineage-specific therapeutic hypotheses

### Tab 3 — Multi-layer Annotation

Plot a background gene set (e.g., a pathway) with specific highlights (e.g., candidate drivers). Useful for:
- Showing how a few hits compare against their pathway context
- Creating layered figures for reviews and presentations

### Experimental Gene × Drug Correlation (currently hidden)

Pick a gene (CRISPR dependency) and a compound (GDSC2 sensitivity); the app aligns cell lines by `ModelID`, computes a **Spearman ρ**, and plots a scatter with a fitted line. Useful for:
- Testing shared-vulnerability hypotheses (e.g. *DUSP6* dependency vs. PARP-inhibitor response)
- Screening whether a dependency tracks with sensitivity to a drug class

**How to read it:**
- Both axes are "lower = more extreme" (more dependent / more sensitive).
- **ρ > 0** → more gene-dependent lines tend to be *more* drug-sensitive (hypothesis-consistent).
- **ρ < 0** → the opposite.
- The result text states the direction and significance for you.

> ⚠️ **Caveat (retained with the experimental module).** The analysis reports association rather than mechanism. Lineage and mutation background — e.g. *BRCA1/2* / HR status for PARP inhibitors — require explicit control when interpreting results.

## 📥 Figure Export

Every plot supports three formats:

| Format | Use Case |
|---|---|
| **PDF** ✨ | Manuscript submission (infinite scalable vector) |
| **PNG** | Presentations, posters (300 DPI scaled) |
| **SVG** | Post-editing in Illustrator / Inkscape |

Exports always use a white background regardless of UI theme.

## 🛠️ Tech Stack

- **Frontend / App framework**: [Streamlit](https://streamlit.io)
- **Visualization**: [Plotly](https://plotly.com/python/)
- **Figure export**: [Kaleido](https://github.com/plotly/Kaleido) (v0.2.1, self-contained)
- **Data hosting**: [Hugging Face Datasets](https://huggingface.co/docs/datasets/)
- **Deployment**: [Streamlit Community Cloud](https://streamlit.io/cloud)

## 📚 Citation

If you use CRISPR Score Analyzer in your research, please cite both the tool and the underlying DepMap data.

### This tool

```bibtex
@software{crispr_score_analyzer_2026,
  author    = {Kan, Changhao},
  title     = {CRISPR Score Analyzer: An interactive platform for DepMap gene essentiality},
  year      = {2026},
  version   = {v1.0.0},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.19607602},
  url       = {https://doi.org/10.5281/zenodo.19607602}
}
```

**APA format:**
> Kan, C. (2026). *CRISPR Score Analyzer: An interactive platform for DepMap gene essentiality* (Version v1.0.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.19607602

### DepMap & GDSC data

> Tsherniak, A., Vazquez, F., Montgomery, P. G., *et al.* (2017). Defining a Cancer Dependency Map. *Cell* 170, 564–576. https://doi.org/10.1016/j.cell.2017.06.010

> Yang, W., Soares, J., Greninger, P., *et al.* (2013). Genomics of Drug Sensitivity in Cancer (GDSC). *Nucleic Acids Research* 41, D955–D961. https://doi.org/10.1093/nar/gks1111

## 🤝 Contributing

Contributions, feature requests, and bug reports are welcome. Please open an [Issue](https://github.com/ChanghaoKan/crispr-score-analyzer/issues) or submit a Pull Request.

For substantial changes, please open an issue first to discuss the scope.

## 📄 License

This project is licensed under the MIT License — see the [LICENSE](License) file for details.

DepMap data is released under the terms of its own license; please consult the [DepMap Portal](https://depmap.org) for details.

## 🙏 Acknowledgements

- **[DepMap Consortium](https://depmap.org)** — for the open release of Chronos-corrected CRISPR dependency data
- **[GDSC / Sanger Institute](https://www.cancerrxgene.org)** — for the GDSC2 drug-response (AUC) dataset
- **[Anthropic Claude](https://www.anthropic.com/claude)** — AI-assisted development
- **[Deng Lab @ Shenzhen Bay Laboratory (SZBL)](https://www.deng-lab.org/resources)** — institutional support

## 📬 Contact

**Changhao Kan**  
Research Assistant, Deng Lab  
Shenzhen Bay Laboratory (SZBL)  
📧 [kch_ynu@163.com] *(add your preferred contact)*  
🌐 [GitHub](https://github.com/ChanghaoKan)

---
