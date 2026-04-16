<div align="center">

# 🧬 CRISPR Score Analyzer

### An interactive platform for gene essentiality analysis using DepMap CRISPR screens

[![Streamlit App](https://img.shields.io/badge/Streamlit-Live%20Demo-FF4B4B?logo=streamlit&logoColor=white)](https://crispr-score-analyzer.streamlit.app)
[![HuggingFace Data](https://img.shields.io/badge/🤗-Dataset-yellow)](https://huggingface.co/datasets/ChanghaoKan/crispr-depmap)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9+-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![DOI](https://zenodo.org/badge/1125701784.svg)](https://doi.org/10.5281/zenodo.19607602)

**[🚀 Launch App](https://crispr-score-analyzer.streamlit.app)** · 
**[📊 Dataset](https://huggingface.co/datasets/ChanghaoKan/crispr-depmap)** · 
**[🐛 Report Issue](https://github.com/ChanghaoKan/crispr-score-analyzer/issues)**

</div>

---

## 📖 Overview

**CRISPR Score Analyzer** is a lightweight, interactive web application for exploring gene essentiality from the [DepMap Cancer Dependency Map](https://depmap.org) (Chronos-corrected CRISPR screens across 1,000+ cancer cell lines). It lets researchers quickly locate genes of interest within a genome-wide dependency ranking, visualize lineage-specific dependencies, and export publication-quality figures — all without writing a line of code.

The tool is built for wet-lab biologists, early-stage computational researchers, and anyone preparing figures for manuscripts involving CRISPR dependency data.

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
| *(add screenshot here)* | *(add screenshot here)* |

</div>

## 🚀 Quick Start

### Option 1: Use the live app (recommended)

Simply visit **[crispr-score-analyzer.streamlit.app](https://crispr-score-analyzer.streamlit.app)** — no installation, no login.

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

The app will open at `http://localhost:8501`.

## 📊 Data Source

This tool uses the **DepMap Public 25Q3 Chronos CRISPR** dependency scores, subsetted for efficient cloud delivery. The full dataset is hosted on Hugging Face:

🔗 **[huggingface.co/datasets/ChanghaoKan/crispr-depmap](https://huggingface.co/datasets/ChanghaoKan/crispr-depmap)**

| Property | Value |
|---|---|
| Source | [DepMap Portal](https://depmap.org/portal/) (Broad Institute) |
| Release | 25Q3 |
| Score type | Chronos (corrected CRISPR dependency) |
| Cell lines | ~1,100 |
| Genes | ~18,000 |
| Score interpretation | Lower = more essential; ≈0 = neutral; >0 = enriched |

## 📋 Expected CSV Format

If uploading your own data, use this schema:

| `cell_line` | `lineage` | `MYC (4609)` | `PTEN (5728)` | ... |
|---|---|---|---|---|
| ACH-000001 | Liver | -1.95 | 0.12 | ... |
| ACH-000002 | Lung | -2.13 | 0.08 | ... |

- **First columns**: metadata (cell line ID, lineage, etc.)
- **Gene columns**: numeric CRISPR scores, headers in `GENE_SYMBOL (ENTREZ_ID)` format
- At least one column containing the word `lineage` is needed for boxplot analysis

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
  author    = {Kan, Changhao and Deng Lab},
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

### DepMap data

> Tsherniak, A., Vazquez, F., Montgomery, P. G., *et al.* (2017). Defining a Cancer Dependency Map. *Cell* 170, 564–576. https://doi.org/10.1016/j.cell.2017.06.010

## 🤝 Contributing

Contributions, feature requests, and bug reports are welcome. Please open an [Issue](https://github.com/ChanghaoKan/crispr-score-analyzer/issues) or submit a Pull Request.

For substantial changes, please open an issue first to discuss the scope.

## 📄 License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

DepMap data is released under the terms of its own license; please consult the [DepMap Portal](https://depmap.org) for details.

## 🙏 Acknowledgements

- **[DepMap Consortium](https://depmap.org)** — for the open release of Chronos-corrected CRISPR dependency data
- **[Anthropic Claude](https://www.anthropic.com/claude)** — AI-assisted development
- **Deng Lab @ Shenzhen Bay Laboratory (SZBL)** — institutional support

## 📬 Contact

**Changhao Kan**  
Research Assistant, Deng Lab  
Shenzhen Bay Laboratory (SZBL)  
📧 [kch_ynu@163.com] *(add your preferred contact)*  
🌐 [GitHub](https://github.com/ChanghaoKan)

---

<div align="center">

**⭐ If you find this tool useful, please consider starring the repository.**

*Made with Streamlit, Plotly, and ☕ in Shenzhen.*

</div>
