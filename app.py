"""
CRISPR Score Analyzer - Public Gene Essentiality Analysis Tool
A multilingual, theme-aware interactive platform for DepMap CRISPR data.

Author: Changhao Kan
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
import copy
import hashlib
from html import escape
import json
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from analysis_core import (analyze_dataset, parse_gene_text, parse_gene_upload,
                           resolve_gene_list, lineage_preference,
                           find_lineage_column, get_lineage_frame,
                           lineage_summary, make_analysis_bundle, read_score_csv,
                           read_builtin_matrix, file_sha256, builtin_provenance_path)
from figure_export import export_fingerprint, render_figure_bytes, MAX_PNG_PIXELS

# =============================================================================
# 页面配置
# =============================================================================
st.set_page_config(
    page_title="CRISPR Score Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="auto"
)

# =============================================================================
# 内置矩阵。打开页面时从仓库读取。
# =============================================================================
BUILTIN_MATRIX = Path(__file__).resolve().parent / "data" / "depmap_chronos.npz"
GENE_ALIAS_FILE = Path(__file__).resolve().parent / "data" / "gene_aliases.json"
DATA_VERSION = "DepMap Public 26Q3"
SCORE_TYPE = "Chronos Gene Effect"
ESSENTIALITY_THRESHOLD = -0.5

# =============================================================================
# Citation / DOI
# 运行中的软件是 1.1.0。Zenodo 上最新的存档仍是 2026-04-16 的 1.0.0。
# =============================================================================
TOOL_VERSION = "1.1.0"
TOOL_RELEASED = "2026-11-21"
TOOL_AUTHORS = "Kan, C."
TOOL_YEAR = "2026"
ZENODO_CONCEPT_DOI = "10.5281/zenodo.19607602"
ZENODO_VERSION_DOI = "10.5281/zenodo.19607603"
ZENODO_VERSION = "1.0.0"
ZENODO_RELEASED = "2026-04-16"
GITHUB_URL = "https://github.com/ChanghaoKan/crispr-score-analyzer"

# =============================================================================
# 国际化
# =============================================================================
TRANSLATIONS = {
    'en': {
        'app_title': 'CRISPR Score Analyzer',
        'app_subtitle': 'Gene dependency analysis and visualization across cancer cell lines',
        'hero_kicker': 'DEPMAP · CRISPR SCREENING',
        'sidebar_settings': 'Settings',
        'language': 'Language',
        'theme': 'Theme',
        'light': 'Light',
        'dark': 'Dark',
        'data_source': 'Data Source',
        'upload_custom': 'Upload custom data (optional)',
        'upload_csv': 'Upload CRISPR Score CSV',
        'reference_genes': 'Reference Genes',
        'essential': 'Reference A',
        'nonessential': 'Reference B',
        'display_settings': 'Display',
        'show_labels': 'Show gene labels',
        'point_size': 'Point size',
        'export_size': 'Export Size',
        'export_height': 'Image height (px)',
        'cell_lines': 'Cell lines',
        'gene_count': 'Genes',
        'essential_genes': 'Mean score < {threshold}',
        'score_range': 'Score range',
        'resources': 'Citation & data source',
        'custom_dataset': 'Custom uploaded dataset',
        'score_guide': 'Lower score = stronger gene dependency',
        'essential_status': 'Mean-score screen',
        'essential_yes': 'Mean below cutoff',
        'essential_no': 'Mean not below cutoff',
        'first_eight_only': 'Showing the first 8 genes to keep the figure readable.',
        'lineage_missing': 'No lineage column was found in this dataset.',
        'tab1': 'Gene ranking',
        'tab2': 'Cancer-type comparison',
        'tab3': 'Gene-set annotation',
        'gene_ranking_title': 'Gene dependency ranking',
        'gene_ranking_desc': 'Locate your genes of interest in genome-wide CRISPR screens.',
        'input_target_genes': 'Input target genes',
        'input_method': 'Input method',
        'input_direct': 'Direct input',
        'input_file': 'Upload gene list',
        'gene_list': 'Gene list',
        'gene_list_help': 'One gene per line, or comma-separated',
        'matched': 'Matched',
        'not_found': 'Not found',
        'export_title': 'Export high-quality figure',
        'gene_details': 'Gene details',
        'boxplot_title': 'CRISPR Score Distribution by Cancer Type',
        'multilayer_title': 'Gene-set annotation',
        'bg_gene_set': 'Background gene set',
        'hl_gene_set': 'Highlight genes',
        'bg_color': 'Background color',
        'hl_color': 'Highlight color',
        'download_pdf': '📄 PDF (vector)',
        'download_png': '🖼️ PNG (300 DPI)',
        'download_svg': '✏️ SVG (editable)',
        'download_csv': '📊 Download data (CSV)',
        'download_hint': '💡 Click a button to generate the file. Recommend **PDF** for publications. Use **SVG** for Illustrator/Inkscape editing.',
        'no_data_warn': '⚠️ Please configure a data source or upload data',
        'loading_upload': 'Loading uploaded data...',
        'loading_hf': 'Loading the reference matrix...',
        'loaded': '✅ Loaded',
        'acknowledgements': 'Acknowledgements',
        'data_from': 'Data Source',
        'dev_with': 'Development Assistance',
        'ai_dev': 'AI assistance with code, interface design, testing and documentation',
        'citation': '📚 How to Cite',
        'cite_this_tool': 'Cite this tool',
    },
    'zh': {
        'app_title': 'CRISPR Score Analyzer',
        'app_subtitle': '癌症细胞系中的基因依赖分析与可视化',
        'hero_kicker': 'DEPMAP · CRISPR 筛选',
        'sidebar_settings': '设置',
        'language': '语言',
        'theme': '主题',
        'light': '浅色',
        'dark': '深色',
        'data_source': '数据来源',
        'upload_custom': '上传自定义数据（可选）',
        'upload_csv': '上传 CRISPR Score CSV',
        'reference_genes': '参考基因',
        'essential': '参考基因 A',
        'nonessential': '参考基因 B',
        'display_settings': '显示设置',
        'show_labels': '显示基因名标签',
        'point_size': '点大小',
        'export_size': '导出尺寸',
        'export_height': '图片高度 (px)',
        'cell_lines': '细胞系',
        'gene_count': '基因数',
        'essential_genes': '平均分 < {threshold}',
        'score_range': '分数范围',
        'resources': '引用与数据来源',
        'custom_dataset': '自定义上传数据',
        'score_guide': '分数越低，基因依赖越强',
        'essential_status': '平均分筛选',
        'essential_yes': '平均分低于阈值',
        'essential_no': '平均分未低于阈值',
        'first_eight_only': '为保证图表清晰，仅展示前 8 个基因。',
        'lineage_missing': '该数据中未找到 lineage 列。',
        'tab1': '基因排名',
        'tab2': '癌种比较',
        'tab3': '基因集标注',
        'gene_ranking_title': '基因依赖排名',
        'gene_ranking_desc': '在全基因组 CRISPR 筛选数据中定位您关注的基因。',
        'input_target_genes': '输入目标基因',
        'input_method': '输入方式',
        'input_direct': '直接输入',
        'input_file': '上传基因列表',
        'gene_list': '基因列表',
        'gene_list_help': '每行一个基因名，或用逗号分隔',
        'matched': '匹配成功',
        'not_found': '未找到',
        'export_title': '导出高质量图片',
        'gene_details': '基因详细信息',
        'boxplot_title': '按癌症类型的 CRISPR Score 分布',
        'multilayer_title': '基因集标注',
        'bg_gene_set': '背景基因集',
        'hl_gene_set': '高亮基因',
        'bg_color': '背景颜色',
        'hl_color': '高亮颜色',
        'download_pdf': '📄 PDF (矢量)',
        'download_png': '🖼️ PNG (300 DPI)',
        'download_svg': '✏️ SVG (可编辑)',
        'download_csv': '📊 下载数据表 (CSV)',
        'download_hint': '💡 点击按钮生成文件。论文投稿推荐 **PDF**（矢量图）。需要编辑选 **SVG**（可在 Illustrator/Inkscape 中修改）。',
        'no_data_warn': '⚠️ 请配置数据源或上传数据文件',
        'loading_upload': '正在加载上传的数据...',
        'loading_hf': '正在读取参考矩阵...',
        'loaded': '✅ 已加载',
        'acknowledgements': '致谢',
        'data_from': '数据来源',
        'dev_with': '开发协助',
        'ai_dev': 'AI 辅助代码开发、界面设计、测试与文档整理',
        'citation': '📚 引用方式',
        'cite_this_tool': '引用本工具',
    }
}


def t(key: str) -> str:
    lang = st.session_state.get('lang', 'en')
    return TRANSLATIONS[lang].get(key, key)


def load_gene_aliases() -> dict:
    """Uppercase previous/alias symbol to one approved symbol. Missing file yields no aliases."""
    if not GENE_ALIAS_FILE.is_file():
        return {}
    payload = json.loads(GENE_ALIAS_FILE.read_text(encoding="utf-8"))
    return {str(key).upper(): str(value).upper() for key, value in payload.get("map", {}).items()}


GENE_ALIASES = load_gene_aliases()


# =============================================================================
# Session State
# =============================================================================
if st.session_state.get('lang') not in TRANSLATIONS:
    st.session_state.lang = 'en'
if st.session_state.get('theme') not in {'light', 'dark'}:
    st.session_state.theme = 'light'


# =============================================================================
# 主题
# =============================================================================
THEMES = {
    'light': {
        'bg': '#f5f7fa',
        'bg_secondary': '#edf2f7',
        'bg_card': '#ffffff',
        'text': '#172231',
        'text_muted': '#607080',
        'border': '#dce4ec',
        'accent': '#245f96',
        'accent_hover': '#1f568d',
        'success': '#1f8a67',
        'danger': '#b54552',
        'shadow': '0 1px 3px rgba(31, 54, 78, 0.04)',
        'shadow_hover': '0 12px 30px rgba(31, 54, 78, 0.10)',
        'plot_bg': '#ffffff',
        'plot_text': '#172231',
        'plot_axis': '#2f3d4b',
        'plot_reference': 'rgba(47,61,75,0.32)',
        'plot_grid': 'rgba(47,61,75,0.06)',
        'plot_scatter_bg': 'rgba(126,139,151,0.23)',
    },
    'dark': {
        'bg': '#0f1720',
        'bg_secondary': '#131e29',
        'bg_card': '#182430',
        'text': '#edf3f8',
        'text_muted': '#9eafbf',
        'border': '#2a3a49',
        'accent': '#76a9da',
        'accent_hover': '#95bee5',
        'success': '#55c39a',
        'danger': '#ee8490',
        'shadow': '0 10px 28px rgba(0,0,0,0.22)',
        'shadow_hover': '0 14px 34px rgba(0,0,0,0.30)',
        'plot_bg': '#182430',
        'plot_text': '#edf3f8',
        'plot_axis': '#b6c3cf',
        'plot_reference': 'rgba(182,195,207,0.32)',
        'plot_grid': 'rgba(255,255,255,0.16)',
        'plot_scatter_bg': 'rgba(198,210,222,0.62)',
    }
}


def get_theme():
    return THEMES.get(st.session_state.get('theme'), THEMES['light'])


def inject_css():
    th = get_theme()
    st.markdown(f"""
    <style>
        html, body, [class*="css"] {{
            font-family: Inter, -apple-system, BlinkMacSystemFont, "Segoe UI",
                         "PingFang SC", "Microsoft YaHei", sans-serif;
        }}
        .stApp {{ background: {th['bg']}; color: {th['text']}; }}
        [data-testid="stMainBlockContainer"] {{
            max-width: 1280px; padding: 1.25rem 2rem 2rem;
        }}
        [data-testid="stVerticalBlock"] {{ gap: 0.85rem; }}
        section[data-testid="stSidebar"] {{
            background: {th['bg_secondary']}; border-right: 1px solid {th['border']};
        }}
        [data-testid="stSidebarUserContent"] {{ padding-top: 1.25rem; }}
        h1, h2, h3, h4, h5, label {{ color: {th['text']} !important; }}
        h2 {{ font-size: 1.15rem !important; }}
        h3 {{ font-size: 1.1rem !important; padding-bottom: 0.25rem !important; }}
        p, .stCaption {{ line-height: 1.5; }}
        [data-testid="stCaptionContainer"] {{ color: {th['text_muted']}; }}
        .hero-shell {{
            position: relative; overflow: hidden;
            display: flex; align-items: center; justify-content: space-between;
            gap: 1.25rem; padding: 1rem 1.25rem;
            background: {th['bg_card']}; border: 1px solid {th['border']};
            border-radius: 14px; box-shadow: 0 5px 18px rgba(31,54,78,0.06);
            margin-bottom: 0.25rem;
        }}
        .hero-shell::after {{
            content: ""; position: absolute; width: 14rem; height: 14rem;
            border-radius: 50%; right: -5rem; top: -9rem;
            background: {th['accent']}09; pointer-events: none;
        }}
        .hero-copy, .hero-version {{ position: relative; z-index: 1; }}
        .hero-kicker {{
            font-size: 0.64rem; font-weight: 650; letter-spacing: 0.12em;
            color: {th['text_muted']}; margin: 0 0 0.35rem;
        }}
        .main-header {{
            font-size: clamp(1.45rem, 2.5vw, 2rem) !important;
            font-weight: 740; letter-spacing: -0.035em; line-height: 1.15;
            margin: 0 !important; padding: 0 !important;
        }}
        .brand-accent {{ color: {th['accent']}; }}
        .sub-header {{ margin: 0.4rem 0 0; color: {th['text_muted']}; font-size: 0.88rem; }}
        .hero-version {{
            flex-shrink: 0; border: 1px solid {th['accent']}25; border-radius: 999px;
            color: {th['accent']}; background: {th['accent']}09;
            font-size: 0.73rem; font-weight: 650; padding: 0.42rem 0.75rem;
        }}
        .metrics-grid {{ display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 1rem; }}
        .metric-card {{
            position: relative; overflow: hidden; background: {th['bg_card']};
            border: 1px solid {th['border']}; border-radius: 12px;
            padding: 0.65rem 1rem; box-shadow: 0 3px 12px rgba(31,54,78,0.04);
        }}
        .metric-card::before {{
            content: ""; position: absolute; left: 0; top: 0; bottom: 0;
            width: 3px; background: {th['accent']}70;
        }}
        .metric-label {{ color: {th['text_muted']}; font-size: 0.76rem; margin-bottom: 0.2rem; }}
        .metric-value {{ color: {th['text']}; font-size: 1.6rem; font-weight: 720; line-height: 1.15; }}
        .scope-note {{ font-size: 0.78rem; color: {th['text_muted']}; padding: 0.15rem 0.1rem; overflow-wrap: anywhere; }}
        .st-key-rank_input, .st-key-box_input, .st-key-multi_input {{
            background: {th['bg_card']}; border: 1px solid {th['border']};
            border-radius: 12px; padding: 1rem; box-shadow: 0 4px 16px rgba(31,54,78,0.04);
        }}
        .status-strip {{
            padding: 0.65rem 0.85rem; border-left: 3px solid {th['success']};
            background: {th['bg_card']}; border-radius: 5px;
            color: {th['text_muted']}; font-size: 0.8rem; overflow-wrap: anywhere;
        }}
        .status-strip strong {{ color: {th['text']}; }}
        .gene-tag {{
            display: inline-block; color: {th['accent']}; background: {th['accent']}0d;
            border: 1px solid {th['accent']}24; border-radius: 5px;
            padding: 0.18rem 0.5rem; margin: 0.12rem 0.3rem 0.12rem 0;
            font-size: 0.82rem; font-weight: 600;
        }}
        [data-testid="stButtonGroup"] {{ margin: 0.1rem 0 0.2rem; }}
        [data-testid="stButtonGroup"] button {{ min-height: 2.45rem; font-weight: 600; }}
        [data-testid="stButtonGroup"] button[aria-checked="true"],
        [data-testid="stButtonGroup"] button[aria-pressed="true"] {{
            background: {th['accent']}14; color: {th['accent']}; border-color: {th['accent']};
        }}
        [data-testid="stVerticalBlockBorderWrapper"] {{
            border-color: {th['border']} !important; border-radius: 10px !important;
        }}
        [data-testid="stVerticalBlockBorderWrapper"] > div {{ background: {th['bg_card']}; }}
        [data-testid="stPlotlyChart"] {{
            background: {th['plot_bg']}; border: 1px solid {th['border']};
            border-radius: 10px; overflow: hidden;
        }}
        button {{ border-radius: 7px !important; }}
        .stButton button, .stDownloadButton button, [data-testid="stPopover"] button {{
            min-height: 2.4rem; font-weight: 600;
        }}
        button[kind="primary"] {{ background: {th['accent']}; border-color: {th['accent']}; color: {th['bg_card']}; }}
        button[kind="primary"]:hover {{ background: {th['accent_hover']}; border-color: {th['accent_hover']}; }}
        .stTextInput input, .stTextArea textarea, .stNumberInput input {{
            background: {th['bg_card']} !important; color: {th['text']} !important;
            -webkit-text-fill-color: {th['text']} !important;
            border-color: {th['border']} !important; line-height: 1.55;
        }}
        .stSelectbox [role="group"], .stMultiSelect [role="group"],
        .stSelectbox div[data-baseweb="select"] > div,
        .stMultiSelect div[data-baseweb="select"] > div {{
            background: {th['bg_card']} !important; color: {th['text']} !important;
            border: 1px solid {th['border']} !important; border-radius: 8px !important;
        }}
        .stSelectbox input[role="combobox"], .stMultiSelect input[role="combobox"] {{
            background: {th['bg_card']} !important; color: {th['text']} !important;
            -webkit-text-fill-color: {th['text']} !important; caret-color: {th['text']};
        }}
        .stSelectbox button[aria-label="Open"], .stMultiSelect button[aria-label="Open"] {{
            background: {th['bg_card']} !important; color: {th['text']} !important;
            border: none !important;
        }}
        .stSelectbox div[data-baseweb="select"] span {{ color: {th['text']}; }}
        div[data-baseweb="popover"], div[role="listbox"],
        div:has(> [role="listbox"]) {{
            background: {th['bg_card']} !important; color: {th['text']} !important;
            border: 1px solid {th['border']} !important;
        }}
        [role="option"] {{
            background: {th['bg_card']} !important; color: {th['text']} !important;
        }}
        [role="option"][aria-selected="true"], [role="option"][data-focused="true"] {{
            background: {th['accent']}33 !important;
        }}
        [data-testid="stExpander"] summary {{ background: {th['bg_card']}; font-weight: 550; }}
        [data-testid="stDataFrame"] {{ border-radius: 8px; overflow: hidden; border: 1px solid {th['border']}; }}
        .theme-grid-wrap {{ overflow-x: auto; border: 1px solid {th['border']}; border-radius: 8px; }}
        .theme-grid {{ width: 100%; border-collapse: collapse; background: {th['bg_card']}; color: {th['text']}; font-size: 0.92rem; }}
        .theme-grid th, .theme-grid td {{ padding: 0.48rem 0.75rem; border-bottom: 1px solid {th['border']}; text-align: right; white-space: nowrap; }}
        .theme-grid th:first-child, .theme-grid td:first-child {{ text-align: left; }}
        .theme-grid th {{ color: {th['text_muted']}; font-weight: 650; background: {th['bg_secondary']}; }}
        .theme-grid tr:last-child td {{ border-bottom: none; }}
        .footer-card {{ padding: 1rem; background: {th['bg_card']}; border: 1px solid {th['border']}; border-radius: 8px; }}
        .footer-card h4 {{ margin: 0 0 0.5rem; }}
        .footer-card p {{ color: {th['text_muted']}; margin: 0; font-size: 0.85rem; }}
        .footer-card a {{ color: {th['accent']}; }}
        .cite-box {{
            padding: 1rem; border-left: 3px solid {th['accent']}; background: {th['bg_card']};
            font-size: 0.85rem; overflow-wrap: anywhere;
        }}
        .doi-badge {{ font-size: 0.8rem; margin-left: 0.5rem; color: {th['accent']}; }}
        #MainMenu, footer {{ visibility: hidden; }}
        header[data-testid="stHeader"] {{ background: transparent; }}
        button:focus-visible, input:focus-visible, textarea:focus-visible {{
            outline: 3px solid {th['accent']}80 !important; outline-offset: 2px !important;
        }}
        :is(.stApp, [data-testid="stPopoverBody"]) button {{
            background-color: {th['bg_card']}; color: {th['text']} !important;
            border-color: {th['border']};
        }}
        :is(.stApp, [data-testid="stPopoverBody"]) button p, :is(.stApp, [data-testid="stPopoverBody"]) button span {{ color: inherit !important; }}
        :is(.stApp, [data-testid="stPopoverBody"]) button[kind="primary"] {{
            background-color: {th['accent']}; color: {th['bg_card']} !important;
            border-color: {th['accent']};
        }}
        :is(.stApp, [data-testid="stPopoverBody"]) button[kind="primary"]:hover {{ background-color: {th['accent_hover']}; }}
        .stApp [data-testid="stButtonGroup"] button[aria-checked="true"],
        .stApp [data-testid="stButtonGroup"] button[aria-pressed="true"] {{
            background-color: {th['accent']}22; color: {th['accent']} !important;
            border-color: {th['accent']};
        }}
        :is(.stApp, [data-testid="stPopoverBody"]) [data-testid="stExpander"] summary,
        :is(.stApp, [data-testid="stPopoverBody"]) [data-testid="stExpander"] summary p,
        :is(.stApp, [data-testid="stPopoverBody"]) [data-testid="stExpander"] summary span {{ color: {th['text']} !important; }}
        .stApp [data-baseweb="select"] > div {{
            background-color: {th['bg_card']} !important;
            color: {th['text']} !important; border-color: {th['border']} !important;
        }}
        .stApp [data-baseweb="select"] span,
        .stApp [data-baseweb="select"] input {{ color: {th['text']} !important; }}
        .stApp [data-testid="stAlert"],
        .stApp [data-testid="stAlert"] > div,
        .stApp [data-testid="stAlert"] [data-baseweb="notification"] {{
            background-color: {th['bg_secondary']} !important;
            color: {th['text']} !important;
            border: 1px solid {th['border']} !important;
            border-left: 3px solid {th['accent']} !important;
        }}
        .stApp [data-testid="stAlert"] p,
        .stApp [data-testid="stAlert"] span,
        .stApp [data-testid="stAlert"] svg {{
            color: {th['text']} !important;
            fill: {th['text']} !important;
        }}
        .theme-grid .score-low {{ color: #9ecbf2; font-weight: 700; }}
        .theme-grid .score-mid {{ color: #d5e4f2; font-weight: 650; }}
        .theme-grid .score-high {{ color: #f0c7b8; font-weight: 650; }}
        [data-testid="stPopoverBody"] {{
            background: {th['bg_card']} !important; color: {th['text']} !important;
            border: 1px solid {th['border']};
        }}
        [data-testid="stPopoverBody"] [data-testid="stCaptionContainer"],
        [data-testid="stPopoverBody"] [data-testid="stCaptionContainer"] p {{ color: {th['text_muted']} !important; }}
        @media (max-width: 900px) {{
            [data-testid="stMainBlockContainer"] {{ padding: 1rem 0.85rem 2rem; }}
            .hero-shell {{ align-items: flex-start; flex-wrap: wrap; gap: 0.65rem; padding: 0.9rem 1rem; }}
            .metrics-grid {{ gap: 0.65rem; }}
            .metric-value {{ font-size: 1.4rem; }}
        }}
        @media (prefers-reduced-motion: reduce) {{
            *, *::before, *::after {{ scroll-behavior: auto !important; transition: none !important; }}
        }}
    </style>
    """, unsafe_allow_html=True)

# =============================================================================
# 数据加载
# =============================================================================
def builtin_matrix_path() -> Path:
    override = os.environ.get("CRISPR_BUILTIN_MATRIX")
    return Path(override) if override else BUILTIN_MATRIX


@st.cache_resource(show_spinner=False)
def load_builtin_data(path: str):
    """Read the bundled matrix once per process. The cache key is the file path."""
    matrix = Path(path)
    provenance_path = builtin_provenance_path(matrix)
    if not matrix.is_file():
        return None, False, f"Built-in matrix not found: {matrix.name}"
    try:
        frame = read_builtin_matrix(matrix)
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return None, False, f"Cannot read the built-in matrix: {exc}"
    frame.attrs["sha256"] = file_sha256(matrix)
    frame.attrs["provenance"] = provenance
    return frame, True, None


@st.cache_resource(show_spinner=False, max_entries=3)
def load_uploaded_data(file_hash: str, _file_content: bytes):
    """按内容摘要缓存上传数据；返回的数据框在应用中只读使用。"""
    if not file_hash:
        raise ValueError("Missing upload content hash")
    return read_score_csv(_file_content)


@st.cache_data(show_spinner=False, max_entries=4)
def compute_gene_rankings(df_hash: str, _df: pd.DataFrame, gene_columns=None):
    return analyze_dataset(_df, gene_columns=gene_columns)


# =============================================================================
# 绘图配置
# =============================================================================
PLOT_COLORS = {
    # Okabe–Ito 色盲友好配色
    'essential': '#D55E00',
    'nonessential': '#009E73',
    'interest': '#0072B2',
    'threshold': '#9BA6B0',
    'boxplot_fill': '#56B4E9',
}
FONT_FAMILY = "Inter, Helvetica Neue, Arial, sans-serif"

PLOT_CONFIG = {
    'displaylogo': False,
    'modeBarButtonsToRemove': [
        'zoom2d', 'pan2d', 'select2d', 'lasso2d',
        'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'resetScale2d',
        'hoverClosestCartesian', 'hoverCompareCartesian', 'toggleSpikelines'
    ],
    'toImageButtonOptions': {
        'format': 'png', 'filename': 'crispr_plot',
        'height': 600, 'width': 1000, 'scale': 3
    }
}


def centered_plot(fig, config=None):
    """Center width-limited figures while allowing them to shrink on narrow screens."""
    render_config = dict(config or PLOT_CONFIG)
    render_config['toImageButtonOptions'] = {
        **render_config.get('toImageButtonOptions', {}),
        'width': int(fig.layout.width or 1000),
        'height': int(fig.layout.height or 600),
    }
    with st.container(horizontal_alignment='center'):
        st.plotly_chart(fig, width=fig.layout.width or 'stretch', config=render_config, theme=None)


def apply_theme_to_fig(fig):
    th = get_theme()
    fig.update_layout(
        plot_bgcolor=th['plot_bg'],
        paper_bgcolor=th['plot_bg'],
        font=dict(family=FONT_FAMILY, size=13, color=th['plot_text']),
        hoverlabel=dict(font=dict(family=FONT_FAMILY, size=13)),
        legend_font_color=th['plot_text'], title_font_color=th['plot_text'],
    )
    fig.update_xaxes(
        linecolor=th['plot_axis'], tickcolor=th['plot_axis'],
        tickfont=dict(size=12, color=th['plot_text']),
        title_font=dict(size=14, color=th['plot_text']),
    )
    fig.update_yaxes(
        linecolor=th['plot_axis'], tickcolor=th['plot_axis'],
        tickfont=dict(size=12, color=th['plot_text']),
        title_font=dict(size=14, color=th['plot_text']),
    )
    fig.update_annotations(font_color=th['plot_text'])
    return fig


# =============================================================================
# 绘图函数
# =============================================================================
def _rank_label_annotations(points, x_range, y_range):
    """Space labels in a conservative viewport, keeping arrows on the true data.

    A small spatial grid bounds collision checks even for a large submitted set.
    Positions use data coordinates so labels and leader lines also remain vector
    objects in exports and scale with the responsive chart.
    """
    plot_width, plot_height = 552, 366
    x_span, y_span = x_range[1] - x_range[0], y_range[1] - y_range[0]
    occupied = {}
    annotations = []

    def cells(rect):
        left, top, right, bottom = rect
        return [(col, row) for col in range(int(left // 64), int(right // 64) + 1)
                for row in range(int(top // 64), int(bottom // 64) + 1)]

    offsets = [(dx, dy) for distance in (24, 48, 72, 96, 120)
               for dx in (0, 48, -48, 96, -96) for dy in (-distance, distance)]
    for point in sorted(points, key=lambda item: (item['x'], item['y'], item['gene'])):
        px = (point['x'] - x_range[0]) / x_span * plot_width
        py = (y_range[1] - point['y']) / y_span * plot_height
        half_width = max(24, len(point['gene']) * 3.7 + 7)
        half_height = 12
        best = None
        for dx, dy in offsets:
            cx = min(plot_width - half_width - 4, max(half_width + 4, px + dx))
            cy = min(plot_height - half_height - 4, max(half_height + 4, py + dy))
            rect = (cx - half_width, cy - half_height, cx + half_width, cy + half_height)
            nearby = {other for cell in cells(rect) for other in occupied.get(cell, [])}
            overlap = sum(max(0, min(rect[2], other[2]) - max(rect[0], other[0]))
                          * max(0, min(rect[3], other[3]) - max(rect[1], other[1]))
                          for other in nearby)
            if best is None or overlap < best[0]:
                best = (overlap, cx, cy, rect)
            if overlap == 0:
                break
        _, cx, cy, rect = best
        for cell in cells(rect):
            occupied.setdefault(cell, []).append(rect)
        annotations.append(dict(
            x=point['x'], y=point['y'], xref='x', yref='y',
            ax=x_range[0] + cx / plot_width * x_span,
            ay=y_range[1] - cy / plot_height * y_span, axref='x', ayref='y',
            text=escape(point['gene']), showarrow=True, arrowhead=0,
            arrowwidth=0.8, arrowcolor=point['color'], standoff=5,
            xanchor='center', yanchor='middle', borderpad=2,
            bgcolor='rgba(0,0,0,0)', borderwidth=0,
            font=dict(size=13, color=point['color'], family=FONT_FAMILY),
        ))
    return annotations


def _ranking_fingerprint(rankings):
    """Identify a ranking curve without hashing every gene name."""
    if len(rankings) == 0:
        return (0, "", "", 0.0)
    scores = rankings["mean_score"].to_numpy(dtype=np.float64, copy=False)
    return (len(rankings), str(rankings["gene"].iloc[0]), str(rankings["gene"].iloc[-1]), round(float(scores.sum()), 5))


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
            y=rankings["mean_score"].to_numpy()[keep],
            mode="markers",
            marker=dict(size=2.5, color=color),
            name=name,
            text=rankings["gene"].to_numpy()[keep],
            customdata=np.column_stack((
                rankings["n_valid"].to_numpy()[keep],
                rankings["percentile"].to_numpy()[keep],
            )),
            hovertemplate=hover,
        )
        if len(cache) >= 8:
            cache.clear()
        cache[key] = cached
    return cached


def build_rank_figure(rankings, layers, references, n_cell_lines, point_size=4,
                      show_reference_labels=True):
    th = get_theme()
    fig = go.Figure()
    highlighted = {g.upper() for layer in layers for g in layer['genes']}
    refs = list(dict.fromkeys(g.upper() for g in references if g))
    hover = ('<b>%{text}</b><br>' + ui('Rank', '排名') + ': %{x:,}<br>'
             + ui('Mean score', '平均分') + ': %{y:.4f}<br>'
             + ui('Valid cell lines', '有效细胞系数') + ': %{customdata[0]}<br>'
             + ui('Rank percentile (lower = stronger)', '排名百分位（越低依赖越强）')
             + ': %{customdata[1]:.2f}%<extra></extra>')
    fig.add_trace(_cached_background_trace(
        rankings, highlighted | set(refs), th['plot_scatter_bg'], ui('All genes', '全部基因'), hover))
    label_points = []

    def collect_labels(part, color):
        label_points.extend({'x': record.rank, 'y': record.mean_score,
                             'gene': record.gene, 'color': color}
                            for record in part.itertuples(index=False))

    for i, gene in enumerate(refs):
        if gene in highlighted:
            continue
        part = rankings[rankings.gene_upper == gene]
        if part.empty:
            continue
        color = PLOT_COLORS['essential'] if i == 0 else PLOT_COLORS['nonessential']
        fig.add_trace(go.Scatter(x=part['rank'], y=part.mean_score, mode='markers',
                                 marker=dict(size=point_size * 2.6, color=color,
                                             symbol='diamond' if i == 0 else 'square',
                                             line=dict(width=1, color=th['plot_bg'])),
                                 text=part.gene, name=ui(f'Reference · {gene}', f'参考基因 · {gene}'),
                                 customdata=part[['n_valid', 'percentile']].values, hovertemplate=hover))
        if show_reference_labels:
            collect_labels(part, color)
    for layer in layers:
        part = rankings[rankings.gene_upper.isin({g.upper() for g in layer['genes']})].sort_values('rank')
        if part.empty:
            continue
        fig.add_trace(go.Scatter(x=part['rank'], y=part.mean_score,
                                 mode='markers',
                                 marker=dict(size=point_size * layer.get('size', 2.5), color=layer['color'],
                                             symbol=layer.get('symbol', 'circle'),
                                             line=dict(width=1, color=th['plot_bg'])),
                                 text=part.gene,
                                 name=layer['name'],
                                 customdata=part[['n_valid', 'percentile']].values, hovertemplate=hover))
        if layer['labels']:
            collect_labels(part, layer['color'])
    fig.add_hline(y=ESSENTIALITY_THRESHOLD, line=dict(dash='dash', color=PLOT_COLORS['threshold'], width=1))
    fig.add_hline(y=0, line=dict(color=th['plot_reference'], width=0.7))
    lo, hi = rankings.mean_score.min(), rankings.mean_score.max()
    span = max(hi - lo, 0.2)
    x_range = [0, max(len(rankings) * 1.02, 1)]
    y_range = [lo - span * 0.12, hi + span * 0.18]
    fig.update_layout(
        xaxis=dict(title=ui('Gene rank', '基因排名'), tickformat=',d', range=x_range,
                   showgrid=False, automargin=True),
        yaxis=dict(title=ui('Mean CRISPR score', '平均 CRISPR 分数'), range=y_range,
                   showgrid=True, gridcolor=th['plot_grid'], automargin=True),
        annotations=_rank_label_annotations(label_points, x_range, y_range),
        title=dict(text=ui(f'{n_cell_lines:,} cell lines · hover for gene-level sample counts',
                           f'{n_cell_lines:,} 个细胞系 · 悬停查看各基因有效样本数'),
                   font=dict(size=13), x=0, xanchor='left'),
        legend=dict(orientation='h', yanchor='top', y=-0.15, xanchor='center', x=0.5,
                    font=dict(size=12), itemsizing='constant'),
        width=640, height=480, margin=dict(l=64, r=24, t=42, b=72))
    return apply_theme_to_fig(fig)


def create_rank_plot(gene_rank_df, genes_of_interest, essential_gene='MYC', nonessential_gene='PTEN',
                     n_cell_lines=0, show_labels=True, point_size=4):
    return build_rank_figure(gene_rank_df,
                            [{'genes': genes_of_interest, 'labels': show_labels, 'color': PLOT_COLORS['interest'],
                              'name': ui('Genes of interest', '目标基因')}],
                            [essential_gene, nonessential_gene], n_cell_lines, point_size,
                            show_reference_labels=show_labels)


def create_multilayer_rank_plot(gene_rank_df, background_genes, highlight_genes, bg_color='#8DA8C0', hl_color='#0072B2',
                               essential_gene='MYC', nonessential_gene='PTEN', n_cell_lines=0, show_labels=True,
                               point_size=None):
    bg_only = [g for g in background_genes if g not in highlight_genes]
    layers = [
        {'genes': bg_only, 'labels': False, 'color': bg_color, 'size': 1.7,
         'name': ui('Background gene set', '背景基因集')},
        {'genes': highlight_genes, 'labels': show_labels, 'color': hl_color, 'symbol': 'circle',
         'name': ui('Highlight genes', '高亮基因')},
    ]
    if point_size is None:
        point_size = st.session_state.get('point_size', 4)
    return build_rank_figure(gene_rank_df, layers, [essential_gene, nonessential_gene], n_cell_lines,
                            point_size, show_reference_labels=show_labels)


def create_lineage_boxplot(lineage_data, genes, sort_mode='median', show_points=False):
    th = get_theme()
    # Horizontal groups keep long cancer-type names readable on narrow screens.
    genes = [g for g in genes if g in set(lineage_data.gene)]
    n_genes = len(genes)
    fig = make_subplots(rows=n_genes, cols=1, vertical_spacing=min(0.08, 0.35 / n_genes),
                        subplot_titles=genes)
    group_max = 1
    score_min, score_max = lineage_data.crispr_score.min(), lineage_data.crispr_score.max()
    pad = max((score_max - score_min) * 0.06, 0.15)
    for row, gene in enumerate(genes, 1):
        sub = lineage_data[lineage_data.gene == gene]
        summary = sub.groupby('lineage').crispr_score.agg(['median', 'count']).reset_index()
        summary = summary.sort_values(['median', 'lineage'] if sort_mode == 'median' else ['lineage'])
        group_max = max(group_max, len(summary))
        labels = []
        for record in summary.itertuples(index=False):
            label = f'{record.lineage} (n={record.count})'
            labels.append(label)
            points = sub[sub.lineage == record.lineage]
            fig.add_trace(go.Box(x=points.crispr_score, y=[label] * len(points), orientation='h',
                                 name=record.lineage, showlegend=False,
                                 marker=dict(color=PLOT_COLORS['interest'], size=4, opacity=0.5),
                                 line=dict(color=PLOT_COLORS['interest'], width=1.2),
                                 fillcolor='rgba(0,114,178,0.18)',
                                 boxpoints='all' if show_points else 'outliers', jitter=0.35, pointpos=0,
                                 customdata=points[['cell_line', 'cell_line_name']].values,
                                 hovertemplate=ui('Cell line', '细胞系') + ': %{customdata[1]} (%{customdata[0]})<br>'
                                               + ui('Score', '分数') + ': %{x:.4f}<extra>%{y}</extra>'), row=row, col=1)
        fig.update_yaxes(categoryorder='array', categoryarray=labels, autorange='reversed',
                         automargin=True, showgrid=False, row=row, col=1)
        fig.update_xaxes(title_text=ui('CRISPR score', 'CRISPR 分数'), range=[score_min - pad, score_max + pad],
                         showgrid=True, gridcolor=th['plot_grid'], row=row, col=1)
        fig.add_vline(x=0, line=dict(color=th['plot_reference'], width=0.7), row=row, col=1)
        fig.add_vline(x=ESSENTIALITY_THRESHOLD, line=dict(dash='dash', color=PLOT_COLORS['threshold'], width=1), row=row, col=1)
    fig.update_annotations(font=dict(size=15, family=FONT_FAMILY, color=th['plot_text']))
    fig.update_layout(height=max(300, group_max * 30 + 90) * n_genes + 50,
                      margin=dict(l=145, r=25, t=45, b=55), showlegend=False)
    return apply_theme_to_fig(fig)


# =============================================================================
# 图片导出（始终白底，方便论文用）
# =============================================================================
# =============================================================================
# Citation 渲染
# =============================================================================
def build_citations():
    """Separate the running app from the Zenodo deposit.

    The version DOI is the 1.0.0 snapshot from 2026-04-16. This checkout is 1.1.0.
    """
    running = (
        f"{TOOL_AUTHORS} ({TOOL_YEAR}). CRISPR Score Analyzer ({TOOL_VERSION}, {TOOL_RELEASED}) "
        f"[Software]. {GITHUB_URL}. Data: {DATA_VERSION}."
    )
    archived = (
        f"{TOOL_AUTHORS} ({TOOL_YEAR}). CRISPR Score Analyzer ({ZENODO_VERSION}) "
        f"[Software]. Zenodo. https://doi.org/{ZENODO_VERSION_DOI}"
    )
    bibtex = f"""@software{{kan_2026_crispr_score_analyzer,
  author  = {{{TOOL_AUTHORS}}},
  title   = {{CRISPR Score Analyzer: An interactive platform for DepMap gene essentiality}},
  year    = {{{TOOL_YEAR}}},
  version = {{{ZENODO_VERSION}}},
  doi     = {{{ZENODO_VERSION_DOI}}},
  url     = {{https://doi.org/{ZENODO_VERSION_DOI}}}
}}"""
    return running, archived, bibtex


def render_citation_section():
    running, archived, bibtex = build_citations()

    st.markdown(f"### {t('citation')}")
    st.markdown(f"**{ui('This deployment', '当前部署')}**")
    st.markdown(f'<div class="cite-box">{escape(running)}</div>', unsafe_allow_html=True)
    st.caption(ui(
        "Record this version and the data release with an analysis from the current app.",
        "记录当前网页上的分析时，写上这个版本和数据版本。"))
    st.markdown(f"**{ui('Archived software', '已归档版本')} (APA)**")
    st.markdown(
        f'<a href="https://doi.org/{ZENODO_VERSION_DOI}" target="_blank" class="doi-badge">'
        f'DOI: {ZENODO_VERSION_DOI}</a>',
        unsafe_allow_html=True)
    st.markdown(f'<div class="cite-box">{escape(archived)}</div>', unsafe_allow_html=True)
    st.caption(ui(
        f"Deposited {ZENODO_RELEASED}. Concept DOI {ZENODO_CONCEPT_DOI} still resolves to version {ZENODO_VERSION}. "
        "That Zenodo record is CC-BY-4.0. The code in this repository is MIT.",
        f"存档日期 {ZENODO_RELEASED}。概念 DOI {ZENODO_CONCEPT_DOI} 目前仍指向 {ZENODO_VERSION}。"
        "Zenodo 记录是 CC-BY-4.0，本仓库代码是 MIT。"))

    st.markdown("**BibTeX:**")
    st.code(bibtex, language="bibtex")

    st.markdown("**DepMap data citation:**")
    st.markdown(
        '<div class="cite-box">Tsherniak, A., et al. (2017). '
        'Defining a Cancer Dependency Map. <i>Cell</i> 170, 564–576. '
        'https://depmap.org/portal/</div>',
        unsafe_allow_html=True
    )


def ui(en, zh):
    return zh if st.session_state.get('lang') == 'zh' else en


def render_download_buttons(fig, filename_base, key_prefix, height=600):
    """Keep at most one current file per format, tied to its exact figure."""
    with st.expander(ui('Figure size', '图形尺寸'), expanded=False):
        c1, c2 = st.columns(2)
        width = c1.number_input(ui('Width (px)', '宽度（px）'), 400, 3000,
                                min(3000, max(400, int(fig.layout.width or 1000))), 50,
                                key=f'{key_prefix}_width')
        out_height = c2.number_input(ui('Height (px)', '高度（px）'), 300, 12000,
                                    min(12000, max(300, int(height))), 50,
                                    key=f'{key_prefix}_height')
    st.caption(ui(f'PDF / SVG: full vector. PNG: {width * 3} × {out_height * 3} px, 300 DPI.',
                  f'PDF / SVG：完整矢量。PNG：{width * 3} × {out_height * 3} 像素，300 DPI。'))
    for col, fmt, mime in zip(st.columns(3), ['pdf', 'png', 'svg'],
                              ['application/pdf', 'image/png', 'image/svg+xml']):
        cache_key = f'_export_{key_prefix}_{fmt}'
        signature = export_fingerprint(fig, fmt, width, out_height)
        cached = st.session_state.get(cache_key)
        if cached and cached['signature'] != signature:
            del st.session_state[cache_key]
            cached = None
        with col:
            if fmt == 'png' and width * out_height * 9 > MAX_PNG_PIXELS:
                st.caption(ui('PNG is too large. Reduce dimensions or use PDF/SVG.',
                              'PNG 尺寸过大，请缩小宽高或使用 PDF/SVG。'))
                continue
            if cached:
                st.download_button(f'↓ {fmt.upper()}', cached['bytes'],
                                   f'{filename_base}.{fmt}', mime,
                                   key=f'{key_prefix}_{fmt}_dl', on_click='ignore', width='stretch')
            elif st.button(fmt.upper(), help=ui(f'Generate {fmt.upper()}', f'生成 {fmt.upper()}'),
                           key=f'{key_prefix}_{fmt}_btn', width='stretch'):
                with st.spinner(ui('Generating figure…', '正在生成图片…')):
                    try:
                        image_bytes = render_figure_bytes(fig, fmt, width, out_height)
                        st.session_state[cache_key] = {'signature': signature, 'bytes': image_bytes}
                        st.rerun()
                    except Exception as exc:
                        st.error(ui('Export failed: ', '导出失败：') + str(exc))


def read_gene_input(prefix, default, method):
    if method == 'text':
        value = st.text_area(t('gene_list'), default, height=110, key=f'{prefix}_genes',
                             help=ui('Separate with newlines, spaces, tabs, commas or semicolons.',
                                     '支持换行、空格、制表符、中英文逗号和分号。'))
        return value, None
    file = st.file_uploader(t('input_file'), type=['csv', 'txt', 'tsv'], key=f'{prefix}_file',
                            help=ui('First column; a gene/symbol header is optional.',
                                    '读取第一列；可有 gene/symbol 表头，也可无表头。'))
    return '', file


def parse_input(text, file):
    return parse_gene_upload(file.getvalue(), file.name) if file else parse_gene_text(text)


def show_matches(parsed, matched, missing, title=None, resolved=None):
    resolved = list(resolved or [])
    alias_note = f" · alias resolved {len(resolved)}" if resolved else ""
    alias_note_zh = f" · 别名解析 {len(resolved)}" if resolved else ""
    if title:
        st.markdown(f'**{title}**')
    st.caption(ui(
        f"Input {parsed['input_count']} · unique {len(parsed['genes'])} · matched {len(matched)}{alias_note} · not matched {len(missing)} · duplicates {len(parsed['duplicates'])}",
        f"输入 {parsed['input_count']} · 去重后 {len(parsed['genes'])} · 匹配 {len(matched)}{alias_note_zh} · 未匹配 {len(missing)} · 重复 {len(parsed['duplicates'])}"))
    if resolved:
        st.caption(ui("Previous or alias symbols: ", "旧符号或别名：")
                   + ", ".join(f"{item['input']} → {item['gene']}" for item in resolved))
    if matched:
        st.markdown(' '.join(f'<span class="gene-tag">{escape(g)}</span>' for g in matched),
                    unsafe_allow_html=True)
    if missing or parsed['duplicates']:
        with st.expander(ui('Input details / unmatched genes', '输入详情 / 未匹配基因')):
            if missing:
                st.write(ui('Not matched (check gene symbols, dataset and score-column mapping):',
                            '未匹配（请检查基因名称、数据集和分数列映射）：'))
                st.code('\n'.join(missing), language=None)
            if parsed['duplicates']:
                st.write(ui('Duplicate entries removed:', '已去除的重复输入：'), ', '.join(parsed['duplicates']))


def build_result(genes, matches, config=None, plot_data=None, draft=None):
    return {
        'genes': genes,
        'draft': copy.deepcopy(draft or {}),
        'input_signature': input_signature(draft or {}),
        'rankings': gene_rankings,
        'matches': matches,
        'config': config or {},
        'plot_data': get_lineage_frame(df_scope, genes, gene_rankings, lineage_col) if plot_data is None else plot_data,
        'metadata': {**dataset_metadata, 'created_at_utc': datetime.now(timezone.utc).isoformat(),
                     'selected_genes': genes, 'parameters': config or {},
                     'input_matching': matches},
    }


def _score_tone(column, value):
    """Color a score by direction. Lower scores are the stronger dependency."""
    if column not in {"mean_score", "mean", "difference"} or pd.isna(value):
        return None
    value = float(value)
    if value < -0.5:
        return "score-low"
    if value < 0:
        return "score-mid"
    return "score-high"


def show_grid(frame, column_config):
    """Use Streamlit's grid in light mode, and a themed HTML table in dark mode.

    The grid draws its header on a canvas from the light Streamlit theme, so a
    Styler cannot recolor that header.
    """
    if st.session_state.get("theme") != "dark":
        st.dataframe(frame, hide_index=True, width="stretch", column_config=column_config)
        return

    def label_of(column):
        return column_config.get(column, {}).get("label", column)

    def text_of(column, value):
        if pd.isna(value):
            return ""
        spec = (column_config.get(column, {}).get("type_config") or {}).get("format")
        if spec == "percent":
            return escape(f"{float(value) * 100:.2f}%")
        if not spec:
            return escape(str(value))
        percent = spec.endswith("%%")
        core = spec[:-2] if percent else spec
        pattern = core[1:] if core.startswith("%") else core
        try:
            rendered = str(int(value)) if pattern in {"d", "i"} else format(value, pattern)
        except (TypeError, ValueError):
            return escape(str(value))
        rendered = escape(rendered + ("%" if percent else ""))
        tone = _score_tone(column, value)
        return f'<span class="{tone}">{rendered}</span>' if tone else rendered

    head = "".join(f"<th>{escape(str(label_of(column)))}</th>" for column in frame.columns)
    body = []
    for _, row in frame.iterrows():
        cells = "".join(f"<td>{text_of(column, row[column])}</td>" for column in frame.columns)
        body.append(f"<tr>{cells}</tr>")
    st.markdown(
        f'<div class="theme-grid-wrap"><table class="theme-grid"><thead><tr>{head}</tr></thead>'
        f"<tbody>{''.join(body)}</tbody></table></div>",
        unsafe_allow_html=True,
    )


def result_downloads(result, fig, prefix):
    selected = result['rankings'].set_index('gene', drop=False).reindex(result['genes']).reset_index(drop=True)
    state = st.session_state
    if "reference_a" in state and "reference_b" in state and "show_labels" in state and "point_size" in state:
        reference_a = str(state.get("reference_a") or "MYC").strip().upper()
        reference_b = str(state.get("reference_b") or "PTEN").strip().upper()
        labels_on = bool(state.get("show_labels"))
        marker_size = int(state.get("point_size") or 4)
    elif all(name in globals() for name in ("essential_gene", "nonessential_gene", "show_labels", "point_size")):
        reference_a, reference_b, labels_on, marker_size = essential_gene, nonessential_gene, show_labels, point_size
    else:
        reference_a, reference_b, labels_on, marker_size = "MYC", "PTEN", True, 4
    title_col, figure_col, data_col = st.columns([2, 1, 1.15], vertical_alignment='center')
    with title_col:
        st.markdown(f"### {ui('Results', '分析结果')}")
    with figure_col:
        with st.popover(ui('Export figure', '导出图形'), width='stretch'):
            render_download_buttons(fig, prefix, prefix, height=int(fig.layout.height or 600))
    with data_col:
        st.download_button(ui('Download results', '下载结果'),
                           selected.drop(columns=['gene_upper'], errors='ignore').to_csv(index=False),
                           f'{prefix}_selected_results.csv', 'text/csv', key=f'{prefix}_csv',
                           help=ui('CSV with original field names and full numeric precision.',
                                   'CSV 保留规范字段名与原始数值精度。'),
                           on_click='ignore', width='stretch')
    display_columns = ['gene', 'mean_score', 'rank', 'percentile', 'n_valid']
    table_config = {
        'gene': st.column_config.TextColumn(ui('Gene', '基因')),
        'mean_score': st.column_config.NumberColumn(ui('Mean score', '平均分'), format='%.3f'),
        'rank': st.column_config.NumberColumn(ui('Rank', '排名'), format='%d'),
        'percentile': st.column_config.NumberColumn(ui('Rank percentile', '排名百分位'), format='%.2f%%'),
        'n_valid': st.column_config.NumberColumn(ui('Valid samples', '有效样本数'), format='%d'),
        'n_missing': st.column_config.NumberColumn(ui('Missing samples', '缺失样本数'), format='%d'),
        'missing_fraction': st.column_config.NumberColumn(ui('Missing fraction', '缺失比例'), format='percent'),
    }
    show_grid(selected[display_columns], table_config)
    st.caption(ui('Lower rank percentile indicates stronger mean dependency. Valid samples count finite scores. CSV downloads retain the original field names and full precision.',
                  '排名百分位越低，平均依赖越强。有效样本数仅统计有限数值。CSV 下载保留规范字段名与原始精度。'))
    if prefix == 'box':
        st.caption(ui('These ranks and means use the sidebar cohort. Group statistics and plot data use the cancer types selected for this comparison.',
                      '上表排名与均值基于侧栏癌种范围；分组统计与绘图数据基于本次比较所选癌种。'))
    metadata = {**result['metadata'], 'display': {
        'language': st.session_state.lang, 'theme': st.session_state.theme,
        'reference_genes': [reference_a, reference_b],
        'show_labels': labels_on, 'point_size': marker_size,
    }, 'export': {'width': st.session_state.get(f'{prefix}_width', int(fig.layout.width or 1000)),
                  'height': st.session_state.get(f'{prefix}_height', int(fig.layout.height or 600)),
                  'png_scale': 3, 'png_dpi': 300, 'full_vector_pdf_svg': True}}
    with st.expander(ui('Details & reproducibility', '详细数据与复现'), expanded=False):
        st.markdown(f"**{ui('Sample completeness', '样本完整性')}**")
        show_grid(selected[['gene', 'n_valid', 'n_missing', 'missing_fraction']], table_config)
        if prefix == 'box':
            st.markdown(f"**{ui('Cancer-type summaries', '癌种分组统计')}**")
            group_config = {
                'gene': table_config['gene'],
                'lineage': st.column_config.TextColumn(ui('Cancer type', '癌种')),
                'n_valid': table_config['n_valid'],
                **{key: st.column_config.NumberColumn(label, format='%.3f') for key, label in [
                    ('median', ui('Median', '中位数')), ('mean', ui('Mean', '均值')),
                    ('std', ui('SD', '标准差')), ('min', ui('Minimum', '最小值')),
                    ('max', ui('Maximum', '最大值')),
                ]},
            }
            show_grid(lineage_summary(result['plot_data']), group_config)
        st.download_button(ui('Plot data CSV', '绘图数据 CSV'),
                           lambda frame=result['plot_data']: frame.to_csv(index=False).encode('utf-8-sig'),
                           f'{prefix}_plot_data.csv', 'text/csv', key=f'{prefix}_data', on_click='ignore')
        bundle_key = f'_export_{prefix}_zip'
        bundle_signature = hashlib.sha256(json.dumps(metadata, sort_keys=True, default=str).encode()).hexdigest()
        bundle = st.session_state.get(bundle_key)
        if bundle and bundle['signature'] != bundle_signature:
            del st.session_state[bundle_key]
            bundle = None
        if bundle is None and st.button(ui('Prepare analysis ZIP', '准备完整分析 ZIP'), key=f'{prefix}_zip_btn'):
            with st.spinner(ui('Preparing analysis bundle…', '正在准备分析文件…')):
                bundle = {'signature': bundle_signature,
                          'bytes': make_analysis_bundle(result['rankings'], result['genes'], result['plot_data'],
                                                        metadata, diagnostics)}
                st.session_state[bundle_key] = bundle
        if bundle:
            st.download_button(ui('Complete analysis ZIP', '完整分析 ZIP'), bundle['bytes'],
                               f'{prefix}_analysis.zip', 'application/zip', key=f'{prefix}_zip', on_click='ignore')
        st.caption(ui('ZIP includes selected results, complete cohort ranking, plot data, input matching, column diagnostics and versioned analysis parameters. Download figures separately.',
                      'ZIP 包含所选结果、当前范围完整排名、绘图数据、输入匹配情况、列校验与版本参数。图形请单独下载。'))
        st.json(metadata, expanded=False)


def classify_input(text, file):
    """Resolve the current draft against the active ranking, including HGNC aliases."""
    parsed = parse_input(text, file)
    matched, missing, resolved = resolve_gene_list(gene_rankings, parsed["genes"], GENE_ALIASES)
    return parsed, matched, missing, resolved


def render_draft_feedback(text, file, draft, result, title=None):
    """Show match status before Run, and again when the draft differs from the displayed result."""
    if result is not None and result.get("input_signature") == input_signature(draft):
        return
    if not text and file is None:
        return
    try:
        parsed, matched, missing, resolved = classify_input(text, file)
    except (ValueError, UnicodeError) as exc:
        st.caption(str(exc))
        return
    if not parsed["genes"] and not parsed["duplicates"]:
        return
    show_matches(parsed, matched, missing, title, resolved)


PREFERENCE_MIN_N = 5


def _preference_scale():
    """Blue is a lower score than the gene's own all-lineage mean."""
    neutral = "#243140" if st.session_state.get("theme") == "dark" else "#f7f4ee"
    return [[0.0, "#2f6fad"], [0.5, neutral], [1.0, "#d4785a"]]


def create_preference_plot(table):
    """One labeled lollipop per gene. Left of zero is more dependent.

    Each gene keeps its own scale, so a small shift is still readable next to
    a gene whose means move much further. The printed number is the difference.
    """
    shown = table[table["n_valid"] >= PREFERENCE_MIN_N]
    if shown.empty:
        return None
    genes = list(dict.fromkeys(shown["gene"]))
    longest = int(shown.groupby("gene").size().max())
    titles = []
    for gene in genes:
        mean = float(shown.loc[shown["gene"] == gene, "panel_mean"].iloc[0])
        titles.append(ui(f"{gene} · all-type mean {mean:.2f}", f"{gene} · 总平均 {mean:.2f}"))
    spacing = 0.12 if len(genes) < 3 else min(0.08, 0.85 / (len(genes) - 1))
    fig = make_subplots(
        rows=len(genes), cols=1, shared_xaxes=False, subplot_titles=titles,
        vertical_spacing=spacing,
    )
    for index, gene in enumerate(genes, start=1):
        _add_preference_lollipop(fig, shown[shown["gene"] == gene], gene, index)
    fig.update_layout(
        height=max(480, len(genes) * (26 * longest + 108)),
        showlegend=False,
        margin=dict(l=8, r=36, t=42, b=56),
    )
    themed = apply_theme_to_fig(fig)
    th = get_theme()
    themed.update_annotations(font=dict(size=15, family=FONT_FAMILY, color=th["plot_text"]))
    return themed


def _add_preference_lollipop(fig, gene_rows, gene, row):
    th = get_theme()
    part = gene_rows.sort_values(["difference", "lineage"], ascending=[False, True])
    order = part["lineage"].tolist()
    diffs = part["difference"].to_numpy(dtype=float)
    limit = float(np.nanmax(np.abs(diffs))) or 0.1
    stems = {"#2f6fad": ([], []), "#d4785a": ([], [])}
    for lineage, diff in zip(order, diffs):
        xs, ys = stems["#2f6fad" if diff < 0 else "#d4785a"]
        xs.extend((0, diff, None))
        ys.extend((lineage, lineage, None))
    for color, (xs, ys) in stems.items():
        if not xs:
            continue
        fig.add_trace(go.Scatter(
            x=xs, y=ys, mode="lines",
            line=dict(color=color, width=2),
            hoverinfo="skip", showlegend=False,
        ), row=row, col=1)
    fig.add_trace(go.Scatter(
        x=diffs, y=order, mode="markers+text",
        text=[f"{value:+.2f}" for value in diffs],
        textposition=["middle left" if value < 0 else "middle right" for value in diffs],
        textfont=dict(size=12, color=th["plot_text"]),
        cliponaxis=False,
        marker=dict(
            size=12, color=diffs, colorscale=_preference_scale(),
            cmid=0, cmin=-limit, cmax=limit, showscale=False, line=dict(width=0),
        ),
        customdata=np.column_stack((part["mean"], part["panel_mean"], part["n_valid"])),
        hovertemplate=(
            f"<b>%{{y}}</b> · {gene}<br>"
            + ui("Difference", "差值") + ": %{x:+.3f}<br>"
            + ui("Cancer-type mean", "该癌种均值") + ": %{customdata[0]:.3f}<br>"
            + ui("All-lineage mean", "全队列均值") + ": %{customdata[1]:.3f}<br>"
            + ui("Valid cell lines", "有效细胞系数") + ": %{customdata[2]:.0f}<extra></extra>"
        ),
        showlegend=False,
    ), row=row, col=1)
    fig.add_vline(x=0, line_width=1.6, line_color=th["plot_axis"], row=row, col=1)
    quiet_grid = "rgba(182,195,207,0.16)" if st.session_state.get("theme") == "dark" else "rgba(47,61,75,0.10)"
    fig.update_yaxes(
        automargin=True, categoryorder="array", categoryarray=order,
        showgrid=True, gridcolor=quiet_grid, gridwidth=1, row=row, col=1,
    )
    fig.update_xaxes(
        range=[-limit * 1.48, limit * 1.58], zeroline=False, showgrid=False,
        title_text=ui("← more dependent · less dependent →", "← 更依赖 · 更不依赖 →"),
        row=row, col=1,
    )


def render_lineage_preference(result):
    """Plot every cancer type against the all-lineage mean."""
    table = result.get("preference")
    if table is None or len(table) == 0:
        return
    st.markdown(f"### {ui('Lineage preference', '癌种偏好')}")
    st.caption(ui(
        "Each gene has its own list and its own scale. A point left of zero means that cancer type scores lower than "
        "this gene's average, so dependency is stronger there. The number is that difference. Stronger types are at the top. "
        f"Types with fewer than {PREFERENCE_MIN_N} cell lines are omitted. This describes the means. It is not a test of selectivity or significance.",
        "每个基因单独一条，刻度也各自独立。点在 0 的左边，表示该癌种分数低于这个基因的总平均，依赖更强。"
        f"数字就是这个差值，靠上的更明显。细胞系少于 {PREFERENCE_MIN_N} 的癌种不显示。这只是均值描述，不是选择性或显著性检验。"))
    figure = create_preference_plot(table)
    if figure is None:
        st.caption(ui("No cancer type has enough cell lines to plot.", "没有癌种达到可绘图的细胞系数。"))
        return
    st.plotly_chart(figure, config=PLOT_CONFIG, width="stretch", theme=None)


def attach_preference(result, genes):
    """Compare every cancer type with the all-lineage mean, independent of the sidebar filter."""
    if not lineage_col:
        return result
    observations = get_lineage_frame(df, genes, gene_rankings, lineage_col)
    table = lineage_preference(observations)
    if table.empty:
        return result
    result["preference"] = table
    result["metadata"]["lineage_preference"] = table.to_dict(orient="records")
    result["metadata"]["lineage_preference_definition"] = (
        "Per cancer type, mean score minus that gene's all-lineage mean. "
        "Negative means a lower score in that type. Descriptive only; not a test of selectivity."
    )
    return result


def figure_controls():
    """Chart-only settings. Kept beside the plot so changing them does not rebuild the page."""
    with st.expander(ui("Figure settings", "图形设置"), expanded=False):
        essential = st.text_input(t("essential"), value="MYC", key="reference_a").strip().upper()
        nonessential = st.text_input(t("nonessential"), value="PTEN", key="reference_b").strip().upper()
        labels = st.checkbox(t("show_labels"), value=True, key="show_labels")
        size = st.slider(t("point_size"), 2, 8, 4, key="point_size")
    return essential, nonessential, labels, size


@st.fragment
def render_curve_result(prefix):
    """Redraw labels and highlight markers without rerunning dataset loading."""
    result = st.session_state.get(f"{prefix}_result")
    if not result:
        return
    essential, nonessential, labels, size = figure_controls()
    n_cell_lines = result["metadata"].get("cohort_rows", 0)
    if prefix == "rank":
        fig = create_rank_plot(result["rankings"], result["genes"], essential, nonessential,
                               n_cell_lines, labels, size)
    else:
        cfg = result["config"]
        fig = create_multilayer_rank_plot(
            result["rankings"], cfg["background"], cfg["highlight"],
            cfg["background_color"], cfg["highlight_color"],
            essential, nonessential, n_cell_lines, labels, size)
    centered_plot(fig)
    render_lineage_preference(result)
    result_downloads(result, fig, prefix)


def input_draft(text='', file=None, **options):
    return {'text': text, 'file_name': file.name if file else None,
            'file_sha256': hashlib.sha256(file.getvalue()).hexdigest() if file else None,
            **options}


def input_signature(draft):
    return hashlib.sha256(json.dumps(draft, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def apply_input_preset(prefix, clear=False):
    presets = {'rank': {'genes': 'E2F1\nE2F2\nE2F3'}, 'box': {'genes': 'E2F1\nE2F2'},
               'multi': {'bg': 'CDK1\nCDK2\nCCNB1\nCCND1\nCCNE1', 'hl': 'PLK1\nAURKA'}}
    for key, value in presets[prefix].items():
        st.session_state[f'{prefix}_{key}'] = '' if clear else value
    if prefix != 'multi':
        st.session_state[f'{prefix}_method'] = 'text'
    if prefix == 'box':
        st.session_state.pop('box_candidates', None)


def input_actions(prefix):
    run_col, example_col, clear_col = st.columns([2, 1, 1])
    with run_col:
        submitted = st.button(ui('Run analysis', 'Run · 生成图表'), type='primary',
                              key=f'{prefix}_run', width='stretch')
    with example_col:
        st.button(ui('Load example', '载入示例'), key=f'{prefix}_example', width='stretch',
                  on_click=apply_input_preset, args=(prefix,))
    with clear_col:
        st.button(ui('Clear', '清空输入'), key=f'{prefix}_clear', width='stretch',
                  on_click=apply_input_preset, args=(prefix, True))
    return submitted


def render_run_status(prefix, draft, result):
    if result is None:
        st.info(ui('Ready to analyze. Enter genes or load an example, then Run.',
                   '准备就绪：输入基因或载入示例，点击 Run 生成结果。'))
        return
    if result.get('input_signature') != input_signature(draft):
        st.warning(ui('Inputs changed. Run to update the result.', '输入已修改，点击 Run 更新结果。'))
    meta = result['metadata']
    scope = ', '.join(meta['cohort_lineages'])
    if scope == 'ALL':
        scope = ui('All cancer types', '全部癌种')
    timestamp = datetime.fromisoformat(meta['created_at_utc']).strftime('%Y-%m-%d %H:%M:%S UTC')
    summary = ui(f"{len(result['genes'])} genes · {meta['cohort_rows']:,} cell lines · {scope}",
                 f"{len(result['genes'])} 个基因 · {meta['cohort_rows']:,} 个细胞系 · {scope}")
    st.markdown(f'<div class="status-strip"><strong>{ui("Displayed result", "当前显示结果")}</strong>'
                f' · {escape(summary)}<br>{escape(timestamp)}</div>', unsafe_allow_html=True)


def update_preferences(preference):
    options = TRANSLATIONS if preference == 'lang' else THEMES
    value = st.session_state.get(f'{preference}_select')
    if value in options:
        st.session_state[preference] = value


def render_preferences():
    with st.sidebar.expander(ui('Preferences', '界面设置')):
        st.selectbox('Language / 语言', ['en', 'zh'], key='lang_select',
                     index=['en', 'zh'].index(st.session_state.lang),
                     format_func={'en': 'English', 'zh': '中文'}.get,
                     on_change=update_preferences, args=('lang',))
        st.selectbox(t('theme'), ['light', 'dark'], key='theme_select',
                     index=['light', 'dark'].index(st.session_state.theme),
                     format_func={x: t(x) for x in ['light', 'dark']}.get,
                     on_change=update_preferences, args=('theme',))


for draft_key in ['rank_genes', 'box_genes', 'multi_bg', 'multi_hl', 'rank_method', 'box_method',
                  'multi_bg_color', 'multi_hl_color', 'box_lineages', 'box_order', 'box_points']:
    if draft_key in st.session_state:
        st.session_state[draft_key] = st.session_state[draft_key]

for preference, options in [('lang', TRANSLATIONS), ('theme', THEMES)]:
    if st.session_state.get(f'{preference}_select') not in options:
        st.session_state[f'{preference}_select'] = st.session_state[preference]
inject_css()
with st.sidebar:
    st.markdown(f"## {t('data_source')}")
    source_mode = st.segmented_control(ui('Dataset', '数据集'), ['builtin', 'custom'], default='builtin',
                                       format_func={'builtin': ui('Built-in', '内置数据'),
                                                    'custom': ui('Custom CSV', '自定义 CSV')}.get,
                                       key='dataset_source', selection_mode='single', width='stretch',
                                       label_visibility='collapsed') or 'builtin'
    uploaded_file = None
    if source_mode == 'custom':
        uploaded_file = st.file_uploader(t('upload_csv'), type=['csv'], key='dataset_file')
    else:
        st.caption(ui('Reference matrix included with the app', '参考矩阵随应用提供'))

crispr_data = None
source_sha = None
matrix_path = builtin_matrix_path()
source_id = ('upload:' + hashlib.sha256(uploaded_file.getvalue()).hexdigest()
             if uploaded_file is not None else ('builtin:' + str(matrix_path) if source_mode == 'builtin' else 'custom:empty'))
if st.session_state.get('_source_id') != source_id:
    for key in list(st.session_state):
        if key in {'rank_result', 'box_result', 'multi_result', 'box_candidates', 'box_display_genes', 'gene_columns', 'custom_mapping', 'box_lineages'} or key.startswith('_export_'):
            del st.session_state[key]
    st.session_state['_source_id'] = source_id
if source_mode == 'custom' and uploaded_file is None:
    render_preferences()
    st.markdown(f'# {t("app_title")}')
    st.info(ui('Upload a CRISPR score CSV in the sidebar to begin, or select Built-in to explore the reference dataset.',
               '请在侧栏上传 CRISPR 分数 CSV，或选择“内置数据”开始探索。'))
    st.stop()
try:
    if uploaded_file is not None:
        uploaded_bytes = uploaded_file.getvalue()
        source_sha = hashlib.sha256(uploaded_bytes).hexdigest()
        crispr_data = load_uploaded_data(source_sha, uploaded_bytes)
    else:
        with st.spinner(t('loading_hf')):
            crispr_data, success, err = load_builtin_data(str(matrix_path))
        if not success:
            st.error(err)
    if crispr_data is None or crispr_data.empty:
        render_preferences()
        st.error(ui('No data rows. Upload a nonempty score matrix.', '数据没有有效行，请上传非空分数矩阵。'))
        st.stop()
except (ValueError, UnicodeError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
    render_preferences()
    st.error(ui('Cannot read this CSV: ', '无法读取 CSV：') + str(exc))
    st.stop()

df = crispr_data
source_sha = source_sha or df.attrs.get('sha256')
provenance = df.attrs.get('provenance') or {}
lineage_col = find_lineage_column(df)
gene_columns = None
with st.sidebar:
    if uploaded_file is None:
        st.caption(f"{provenance.get('release', DATA_VERSION)} · {provenance.get('score_type', SCORE_TYPE)}")
    if uploaded_file is not None:
        with st.expander(ui('Gene column mapping', '基因列映射')):
            st.caption(ui('Automatically recognize SYMBOL (ENTREZ_ID). For plain symbols, explicitly select score columns. Metadata is never inferred from score ranges.',
                          '自动识别 SYMBOL (ENTREZ_ID)。仅有基因名时，请手动选择分数列；不会依据数值范围把元数据当作基因。'))
            if st.checkbox(ui('Select columns manually', '手动选择基因列'), key='custom_mapping'):
                gene_columns = st.multiselect(ui('Score columns', '分数列'), list(dict.fromkeys(df.columns)), key='gene_columns')
    st.markdown(ui('### Analysis scope', '### 分析范围'))
    if lineage_col:
        labels = df[lineage_col].fillna('Unknown').astype(str).str.strip().replace('', 'Unknown')
        all_lineages = sorted(labels.unique().tolist())
        if any(x not in all_lineages for x in st.session_state.get('cohort', [])):
            st.session_state.cohort = []
        cohort = st.multiselect(ui('Cancer types', '癌种'), all_lineages, key='cohort',
                                placeholder=ui('All cancer types', '全部癌种'),
                                help=ui('Leave empty to include all cancer types. Changing the scope clears previous results.',
                                        '留空表示全部癌种；修改分析范围后，需要重新运行。'))
        df_scope = df.loc[labels.isin(cohort)] if cohort else df
    else:
        cohort, all_lineages, df_scope = [], [], df
        st.caption(t('lineage_missing'))

context = json.dumps({'sha256': source_sha, 'cohort': sorted(cohort),
                      'gene_columns': gene_columns, 'source': uploaded_file.name if uploaded_file else matrix_path.name}, sort_keys=True)
context_hash = hashlib.sha256(context.encode()).hexdigest()
if st.session_state.get('_analysis_context') != context_hash:
    for key in list(st.session_state):
        if key in {'rank_result', 'box_result', 'multi_result', 'box_candidates', 'box_display_genes', 'box_lineages'} or key.startswith('_export_'):
            del st.session_state[key]
    st.session_state['_analysis_context'] = context_hash

gene_rankings, diagnostics = compute_gene_rankings(context_hash, df_scope, gene_columns)
n_cell_lines = len(df_scope)
hero_version = t('custom_dataset') if uploaded_file is not None else provenance.get('release', DATA_VERSION)
st.markdown(f'<section class="hero-shell"><div class="hero-copy">'
            f'<div class="hero-kicker">GENE DEPENDENCY EXPLORER</div>'
            f'<h1 class="main-header"><span class="brand-accent">CRISPR</span> Score Analyzer</h1>'
            f'<p class="sub-header">{t("app_subtitle")}</p></div>'
            f'<div class="hero-version">{escape(hero_version)}</div></section>', unsafe_allow_html=True)
scope_name = ', '.join(cohort) if cohort else ui('All cancer types', '全部癌种')
if gene_rankings.empty:
    render_preferences()
    st.warning(ui('No analyzable gene columns. For custom CSVs, select score columns in the sidebar. Gene symbols must be unique, with at least one finite numeric score per gene.',
                  '没有可分析的基因列。自定义 CSV 请在侧栏选择分数列；基因名不能重复，且每个基因至少需要一个有效数值分数。'))
    st.download_button(ui('Download validation report', '下载数据校验报告'),
                       lambda frame=diagnostics: frame.to_csv(index=False).encode('utf-8-sig'),
                       'column_diagnostics.csv', 'text/csv', key='validation_report', on_click='ignore')
    st.stop()

st.markdown(f'<div class="metrics-grid">'
            f'<div class="metric-card"><div class="metric-label">{t("cell_lines")}</div>'
            f'<div class="metric-value">{n_cell_lines:,}</div></div>'
            f'<div class="metric-card"><div class="metric-label">{t("gene_count")}</div>'
            f'<div class="metric-value">{len(gene_rankings):,}</div></div></div>'
            f'<div class="scope-note">{ui("Analysis scope", "分析范围")} · {escape(scope_name)}</div>',
            unsafe_allow_html=True)
with st.sidebar.expander(ui('Data overview', '数据概况')):
    st.metric(t('essential_genes').format(threshold=ESSENTIALITY_THRESHOLD),
              f"{(gene_rankings.mean_score < ESSENTIALITY_THRESHOLD).sum():,}")
    st.caption(ui('Mean-score range: ', '平均分范围：') +
               f'{gene_rankings.mean_score.min():.2f} – {gene_rankings.mean_score.max():.2f}')
    st.caption(ui('Lower score means stronger dependency. Rankings use means within the analysis scope; they do not establish cancer-type selectivity. The −0.5 line is a screening reference.',
                  '分数越低，依赖越强。排名使用分析范围内的平均分，本身不能证明癌种特异性。−0.5 是筛选参考线。'))
render_preferences()
source_files = (provenance.get('files') or {}) if uploaded_file is None else {}
dataset_metadata = {
    'data_source': uploaded_file.name if uploaded_file else 'bundled matrix',
    'data_file': uploaded_file.name if uploaded_file else matrix_path.name,
    'data_release': 'custom' if uploaded_file else provenance.get('release', DATA_VERSION),
    'data_revision': None if uploaded_file else (source_files.get('CRISPRGeneEffect.csv') or {}).get('sha256'),
    'data_sha256': source_sha, 'source_files': source_files or None,
    'cohort_lineages': cohort or ['ALL'],
    'lineage_column': lineage_col, 'gene_column_mapping': gene_columns,
    'cohort_rows': n_cell_lines, 'ranked_genes': len(gene_rankings),
    'mean_score_cutoff': ESSENTIALITY_THRESHOLD,
    'tool_version': TOOL_VERSION,
    'tool_released': TOOL_RELEASED,
    'zenodo_version': ZENODO_VERSION,
    'zenodo_version_doi': ZENODO_VERSION_DOI,
    'software': {name: version(name) for name in ['streamlit', 'pandas', 'numpy', 'plotly', 'kaleido', 'pillow']},
    'analysis_schema': 4,
}

view = st.segmented_control(ui('Analysis', '分析'), ['rank', 'box', 'multi'], default='rank',
                            format_func={'rank': t('tab1'), 'box': t('tab2'), 'multi': t('tab3')}.get,
                            key='analysis_tab', selection_mode='single', label_visibility='collapsed', width='stretch') or 'rank'
if view == 'rank':
    st.caption(ui('Locate genes in the dependency ranking. Lower mean scores indicate stronger dependency.',
                  '定位关注基因的依赖排名；平均分越低，依赖越强。'))
    with st.container(border=False, key='rank_input'):
        method = st.radio(t('input_method'), ['text', 'file'], horizontal=True, key='rank_method',
                          format_func={'text': t('input_direct'), 'file': t('input_file')}.get,
                          label_visibility='collapsed')
        gene_text, gene_file = read_gene_input('rank', 'E2F1\nE2F2\nE2F3', method)
        draft = input_draft(gene_text, gene_file, method=method)
        submitted = input_actions('rank')
    if submitted:
        try:
            parsed, matched, missing, resolved = classify_input(gene_text, gene_file)
            if matched:
                st.session_state.rank_result = attach_preference(build_result(
                    matched, {'targets': {'parsed': parsed, 'matched': matched, 'missing': missing, 'resolved': resolved}},
                    draft=draft), matched)
            else:
                st.warning(ui('No matched genes. The previous result is unchanged.', '没有匹配的基因，保留上次运行结果。'))
        except (ValueError, UnicodeError) as exc:
            st.error(str(exc))
    result = st.session_state.get('rank_result')
    render_run_status('rank', draft, result)
    render_draft_feedback(gene_text, gene_file, draft, result)
    if result:
        match = result['matches']['targets']
        show_matches(match['parsed'], match['matched'], match['missing'], resolved=match.get('resolved'))
        render_curve_result('rank')

elif view == 'box':
    st.caption(ui('Compare score distributions across cancer types for each gene.', '比较每个基因在不同癌种中的分数分布。'))
    if not lineage_col:
        st.info(t('lineage_missing'))
    else:
        with st.container(border=False, key='box_input'):
            method = st.radio(t('input_method'), ['text', 'file'], horizontal=True, key='box_method',
                              format_func={'text': t('input_direct'), 'file': t('input_file')}.get,
                              label_visibility='collapsed')
            gene_text, gene_file = read_gene_input('box', 'E2F1\nE2F2', method)
            with st.expander(ui('Comparison options', '比较选项')):
                box_lineages = st.multiselect(ui('Cancer types shown in the chart', '图中展示的癌种'),
                                              cohort or all_lineages, key='box_lineages',
                                              placeholder=ui('All types in the analysis scope', '分析范围内的全部癌种'),
                                              help=ui('This only limits the chart. Gene ranks still use the analysis scope in the sidebar.',
                                                      '此处仅限制图中分组；基因排名仍按侧栏的分析范围计算。'))
                sort_mode = st.selectbox(ui('Cancer-type order', '癌种顺序'), ['median', 'alphabetical'], key='box_order',
                                         format_func={'median': ui('Median score', '中位数'), 'alphabetical': ui('Alphabetical', '字母顺序')}.get)
                all_points = st.checkbox(ui('Show all cell lines', '显示全部细胞系散点'), key='box_points')
            draft = input_draft(gene_text, gene_file, method=method, lineages=box_lineages, sort=sort_mode, all_points=all_points)
            submitted = input_actions('box')
        render_draft_feedback(gene_text, gene_file, draft, st.session_state.get('box_result'))
        candidate = st.session_state.get('box_candidates')
        if candidate and candidate['input_signature'] != input_signature(draft):
            st.session_state.pop('box_candidates', None)
        if submitted:
            try:
                parsed, matched, missing, resolved = classify_input(gene_text, gene_file)
                match = {'targets': {'parsed': parsed, 'matched': matched, 'missing': missing, 'resolved': resolved}}
                config = {'lineages': box_lineages, 'sort': sort_mode, 'all_points': all_points}
                st.session_state.pop('box_candidates', None)
                if len(matched) > 8:
                    st.session_state.box_candidates = {'genes': matched, 'matches': match, 'config': config,
                                                       'draft': copy.deepcopy(draft), 'input_signature': input_signature(draft)}
                    st.session_state.pop('box_display_genes', None)
                elif matched:
                    frame = get_lineage_frame(df_scope, matched, gene_rankings, lineage_col)
                    if box_lineages:
                        frame = frame[frame.lineage.isin(box_lineages)]
                    st.session_state.box_result = attach_preference(build_result(matched, match, config, frame, draft), matched)
                else:
                    st.warning(ui('No matched genes. The previous result is unchanged.', '没有匹配的基因，保留上次运行结果。'))
            except (ValueError, UnicodeError) as exc:
                st.error(str(exc))
        candidate = st.session_state.get('box_candidates')
        if candidate:
            st.info(ui('More than 8 genes matched. Choose up to 8 to plot; the full input stays in the analysis metadata.',
                       '匹配超过 8 个基因，请选择最多 8 个绘图；完整输入仍保存在分析参数中。'))
            with st.form('box_selection'):
                chosen = st.multiselect(ui('Genes to display', '选择展示基因'), candidate['genes'],
                                        default=candidate['genes'][:8], max_selections=8, key='box_display_genes')
                apply_selection = st.form_submit_button(ui('Run selected genes', '运行所选基因'), key='box_selection_run', type='primary')
            if apply_selection and chosen:
                frame = get_lineage_frame(df_scope, chosen, gene_rankings, lineage_col)
                if candidate['config']['lineages']:
                    frame = frame[frame.lineage.isin(candidate['config']['lineages'])]
                st.session_state.box_result = attach_preference(
                    build_result(chosen, candidate['matches'], candidate['config'], frame, candidate['draft']), chosen)
                del st.session_state.box_candidates
                st.rerun()
        result = st.session_state.get('box_result')
        render_run_status('box', draft, result)
        if result:
            match = result['matches']['targets']
            show_matches(match['parsed'], match['matched'], match['missing'], resolved=match.get('resolved'))
            st.caption(ui('Displayed genes: ', '实际展示基因：') + ', '.join(result['genes']))
            valid_data = result['plot_data'].dropna(subset=['crispr_score'])
            if valid_data.empty:
                st.warning(ui('No finite scores in the selected groups.', '所选分组没有有效分数。'))
            else:
                absent = [g for g in result['genes'] if g not in set(valid_data.gene)]
                if absent:
                    st.info(ui('No valid scores in these groups: ', '这些基因在所选分组没有有效分数：') + ', '.join(absent))
                fig = create_lineage_boxplot(valid_data, result['genes'], result['config']['sort'], result['config']['all_points'])
                st.plotly_chart(fig, config=PLOT_CONFIG, width='stretch', theme=None)
                render_lineage_preference(result)
                result_downloads(result, fig, 'box')

elif view == 'multi':
    st.caption(ui('Mark a background gene set and highlight selected genes within the same ranking.',
                  '在同一排名中标注背景基因集，并突出显示关注基因。'))
    with st.container(border=False, key='multi_input'):
        c1, c2 = st.columns(2)
        with c1:
            bg_text = st.text_area(t('bg_gene_set'), 'CDK1\nCDK2\nCCNB1\nCCND1\nCCNE1', key='multi_bg', height=110)
        with c2:
            hl_text = st.text_area(t('hl_gene_set'), 'PLK1\nAURKA', key='multi_hl', height=110)
        with st.expander(ui('Annotation colors', '标注颜色')):
            color1, color2 = st.columns(2)
            with color1:
                bg_color = st.color_picker(t('bg_color'), '#8DA8C0', key='multi_bg_color')
            with color2:
                hl_color = st.color_picker(t('hl_color'), '#0072B2', key='multi_hl_color')
        draft = {'background': bg_text, 'highlight': hl_text, 'background_color': bg_color, 'highlight_color': hl_color}
        submitted = input_actions('multi')
    if submitted:
        bg_parsed, hl_parsed = parse_gene_text(bg_text), parse_gene_text(hl_text)
        bg_matched, bg_missing, bg_resolved = resolve_gene_list(gene_rankings, bg_parsed['genes'], GENE_ALIASES)
        hl_matched, hl_missing, hl_resolved = resolve_gene_list(gene_rankings, hl_parsed['genes'], GENE_ALIASES)
        matches = {'background': {'parsed': bg_parsed, 'matched': bg_matched, 'missing': bg_missing, 'resolved': bg_resolved},
                   'highlight': {'parsed': hl_parsed, 'matched': hl_matched, 'missing': hl_missing, 'resolved': hl_resolved}}
        if bg_matched or hl_matched:
            genes = list(dict.fromkeys(bg_matched + hl_matched))
            st.session_state.multi_result = attach_preference(build_result(
                genes, matches, {'background': bg_matched, 'highlight': hl_matched,
                                 'background_color': bg_color, 'highlight_color': hl_color}, draft=draft), genes)
        else:
            st.warning(ui('No matched genes. The previous result is unchanged.', '没有匹配的基因，保留上次运行结果。'))
    result = st.session_state.get('multi_result')
    render_run_status('multi', draft, result)
    render_draft_feedback(bg_text, None, draft, result, t('bg_gene_set'))
    render_draft_feedback(hl_text, None, draft, result, t('hl_gene_set'))
    if result:
        for key, title in [('background', t('bg_gene_set')), ('highlight', t('hl_gene_set'))]:
            match = result['matches'][key]
            show_matches(match['parsed'], match['matched'], match['missing'], title, match.get('resolved'))
        render_curve_result('multi')


# =============================================================================
# 页脚：致谢 + 引用
# =============================================================================
st.markdown("---")
with st.expander(f"📚 {t('resources')}", expanded=False):
    st.markdown(f"""
    <div class="footer-card">
        <h4>🙏 {t('acknowledgements')}</h4>
        <div style="display: flex; flex-wrap: wrap; gap: 1.5rem;">
            <div style="flex: 1; min-width: 250px;">
                <p>
                    <strong>{t('data_from')}</strong><br>
                    <a href="https://depmap.org" target="_blank">DepMap Portal (Broad Institute)</a><br>
                    <span style="font-size: 0.8rem;">CRISPR Chronos dependency scores</span>
                </p>
            </div>
            <div style="flex: 1; min-width: 250px;">
                <p>
                    <strong>{t('dev_with')}</strong><br>
                    <a href="https://www.anthropic.com/claude" target="_blank">Claude (Anthropic)</a><br>
                    <a href="https://openai.com/codex/" target="_blank">Codex (OpenAI)</a><br>
                    <span style="font-size: 0.8rem;">{t('ai_dev')}</span>
                </p>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)

    render_citation_section()

footer_release = f' · {escape(provenance.get("release", DATA_VERSION))}' if uploaded_file is None else ''
st.markdown(
    f'<div style="text-align:center; color:{get_theme()["text_muted"]}; '
    f'font-size:0.8rem; padding:1rem;">'
    f'CRISPR Score Analyzer {TOOL_VERSION}{footer_release} | Changhao Kan | '
    f'<a href="{GITHUB_URL}" target="_blank" style="color:{get_theme()["accent"]};">GitHub</a>'
    f'</div>',
    unsafe_allow_html=True
)
