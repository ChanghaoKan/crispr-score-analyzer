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
from crispr_score_analyzer.analysis_core import (
    analyze_dataset,
    builtin_provenance_path,
    default_builtin_matrix_path,
    file_sha256,
    filter_gene_suggestions,
    find_disease_column,
    find_lineage_column,
    gene_correlation_matrix,
    get_lineage_frame,
    lineage_preference,
    lineage_summary,
    load_gene_aliases as _load_package_aliases,
    make_analysis_bundle,
    parse_gene_text,
    parse_gene_upload,
    read_builtin_matrix,
    read_score_csv,
    resolve_gene_list,
    top_cell_lines,
)
from crispr_score_analyzer.codependency import (
    create_codependency_scan_plot,
    create_gene_pair_scatter,
    gene_pair_scatter_frame,
    genomewide_codependency,
)
from crispr_score_analyzer.figure_export import (
    MAX_PNG_PIXELS,
    export_fingerprint,
    render_figure_bytes,
)
from crispr_score_analyzer import plotting as _plotting
from crispr_score_analyzer.plotting import (
    ESSENTIALITY_THRESHOLD,
    FONT_FAMILY,
    PLOT_COLORS,
    PREFERENCE_MAX_LOLLIPOP_GENES,
    PREFERENCE_MIN_N,
)

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
BUILTIN_MATRIX = default_builtin_matrix_path()
GENE_ALIAS_FILE = Path(__file__).resolve().parent / "data" / "gene_aliases.json"
DATA_VERSION = "DepMap Public 26Q3"
SCORE_TYPE = "Chronos Gene Effect"

# =============================================================================
# Citation / DOI
# 1.2.0 的版本 DOI 是 10.5281/zenodo.24458535。概念 DOI 指向这一版。
# 2026-11-21 的 1.1.0 (10.5281/zenodo.24424595, MIT) 与 2026-04-16 的 1.0.0
# (10.5281/zenodo.19607603, CC-BY-4.0) 仍可按各自版本 DOI 引用。
# =============================================================================
TOOL_VERSION = "1.2.0"
TOOL_RELEASED = "2026-11-23"
TOOL_AUTHORS = "Kan, C."
TOOL_YEAR = "2026"
ZENODO_CONCEPT_DOI = "10.5281/zenodo.19607602"
ZENODO_VERSION_DOI = "10.5281/zenodo.24458535"
ZENODO_VERSION = "1.2.0"
ZENODO_RELEASED = "2026-11-23"
ZENODO_PREVIOUS_DOI = "10.5281/zenodo.19607603"
ZENODO_PREVIOUS_VERSION = "1.0.0"
ZENODO_PREVIOUS_RELEASED = "2026-04-16"
GITHUB_URL = "https://github.com/ChanghaoKan/crispr-score-analyzer"

# =============================================================================
# 国际化
# =============================================================================
TRANSLATIONS = {
    'en': {
        'app_title': 'CRISPR Score Analyzer',
        'app_subtitle': 'Gene dependency analysis and visualization across cancer cell lines',
        'language': 'Language',
        'theme': 'Theme',
        'light': 'Light',
        'dark': 'Dark',
        'data_source': 'Data Source',
        'upload_csv': 'Upload CRISPR Score CSV',
        'essential': 'Reference A',
        'nonessential': 'Reference B',
        'show_labels': 'Show gene labels',
        'point_size': 'Point size',
        'cell_lines': 'Cell lines',
        'gene_count': 'Genes',
        'essential_genes': 'Mean score < {threshold}',
        'resources': 'Citation & data source',
        'custom_dataset': 'Custom uploaded dataset',
        'lineage_missing': 'No lineage column was found in this dataset.',
        'tab1': 'Gene ranking',
        'tab2': 'Cancer-type comparison',
        'tab3': 'Gene-set annotation',
        'tab4': 'Co-dependency',
        'input_method': 'Input method',
        'input_direct': 'Direct input',
        'input_file': 'Upload gene list',
        'gene_list': 'Gene list',
        'bg_gene_set': 'Background gene set',
        'hl_gene_set': 'Highlight genes',
        'bg_color': 'Background color',
        'hl_color': 'Highlight color',
        'loading_hf': 'Loading the reference matrix...',
        'acknowledgements': 'Acknowledgements',
        'data_from': 'Data Source',
        'dev_with': 'Development Assistance',
        'ai_dev': 'AI assistance with code, interface design, testing and documentation',
        'citation': 'How to Cite',
    },
    'zh': {
        'app_title': 'CRISPR Score Analyzer',
        'app_subtitle': '癌症细胞系中的基因依赖分析与可视化',
        'language': '语言',
        'theme': '主题',
        'light': '浅色',
        'dark': '深色',
        'data_source': '数据来源',
        'upload_csv': '上传 CRISPR Score CSV',
        'essential': '参考基因 A',
        'nonessential': '参考基因 B',
        'show_labels': '显示基因名标签',
        'point_size': '点大小',
        'cell_lines': '细胞系',
        'gene_count': '基因数',
        'essential_genes': '平均分 < {threshold}',
        'resources': '引用与数据来源',
        'custom_dataset': '自定义上传数据',
        'lineage_missing': '该数据中未找到 lineage 列。',
        'tab1': '基因排名',
        'tab2': '癌种比较',
        'tab3': '基因集标注',
        'tab4': '基因共依赖',
        'input_method': '输入方式',
        'input_direct': '直接输入',
        'input_file': '上传基因列表',
        'gene_list': '基因列表',
        'bg_gene_set': '背景基因集',
        'hl_gene_set': '高亮基因',
        'bg_color': '背景颜色',
        'hl_color': '高亮颜色',
        'loading_hf': '正在读取参考矩阵...',
        'acknowledgements': '致谢',
        'data_from': '数据来源',
        'dev_with': '开发协助',
        'ai_dev': 'AI 辅助代码开发、界面设计、测试与文档整理',
        'citation': '引用方式',
    }
}


def t(key: str) -> str:
    lang = st.session_state.get('lang', 'en')
    return TRANSLATIONS[lang].get(key, key)


@st.cache_resource(show_spinner=False)
def load_gene_aliases(path: str) -> dict:
    """Uppercase previous/alias symbol to one approved symbol. Missing file yields no aliases."""
    return _load_package_aliases(path) or _load_package_aliases()


@st.cache_resource(show_spinner=False)
def get_software_versions() -> dict:
    return {name: version(name) for name in ['streamlit', 'pandas', 'numpy', 'plotly', 'kaleido', 'pillow']}


GENE_ALIASES = load_gene_aliases(str(GENE_ALIAS_FILE))


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
        'bg': '#f5f5f7',
        'bg_ambient': 'radial-gradient(ellipse 88% 320px at 50% -70px, rgba(0, 113, 227, 0.065), transparent 72%)',
        'bg_secondary': '#fbfbfd',
        'bg_card': '#ffffff',
        'bg_field': '#f5f5f7',
        'seg_track': 'rgba(118, 118, 128, 0.11)',
        'seg_active': '#ffffff',
        'glass_bg': 'rgba(255, 255, 255, 0.72)',
        'glass_bg_strong': 'rgba(255, 255, 255, 0.84)',
        'glass_sidebar': 'rgba(251, 251, 253, 0.78)',
        'glass_nav': 'rgba(240, 242, 247, 0.72)',
        'glass_nav_pill': 'rgba(255, 255, 255, 0.96)',
        'glass_th': 'rgba(251, 251, 253, 0.84)',
        'glass_rim': 'rgba(15, 23, 42, 0.085)',
        'glass_highlight': 'inset 0 1px 0 rgba(255, 255, 255, 0.92), inset 0 -0.5px 0 rgba(0, 0, 0, 0.04)',
        'glass_shadow': '0 8px 26px rgba(15, 23, 42, 0.065), 0 1.5px 4px rgba(15, 23, 42, 0.04)',
        'glass_popover_shadow': '0 18px 46px rgba(15, 23, 42, 0.13), 0 2px 8px rgba(15, 23, 42, 0.05)',
        'glass_pending_bg': 'rgba(255, 250, 240, 0.78)',
        'glass_pending_rim': 'rgba(255, 159, 10, 0.30)',
        'text': '#1d1d1f',
        'text_muted': '#6e6e73',
        'border': '#e5e5ea',
        'hairline': 'rgba(0, 0, 0, 0.06)',
        'accent': '#0071e3',
        'accent_hover': '#0077ed',
        'accent_text': '#ffffff',
        'success': '#248a3d',
        'warning': '#ff9f0a',
        'danger': '#d70015',
        'shadow': '0 1px 2px rgba(0, 0, 0, 0.03), 0 6px 20px rgba(0, 0, 0, 0.025)',
        'shadow_pill': '0 2px 6px rgba(0, 0, 0, 0.08), 0 1px 1.5px rgba(0, 0, 0, 0.04)',
        'plot_bg': '#ffffff',
        'plot_text': '#1d1d1f',
        'plot_axis': '#424245',
        'plot_reference': 'rgba(29,29,31,0.24)',
        'plot_grid': 'rgba(29,29,31,0.06)',
        'plot_scatter_bg': 'rgba(110,110,115,0.22)',
        'score_low': '#2166ac',
        'score_mid': '#3c5488',
        'score_high': '#b2182b',
    },
    'dark': {
        'bg': '#0f0f11',
        'bg_ambient': 'radial-gradient(ellipse 88% 360px at 50% -70px, rgba(10, 132, 255, 0.095), transparent 72%)',
        'bg_secondary': '#151518',
        'bg_card': '#1c1c1e',
        'bg_field': '#141417',
        'seg_track': 'rgba(118, 118, 128, 0.24)',
        'seg_active': '#2c2c2e',
        'glass_bg': 'rgba(30, 30, 34, 0.72)',
        'glass_bg_strong': 'rgba(34, 34, 38, 0.85)',
        'glass_sidebar': 'rgba(21, 21, 24, 0.80)',
        'glass_nav': 'rgba(34, 34, 39, 0.70)',
        'glass_nav_pill': 'rgba(255, 255, 255, 0.14)',
        'glass_th': 'rgba(25, 25, 28, 0.86)',
        'glass_rim': 'rgba(255, 255, 255, 0.12)',
        'glass_highlight': 'inset 0 1px 0 rgba(255, 255, 255, 0.18), inset 0 -0.5px 0 rgba(0, 0, 0, 0.42)',
        'glass_shadow': '0 10px 30px rgba(0, 0, 0, 0.50), 0 2px 6px rgba(0, 0, 0, 0.34)',
        'glass_popover_shadow': '0 22px 54px rgba(0, 0, 0, 0.66), 0 3px 10px rgba(0, 0, 0, 0.44)',
        'glass_pending_bg': 'rgba(48, 38, 18, 0.74)',
        'glass_pending_rim': 'rgba(255, 214, 10, 0.30)',
        'text': '#f5f5f7',
        'text_muted': '#98989d',
        'border': '#2c2c2e',
        'hairline': 'rgba(255, 255, 255, 0.09)',
        'accent': '#0a84ff',
        'accent_hover': '#409cff',
        'accent_text': '#ffffff',
        'success': '#30d158',
        'warning': '#ffd60a',
        'danger': '#ff453a',
        'shadow': '0 2px 10px rgba(0, 0, 0, 0.45)',
        'shadow_pill': '0 2px 8px rgba(0, 0, 0, 0.42)',
        'plot_bg': '#1c1c1e',
        'plot_text': '#f5f5f7',
        'plot_axis': '#aeaeb2',
        'plot_reference': 'rgba(245,245,247,0.24)',
        'plot_grid': 'rgba(255,255,255,0.07)',
        'plot_scatter_bg': 'rgba(174,174,178,0.34)',
        'score_low': '#67a9cf',
        'score_mid': '#92c5de',
        'score_high': '#e06b64',
    }
}


def get_theme():
    return THEMES.get(st.session_state.get('theme'), THEMES['light'])


def inject_css():
    th = get_theme()
    st.markdown(f"""
    <style>
        html, body, [class*="css"] {{
            font-family: -apple-system, BlinkMacSystemFont, "SF Pro Display", "SF Pro Text",
                         "Inter", "PingFang SC", "Hiragino Sans GB", "Helvetica Neue", sans-serif;
            -webkit-font-smoothing: antialiased;
            -moz-osx-font-smoothing: grayscale;
        }}
        .stApp {{
            background-color: {th['bg']};
            background-image: {th['bg_ambient']};
            background-attachment: fixed;
            color: {th['text']};
        }}
        [data-testid="stMainBlockContainer"] {{
            max-width: 1200px; padding: 1.05rem 2.2rem 2.6rem;
        }}
        [data-testid="stVerticalBlock"] {{ gap: 0.72rem; }}
        section[data-testid="stSidebar"] {{
            background: {th['glass_sidebar']};
            backdrop-filter: blur(22px) saturate(180%);
            -webkit-backdrop-filter: blur(22px) saturate(180%);
            border-right: 1px solid {th['glass_rim']};
        }}
        [data-testid="stSidebarUserContent"] {{ padding-top: 1rem; }}
        h1, h2, h3, h4, h5, label {{ color: {th['text']} !important; }}
        h2 {{ font-size: 1.04rem !important; font-weight: 620 !important; letter-spacing: -0.015em; }}
        h3 {{ font-size: 1.0rem !important; font-weight: 620 !important; letter-spacing: -0.015em; padding-bottom: 0.1rem !important; }}
        h4 {{ font-size: 0.94rem !important; font-weight: 600 !important; letter-spacing: -0.01em; }}
        p, .stCaption {{ line-height: 1.5; }}
        [data-testid="stCaptionContainer"],
        [data-testid="stCaptionContainer"] p {{
            color: {th['text_muted']} !important;
            opacity: 1 !important;
            font-size: 0.8rem;
        }}
        .hero-shell {{
            display: flex; align-items: flex-end; justify-content: space-between;
            flex-wrap: wrap; gap: 0.8rem 1.5rem;
            padding: 0.15rem 0 0.55rem; margin-bottom: 0.05rem;
        }}
        .hero-kicker {{
            font-size: 0.66rem; font-weight: 650; letter-spacing: 0.12em;
            text-transform: uppercase; color: {th['accent']}; margin: 0 0 0.22rem;
        }}
        .main-header {{
            font-size: clamp(1.5rem, 2.3vw, 1.92rem) !important;
            font-weight: 680; letter-spacing: -0.032em; line-height: 1.1;
            margin: 0 !important; padding: 0 !important;
            color: {th['text']} !important;
        }}
        .hero-shell [data-testid="stHeaderActionElements"],
        .hero-shell h1 a {{ display: none !important; }}
        .sub-header {{ margin: 0.28rem 0 0; color: {th['text_muted']}; font-size: 0.86rem; font-weight: 420; }}
        .hero-stats {{
            display: flex; align-items: center; flex-wrap: wrap; gap: 0.45rem 0.95rem;
            font-variant-numeric: tabular-nums;
        }}
        .hero-metric {{
            display: inline-flex; align-items: baseline; gap: 0.3rem;
            font-size: 0.78rem; color: {th['text_muted']};
        }}
        .hero-metric strong {{
            color: {th['text']}; font-size: 0.94rem; font-weight: 640;
        }}
        .hero-sep {{ color: {th['border']}; font-size: 0.85rem; user-select: none; }}
        .hero-version {{
            flex-shrink: 0; border: 1px solid {th['glass_rim']}; border-radius: 999px;
            color: {th['text_muted']}; background: {th['glass_bg']};
            backdrop-filter: blur(14px) saturate(170%);
            -webkit-backdrop-filter: blur(14px) saturate(170%);
            font-size: 0.72rem; font-weight: 550; padding: 0.24rem 0.72rem;
            box-shadow: {th['glass_highlight']}, 0 2px 8px rgba(0, 0, 0, 0.04);
        }}
        .scope-note {{
            font-size: 0.79rem; color: {th['accent']}; padding: 0.02rem 0.1rem;
            font-weight: 560; overflow-wrap: anywhere;
        }}
        .input-title {{
            font-size: 0.84rem; font-weight: 620; color: {th['text']};
            letter-spacing: -0.01em; line-height: 1.3;
        }}
        .section-gap {{ height: 0.45rem; }}
        .st-key-rank_input, .st-key-box_input, .st-key-multi_input, .st-key-codep_input {{
            background: {th['bg_card']}; border: 1px solid {th['hairline']};
            border-radius: 14px; padding: 0.95rem 1.15rem; box-shadow: {th['shadow']};
        }}
        .status-strip {{
            display: flex; align-items: center; gap: 0.5rem; flex-wrap: wrap;
            padding: 0.52rem 0.86rem; margin: 0.02rem 0;
            background: {th['glass_bg']};
            backdrop-filter: blur(16px) saturate(175%);
            -webkit-backdrop-filter: blur(16px) saturate(175%);
            border: 1px solid {th['glass_rim']};
            border-radius: 12px;
            box-shadow: {th['glass_highlight']}, 0 2px 10px rgba(0, 0, 0, 0.03);
            color: {th['text_muted']}; font-size: 0.8rem;
            font-variant-numeric: tabular-nums; overflow-wrap: anywhere;
            transition: background-color 0.2s ease, border-color 0.2s ease;
        }}
        .status-strip:has(.status-dot.pending) {{
            background: {th['glass_pending_bg']};
            border-color: {th['glass_pending_rim']};
        }}
        .status-dot {{
            width: 7px; height: 7px; border-radius: 50%; flex-shrink: 0;
            background: {th['success']};
        }}
        .status-dot.pending {{
            background: {th['warning']};
            box-shadow: 0 0 0 3px {th['warning']}2e;
        }}
        .status-strip strong {{ color: {th['text']}; font-weight: 590; }}
        .status-pending {{
            color: {th['warning']}; font-size: 0.75rem; font-weight: 600;
            letter-spacing: 0.005em;
        }}
        .status-time {{ opacity: 0.68; font-size: 0.73rem; margin-left: auto; }}
        .ready-note {{
            padding: 0.55rem 0.2rem; color: {th['text_muted']};
            font-size: 0.83rem; letter-spacing: 0.005em;
        }}
        /* Sticky floating Apple frosted glass navigation dock with spring-sliding glass thumb */
        [data-testid="stLayoutWrapper"]:has(> .st-key-analysis_nav) {{
            position: sticky;
            top: 0.75rem;
            z-index: 40;
        }}
        .st-key-analysis_nav [data-testid="stButtonGroup"] {{
            margin: 0 !important;
        }}
        .st-key-analysis_nav [data-testid="stButtonGroup"] > div {{
            position: relative !important;
            display: grid !important;
            grid-template-columns: repeat(4, minmax(0, 1fr)) !important;
            gap: 0 !important;
            background: {th['glass_nav']} !important;
            backdrop-filter: blur(24px) saturate(190%);
            -webkit-backdrop-filter: blur(24px) saturate(190%);
            border: 1px solid {th['glass_rim']};
            border-radius: 15px !important;
            padding: 4px !important;
            box-shadow: {th['glass_highlight']}, {th['glass_shadow']};
            overflow: hidden;
        }}
        .st-key-analysis_nav [data-testid="stButtonGroup"] > div::before {{
            content: "";
            position: absolute;
            top: 4px;
            bottom: 4px;
            left: 4px;
            width: calc((100% - 8px) / 4);
            border-radius: 11px;
            background: {th['glass_nav_pill']};
            border: 1px solid {th['glass_rim']};
            box-shadow: {th['glass_highlight']}, {th['shadow_pill']};
            pointer-events: none;
            z-index: 0;
            will-change: transform;
            transform: translate3d(0%, 0, 0);
            transition: transform 0.36s cubic-bezier(0.22, 1.1, 0.36, 1);
        }}
        .st-key-analysis_nav [data-testid="stButtonGroup"] > div:has(> button:nth-of-type(1)[aria-checked="true"])::before,
        .st-key-analysis_nav [data-testid="stButtonGroup"] > div:has(> button:nth-of-type(1)[aria-pressed="true"])::before {{
            transform: translate3d(0%, 0, 0);
        }}
        .st-key-analysis_nav [data-testid="stButtonGroup"] > div:has(> button:nth-of-type(2)[aria-checked="true"])::before,
        .st-key-analysis_nav [data-testid="stButtonGroup"] > div:has(> button:nth-of-type(2)[aria-pressed="true"])::before {{
            transform: translate3d(100%, 0, 0);
        }}
        .st-key-analysis_nav [data-testid="stButtonGroup"] > div:has(> button:nth-of-type(3)[aria-checked="true"])::before,
        .st-key-analysis_nav [data-testid="stButtonGroup"] > div:has(> button:nth-of-type(3)[aria-pressed="true"])::before {{
            transform: translate3d(200%, 0, 0);
        }}
        .st-key-analysis_nav [data-testid="stButtonGroup"] > div:has(> button:nth-of-type(4)[aria-checked="true"])::before,
        .st-key-analysis_nav [data-testid="stButtonGroup"] > div:has(> button:nth-of-type(4)[aria-pressed="true"])::before {{
            transform: translate3d(300%, 0, 0);
        }}
        [data-testid="stButtonGroup"] {{ margin: 0.05rem 0 0.1rem; }}
        [data-testid="stButtonGroup"] > div {{
            background: {th['seg_track']}; padding: 3px; border-radius: 11px;
            gap: 2px !important;
        }}
        [data-testid="stButtonGroup"] button {{
            min-height: 2.15rem; font-weight: 540; font-size: 0.84rem;
            letter-spacing: -0.008em; border-radius: 8px !important;
            border: 1px solid transparent !important;
            background-color: transparent !important;
            color: {th['text_muted']} !important;
            transition: background-color 0.16s ease, color 0.16s ease, box-shadow 0.16s ease;
        }}
        [data-testid="stButtonGroup"] button:hover {{
            color: {th['text']} !important;
        }}
        .stApp [data-testid="stButtonGroup"] button[aria-checked="true"],
        .stApp [data-testid="stButtonGroup"] button[aria-pressed="true"] {{
            background-color: {th['seg_active']} !important;
            color: {th['text']} !important;
            border-color: {th['hairline']} !important;
            box-shadow: {th['shadow_pill']};
            font-weight: 610;
        }}
        .stApp .st-key-analysis_nav [data-testid="stButtonGroup"] button {{
            position: relative;
            z-index: 1;
            min-height: 2.3rem;
            width: 100%;
            border-radius: 11px !important;
            border: 1px solid transparent !important;
            background-color: transparent !important;
            box-shadow: none !important;
            font-weight: 550;
            transition: color 0.2s ease, background-color 0.18s ease;
        }}
        .stApp .st-key-analysis_nav [data-testid="stButtonGroup"] button:hover:not([aria-checked="true"]):not([aria-pressed="true"]) {{
            background-color: {th['seg_track']} !important;
            color: {th['text']} !important;
        }}
        .stApp .st-key-analysis_nav [data-testid="stButtonGroup"] button[aria-checked="true"],
        .stApp .st-key-analysis_nav [data-testid="stButtonGroup"] button[aria-pressed="true"] {{
            background-color: transparent !important;
            border-color: transparent !important;
            box-shadow: none !important;
            color: {th['text']} !important;
            font-weight: 620;
        }}
        [data-testid="stVerticalBlockBorderWrapper"] {{
            border-color: {th['hairline']} !important; border-radius: 14px !important;
        }}
        [data-testid="stVerticalBlockBorderWrapper"] > div {{ background: {th['bg_card']}; }}
        [data-testid="stPlotlyChart"] {{
            background: {th['plot_bg']}; border: 1px solid {th['hairline']};
            border-radius: 14px; overflow: hidden; box-shadow: {th['shadow']};
        }}
        button {{
            border-radius: 9px !important;
            transition: border-color 0.16s ease, background-color 0.16s ease,
                        color 0.16s ease, box-shadow 0.18s ease, transform 0.16s ease;
        }}
        .stButton button, .stDownloadButton button, [data-testid="stPopover"] button {{
            min-height: 2.15rem; font-weight: 540; font-size: 0.83rem;
        }}
        button[kind="primary"] {{
            min-height: 2.38rem !important;
            background: {th['accent']} !important; border-color: {th['accent']} !important;
            color: {th['accent_text']} !important; font-weight: 600 !important;
            letter-spacing: -0.005em;
            box-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.28), 0 2px 8px {th['accent']}38 !important;
        }}
        button[kind="primary"]:hover {{
            background: {th['accent_hover']} !important; border-color: {th['accent_hover']} !important;
            box-shadow: inset 0 1px 0 rgba(255, 255, 255, 0.34), 0 4px 14px {th['accent']}48 !important;
            transform: translateY(-0.5px);
        }}
        .stTextArea [data-baseweb="textarea"],
        .stTextArea [data-baseweb="base-input"],
        .stTextInput [data-baseweb="input"],
        .stTextInput [data-baseweb="base-input"],
        .stNumberInput [data-baseweb="input"],
        .stNumberInput [data-baseweb="base-input"] {{
            background: {th['bg_field']} !important;
            border-color: {th['border']} !important;
            border-radius: 10px !important;
        }}
        .stTextArea [data-baseweb="textarea"]:focus-within,
        .stTextInput [data-baseweb="input"]:focus-within,
        .stNumberInput [data-baseweb="input"]:focus-within {{
            background: {th['bg_card']} !important;
            border-color: {th['accent']} !important;
            box-shadow: 0 0 0 3px {th['accent']}24 !important;
        }}
        .stTextInput input, .stTextArea textarea, .stNumberInput input {{
            background: {th['bg_field']} !important; color: {th['text']} !important;
            -webkit-text-fill-color: {th['text']} !important;
            border-color: transparent !important; border-radius: 10px !important;
            line-height: 1.55;
            font-family: "SFMono-Regular", ui-monospace, "JetBrains Mono", Menlo, monospace;
            font-size: 0.85rem !important;
        }}
        .stTextInput input:focus, .stTextArea textarea:focus, .stNumberInput input:focus {{
            background: {th['bg_card']} !important;
        }}
        [data-testid="stCheckbox"] [data-baseweb="checkbox"] > div:first-child {{
            background-color: {th['seg_track']};
            border: 1px solid {th['border']};
        }}
        .stSelectbox [role="group"], .stMultiSelect [role="group"],
        .stSelectbox div[data-baseweb="select"] > div,
        .stMultiSelect div[data-baseweb="select"] > div {{
            background: {th['bg_card']} !important; color: {th['text']} !important;
            border: 1px solid {th['border']} !important; border-radius: 10px !important;
        }}
        .stSelectbox input[role="combobox"], .stMultiSelect input[role="combobox"] {{
            background: transparent !important; color: {th['text']} !important;
            -webkit-text-fill-color: {th['text']} !important; caret-color: {th['text']};
        }}
        .stSelectbox button[aria-label="Open"], .stMultiSelect button[aria-label="Open"] {{
            background: transparent !important; color: {th['text']} !important;
            border: none !important;
        }}
        .stSelectbox div[data-baseweb="select"] span {{ color: {th['text']}; }}
        /* Frosted translucent dropdowns & popovers */
        div[data-baseweb="popover"], div[role="listbox"],
        div:has(> [role="listbox"]) {{
            background: {th['glass_bg_strong']} !important; color: {th['text']} !important;
            backdrop-filter: blur(24px) saturate(190%);
            -webkit-backdrop-filter: blur(24px) saturate(190%);
            border: 1px solid {th['glass_rim']} !important; border-radius: 14px !important;
            box-shadow: {th['glass_highlight']}, {th['glass_popover_shadow']} !important;
        }}
        [role="option"] {{
            background: transparent !important; color: {th['text']} !important;
            border-radius: 8px;
        }}
        [role="option"][aria-selected="true"], [role="option"][data-focused="true"] {{
            background: {th['accent']}22 !important;
        }}
        [data-testid="stFileUploaderDropzone"] {{
            background: {th['bg_field']} !important; color: {th['text']} !important;
            border-color: {th['border']} !important; border-radius: 10px !important;
        }}
        [data-testid="stCode"] pre {{
            background: {th['bg_field']} !important; color: {th['text']} !important;
            border: 1px solid {th['hairline']}; border-radius: 10px;
        }}
        .stMarkdown :not(pre) > code {{
            background: {th['bg_field']} !important; color: {th['accent']} !important;
            border: 1px solid {th['hairline']}; border-radius: 6px; padding: 0.08rem 0.36rem;
            font-family: "SFMono-Regular", ui-monospace, "JetBrains Mono", Menlo, monospace;
            font-size: 0.85em;
        }}
        [data-testid="stTooltipIcon"] svg {{
            color: {th['text_muted']} !important; stroke: {th['text_muted']} !important;
        }}
        [data-testid="stTabs"] [role="tab"],
        [data-testid="stTabs"] [role="tab"] p {{
            color: {th['text_muted']} !important; font-weight: 560; font-size: 0.86rem;
        }}
        [data-testid="stTabs"] [role="tab"]:hover,
        [data-testid="stTabs"] [role="tab"]:hover p {{
            color: {th['text']} !important;
        }}
        [data-testid="stTabs"] [role="tab"][aria-selected="true"],
        [data-testid="stTabs"] [role="tab"][aria-selected="true"] p {{
            color: {th['accent']} !important; font-weight: 620;
        }}
        [data-testid="stTabs"] [data-baseweb="tab-highlight"] {{
            background-color: {th['accent']} !important;
        }}
        [data-testid="stTabs"] [data-baseweb="tab-border"] {{
            background-color: {th['hairline']} !important;
        }}
        [data-testid="stExpander"] {{
            border-color: {th['hairline']} !important; border-radius: 12px !important;
            background: {th['bg_card']};
        }}
        [data-testid="stExpander"] summary {{ background: transparent; font-weight: 560; font-size: 0.85rem; }}
        [data-testid="stDataFrame"] {{ border-radius: 12px; overflow: hidden; border: 1px solid {th['hairline']}; }}
        .theme-grid-wrap {{
            overflow: auto; max-height: 420px;
            border: 1px solid {th['hairline']}; border-radius: 12px;
        }}
        .theme-grid {{
            width: 100%; border-collapse: collapse; background: {th['bg_card']};
            color: {th['text']}; font-size: 0.87rem; font-variant-numeric: tabular-nums;
        }}
        .theme-grid th, .theme-grid td {{ padding: 0.5rem 0.82rem; border-bottom: 1px solid {th['hairline']}; text-align: right; white-space: nowrap; }}
        .theme-grid th:first-child, .theme-grid td:first-child {{ text-align: left; font-weight: 560; }}
        .theme-grid th {{
            position: sticky; top: 0; z-index: 1;
            color: {th['text_muted']}; font-weight: 600; font-size: 0.77rem;
            letter-spacing: 0.015em;
            background: {th['glass_th']};
            backdrop-filter: blur(14px) saturate(165%);
            -webkit-backdrop-filter: blur(14px) saturate(165%);
        }}
        .theme-grid tbody tr:hover {{ background: {th['accent']}12; }}
        .theme-grid tr:last-child td {{ border-bottom: none; }}
        .footer-card {{ padding: 1.05rem 1.15rem; background: {th['bg_card']}; border: 1px solid {th['hairline']}; border-radius: 12px; }}
        .footer-card h4 {{ margin: 0 0 0.5rem; }}
        .footer-card p {{ color: {th['text_muted']}; margin: 0; font-size: 0.84rem; }}
        .footer-card a {{ color: {th['accent']}; text-decoration: none; }}
        .footer-card a:hover {{ text-decoration: underline; }}
        .cite-box {{
            padding: 0.95rem 1.05rem; border-left: 3px solid {th['accent']}; background: {th['bg_field']};
            border-radius: 0 10px 10px 0; font-size: 0.84rem; overflow-wrap: anywhere;
        }}
        .doi-badge {{ font-size: 0.8rem; margin-left: 0.5rem; color: {th['accent']}; text-decoration: none; }}
        #MainMenu, footer, [data-testid="stAppDeployButton"] {{ visibility: hidden; display: none; }}
        header[data-testid="stHeader"] {{ background: transparent; }}
        button:focus-visible, input:focus-visible, textarea:focus-visible {{
            outline: 2px solid {th['accent']}99 !important; outline-offset: 2px !important;
        }}
        :is(.stApp, [data-testid="stPopoverBody"]) button {{
            background-color: {th['bg_card']}; color: {th['text']} !important;
            border-color: {th['border']};
        }}
        :is(.stApp, [data-testid="stPopoverBody"]) button:hover {{
            border-color: {th['accent']}66;
        }}
        :is(.stApp, [data-testid="stPopoverBody"]) button p, :is(.stApp, [data-testid="stPopoverBody"]) button span {{ color: inherit !important; }}
        :is(.stApp, [data-testid="stPopoverBody"]) button[kind="primary"] {{
            background-color: {th['accent']} !important; color: {th['accent_text']} !important;
            border-color: {th['accent']} !important;
        }}
        :is(.stApp, [data-testid="stPopoverBody"]) button[kind="primary"]:hover {{ background-color: {th['accent_hover']} !important; }}
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
            background-color: {th['bg_card']} !important;
            color: {th['text']} !important;
            border: 1px solid {th['hairline']} !important;
            border-left: 3px solid {th['accent']} !important;
            border-radius: 10px !important;
        }}
        .stApp [data-testid="stAlert"] p,
        .stApp [data-testid="stAlert"] span,
        .stApp [data-testid="stAlert"] svg {{
            color: {th['text']} !important;
            fill: {th['text']} !important;
        }}
        .theme-grid .score-low {{ color: {th['score_low']}; font-weight: 660; }}
        .theme-grid .score-mid {{ color: {th['score_mid']}; font-weight: 580; }}
        .theme-grid .score-high {{ color: {th['score_high']}; font-weight: 620; }}
        [data-testid="stPopoverBody"] {{
            background: {th['glass_bg_strong']} !important; color: {th['text']} !important;
            backdrop-filter: blur(24px) saturate(190%);
            -webkit-backdrop-filter: blur(24px) saturate(190%);
            border: 1px solid {th['glass_rim']} !important; border-radius: 14px !important;
            box-shadow: {th['glass_highlight']}, {th['glass_popover_shadow']} !important;
        }}
        [data-testid="stPopoverBody"] button {{
            background-color: {th['glass_bg']} !important;
            border-color: {th['glass_rim']} !important;
        }}
        [data-testid="stPopoverBody"] button:hover {{
            background-color: {th['glass_nav_pill']} !important;
            border-color: {th['accent']}66 !important;
        }}
        [data-testid="stPopoverBody"] [data-testid="stCaptionContainer"],
        [data-testid="stPopoverBody"] [data-testid="stCaptionContainer"] p {{ color: {th['text_muted']} !important; }}
        @media (max-width: 900px) {{
            [data-testid="stMainBlockContainer"] {{ padding: 0.75rem 0.9rem 2rem; }}
            .hero-shell {{ align-items: flex-start; gap: 0.55rem; padding: 0.1rem 0 0.45rem; }}
            .status-time {{ margin-left: 0; }}
            [data-testid="stLayoutWrapper"]:has(> .st-key-analysis_nav) {{ top: 0.4rem; }}
        }}
        @media (max-width: 720px) {{
            .st-key-analysis_nav [data-testid="stButtonGroup"] > div {{
                grid-template-columns: repeat(2, minmax(0, 1fr)) !important;
                gap: 3px !important;
            }}
            .st-key-analysis_nav [data-testid="stButtonGroup"] > div::before {{
                display: none !important;
            }}
            .stApp .st-key-analysis_nav [data-testid="stButtonGroup"] button[aria-checked="true"],
            .stApp .st-key-analysis_nav [data-testid="stButtonGroup"] button[aria-pressed="true"] {{
                background-color: {th['glass_nav_pill']} !important;
                border-color: {th['glass_rim']} !important;
                box-shadow: {th['glass_highlight']}, {th['shadow_pill']} !important;
            }}
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
    frame.attrs["sha256"] = provenance.get("matrix_sha256") or file_sha256(matrix)
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


@st.cache_data(show_spinner=False, max_entries=8)
def compute_codependency_scan(df_hash: str, _df: pd.DataFrame, query_gene: str, _rankings: pd.DataFrame, method: str):
    return genomewide_codependency(_df, query_gene, _rankings, method=method)


# =============================================================================
# 绘图配置（核心绘图函数已抽离至 crispr_score_analyzer.plotting）
# =============================================================================
PLOT_CONFIG = {
    'displaylogo': False,
    'modeBarButtonsToRemove': [
        'select2d', 'lasso2d', 'autoScale2d',
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


def _active_plot_state():
    return st.session_state.get('theme', 'light'), st.session_state.get('lang', 'en')


def apply_theme_to_fig(fig):
    theme_name, _ = _active_plot_state()
    return _plotting.apply_theme_to_fig(fig, theme=theme_name)


def build_rank_figure(rankings, layers, references, n_cell_lines, point_size=4,
                      show_reference_labels=True, width=640, height=480):
    theme_name, lang = _active_plot_state()
    return _plotting.build_rank_figure(
        rankings, layers, references, n_cell_lines,
        point_size=point_size, show_reference_labels=show_reference_labels,
        width=width, height=height, theme=theme_name, lang=lang,
    )


def create_rank_plot(gene_rank_df, genes_of_interest, essential_gene='MYC', nonessential_gene='PTEN',
                     n_cell_lines=0, show_labels=True, point_size=4, width=640, height=480):
    theme_name, lang = _active_plot_state()
    return _plotting.create_rank_plot(
        gene_rank_df, genes_of_interest, essential_gene=essential_gene,
        nonessential_gene=nonessential_gene, n_cell_lines=n_cell_lines,
        show_labels=show_labels, point_size=point_size, width=width, height=height,
        theme=theme_name, lang=lang,
    )


def create_multilayer_rank_plot(gene_rank_df, background_genes, highlight_genes, bg_color='#8E99A2', hl_color='#3C5488',
                               essential_gene='MYC', nonessential_gene='PTEN', n_cell_lines=0, show_labels=True,
                               point_size=None, width=640, height=480):
    theme_name, lang = _active_plot_state()
    if point_size is None:
        point_size = st.session_state.get('point_size', 4)
    return _plotting.create_multilayer_rank_plot(
        gene_rank_df, background_genes, highlight_genes, bg_color=bg_color, hl_color=hl_color,
        essential_gene=essential_gene, nonessential_gene=nonessential_gene,
        n_cell_lines=n_cell_lines, show_labels=show_labels, point_size=point_size,
        width=width, height=height, theme=theme_name, lang=lang,
    )


def create_lineage_boxplot(lineage_data, genes, sort_mode='median', show_points=False,
                           columns=1, compare_mode='by_lineage'):
    theme_name, lang = _active_plot_state()
    return _plotting.create_lineage_boxplot(
        lineage_data, genes, sort_mode=sort_mode, show_points=show_points,
        columns=columns, compare_mode=compare_mode, theme=theme_name, lang=lang,
    )


# =============================================================================
# 图片导出（始终白底，方便论文用）
# =============================================================================
# =============================================================================
# Citation 渲染
# =============================================================================
def build_citations():
    """Cite the 1.2.0 Zenodo deposit. Keep the April 2026 snapshot addressable."""
    current = (
        f"{TOOL_AUTHORS} ({TOOL_YEAR}). CRISPR Score Analyzer ({ZENODO_VERSION}) "
        f"[Software]. Zenodo. https://doi.org/{ZENODO_VERSION_DOI}"
    )
    previous = (
        f"{TOOL_AUTHORS} ({TOOL_YEAR}). CRISPR Score Analyzer ({ZENODO_PREVIOUS_VERSION}) "
        f"[Software]. Zenodo. https://doi.org/{ZENODO_PREVIOUS_DOI}"
    )
    bibtex = f"""@software{{kan_2026_crispr_score_analyzer,
  author  = {{{TOOL_AUTHORS}}},
  title   = {{CRISPR Score Analyzer: an interactive, cohort-scoped browser platform for DepMap gene essentiality and co-dependency analysis}},
  year    = {{{TOOL_YEAR}}},
  version = {{{ZENODO_VERSION}}},
  doi     = {{{ZENODO_VERSION_DOI}}},
  url     = {{https://doi.org/{ZENODO_VERSION_DOI}}}
}}"""
    return current, previous, bibtex


def render_citation_section():
    current, previous, bibtex = build_citations()

    st.markdown(f"### {t('citation')}")
    st.markdown(f"**{ui('This version', '这一版')} (APA)**")
    st.markdown(
        f'<a href="https://doi.org/{ZENODO_VERSION_DOI}" target="_blank" class="doi-badge">'
        f'DOI: {ZENODO_VERSION_DOI}</a>',
        unsafe_allow_html=True)
    st.markdown(f'<div class="cite-box">{escape(current)}</div>', unsafe_allow_html=True)
    st.caption(ui(
        f"Deposited {ZENODO_RELEASED}. Concept DOI {ZENODO_CONCEPT_DOI} resolves to version {ZENODO_VERSION}. "
        f"This record and the repository code are MIT. Record {DATA_VERSION} with an analysis. "
        f"The {ZENODO_PREVIOUS_RELEASED} snapshot remains {previous}. Its Zenodo license is CC-BY-4.0.",
        f"存档日期 {ZENODO_RELEASED}。概念 DOI {ZENODO_CONCEPT_DOI} 指向 {ZENODO_VERSION}。"
        f"这条记录和仓库代码都是 MIT。分析请同时写上 {DATA_VERSION}。"
        f"{ZENODO_PREVIOUS_RELEASED} 的快照仍是 {previous}，当时的 Zenodo 许可证是 CC-BY-4.0。"))

    st.markdown("**BibTeX:**")
    st.code(bibtex, language="bibtex")

    st.markdown("**DepMap data citation & terms:**")
    st.markdown(
        '<div class="cite-box">Tsherniak, A., et al. (2017). '
        'Defining a Cancer Dependency Map. <i>Cell</i> 170, 564–576. '
        'DepMap Public 26Q3 (https://depmap.org/portal/data_page/).</div>',
        unsafe_allow_html=True
    )
    st.caption(ui(
        "The application code is MIT-licensed. Built-in DepMap Public 26Q3 data are provided for non-commercial "
        "research use under the Broad Institute DepMap Terms and Conditions (commercial use requires a separate "
        "license from Broad Institute or its contributors).",
        "本应用源码采用 MIT 许可证；内置 DepMap Public 26Q3 数据仅限非商业科研使用，受 Broad Institute DepMap "
        "使用条款约束（商业用途须另行获得 Broad Institute 或数据贡献方授权）。"))
    terms_path = Path(__file__).resolve().parent / "data" / "DEPMAP_TERMS.md"
    if terms_path.exists():
        with st.expander(ui("DepMap Terms and Conditions (Full Text)", "DepMap 数据使用条款（全文）"), expanded=False):
            st.markdown(terms_path.read_text(encoding="utf-8"))



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
        cached = st.session_state.get(cache_key)
        signature = None
        if cached:
            signature = export_fingerprint(fig, fmt, width, out_height)
            if cached['signature'] != signature:
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
                        if signature is None:
                            signature = export_fingerprint(fig, fmt, width, out_height)
                        image_bytes = render_figure_bytes(fig, fmt, width, out_height)
                        st.session_state[cache_key] = {'signature': signature, 'bytes': image_bytes}
                        st.rerun()
                    except Exception as exc:
                        st.error(ui('Export failed: ', '导出失败：') + str(exc))


def _apply_picked_genes(target_key, picker_key, query_key, single_select=False, fallback_gene=None):
    raw = st.session_state.get(picker_key)
    if isinstance(raw, str):
        chosen = [raw]
    else:
        chosen = list(raw or [])
    if not chosen and fallback_gene:
        chosen = [fallback_gene]
    if not chosen:
        return
    if single_select:
        st.session_state[target_key] = chosen[0]
    else:
        current = str(st.session_state.get(target_key) or '').strip()
        existing = parse_gene_text(current)['genes']
        additions = [g for g in chosen if g.upper() not in set(existing)]
        if additions:
            st.session_state[target_key] = (current + '\n' if current else '') + '\n'.join(additions)
    st.session_state[picker_key] = None if single_select else []
    st.session_state[query_key] = ''


def render_gene_quick_add(prefix, target_key, label=None, icon=':material/search:', single_select=False):
    """Native text-filtered symbol/alias search helper without fragile combobox hijacking."""
    if gene_rankings.empty:
        return
    picker_key = f'{prefix}_{target_key}_picker'
    query_key = f'{picker_key}_q'
    with st.popover(label or ui('Search', '搜索'), icon=icon):
        query = st.text_input(
            ui('Search symbols or HGNC aliases', '搜索基因符号或别名'),
            key=query_key,
            autocomplete='off',
            live=True,
            placeholder=ui('Type symbol or alias (e.g. TEAD1, SOX10, HER2)…',
                           '输入基因符号或别名（如 TEAD1、SOX10、HER2）…'),
        )
        suggestions = filter_gene_suggestions(query, all_gene_options, GENE_ALIASES, limit=24)
        fallback_gene = None
        if query and not suggestions:
            st.caption(ui('No matching symbols in the active dataset.', '当前数据集中无匹配基因。'))
        elif not query:
            st.caption(ui('Type at least 1 character to filter symbols or aliases.',
                          '输入字符即可实时筛选官方符号或 HGNC 别名。'))
        else:
            options = [item['gene'] for item in suggestions]
            via_map = {item['gene']: item['matched_via'] for item in suggestions if item.get('matched_via')}
            # Drop stale pill selections not present in the current filtered slice.
            cur_val = st.session_state.get(picker_key)
            if single_select:
                if cur_val not in options:
                    st.session_state[picker_key] = None
            elif cur_val:
                st.session_state[picker_key] = [g for g in cur_val if g in options]
            st.pills(
                ui('Matching genes', '匹配基因'),
                options,
                selection_mode='single' if single_select else 'multi',
                format_func=lambda g: f"{g} ({via_map[g]})" if g in via_map else g,
                key=picker_key,
                label_visibility='collapsed',
            )
            if len(options) == 1:
                fallback_gene = options[0]
        btn_label = ui('Use gene', '使用该基因') if single_select else ui('Add to list', '添加到列表')
        st.button(
            btn_label,
            key=f'{picker_key}_apply',
            on_click=_apply_picked_genes,
            args=(target_key, picker_key, query_key, single_select, fallback_gene),
            disabled=not suggestions,
            width='stretch',
        )


def read_gene_input(prefix, default, method, placeholder=None):
    if method == 'text':
        key = f'{prefix}_genes'
        st.session_state.setdefault(key, default)
        value = st.text_area(
            t('gene_list'), height=68, key=key,
            label_visibility='collapsed',
            placeholder=placeholder or ui(
                'Enter gene symbols (e.g. KRAS, SOX10, BRAF) — separated by newlines, spaces, or commas',
                '输入基因符号（如 KRAS、SOX10、BRAF）— 支持换行、空格或逗号分隔',
            ),
            help=ui('Separate with newlines, spaces, tabs, commas or semicolons. DepMap headers like "MYC (4609)" are also accepted.',
                    '支持换行、空格、制表符、中英文逗号和分号，也兼容直接粘贴 "MYC (4609)" 格式。'),
        )
        return value, None
    file = st.file_uploader(
        t('input_file'), type=['csv', 'txt', 'tsv'], key=f'{prefix}_file',
        label_visibility='collapsed',
        help=ui('First column; a gene/symbol header is optional.',
                '读取第一列；可有 gene/symbol 表头，也可无表头。'),
    )
    return '', file


def parse_input(text, file):
    return parse_gene_upload(file.getvalue(), file.name) if file else parse_gene_text(text)


def _has_match_issues(parsed, missing, resolved=None) -> bool:
    return bool(missing or (parsed or {}).get('duplicates') or resolved)


def show_matches(parsed, matched, missing, title=None, resolved=None):
    """Apple-minimalist exception-only feedback: stay silent when every gene matches cleanly."""
    resolved = list(resolved or [])
    duplicates = (parsed or {}).get('duplicates') or []
    if not missing and not duplicates and not resolved:
        return
    prefix_text = f"{title} · " if title else ""
    notes_en, notes_zh = [], []
    if resolved:
        notes_en.append(f"alias resolved {len(resolved)}")
        notes_zh.append(f"别名解析 {len(resolved)}")
    if missing:
        notes_en.append(f"not matched {len(missing)}")
        notes_zh.append(f"未匹配 {len(missing)}")
    if duplicates:
        notes_en.append(f"duplicates removed {len(duplicates)}")
        notes_zh.append(f"去重 {len(duplicates)}")
    st.caption(prefix_text + ui(" · ".join(notes_en), " · ".join(notes_zh)))
    if resolved:
        st.caption(ui("Previous or alias symbols: ", "旧符号或别名：")
                   + ", ".join(f"{item['input']} → {item['gene']}" for item in resolved))
    if missing or duplicates:
        with st.expander(ui('Unmatched / duplicate details', '未匹配 / 重复详情')):
            if missing:
                st.write(ui('Not matched (check gene symbols, dataset and score-column mapping):',
                            '未匹配（请检查基因名称、数据集和分数列映射）：'))
                st.code('\n'.join(missing), language=None)
            if duplicates:
                st.write(ui('Duplicate entries removed:', '已去除的重复输入：'), ', '.join(duplicates))


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
    if column not in {"mean_score", "mean", "difference", "crispr_score"} or pd.isna(value):
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
        tone = _score_tone(column, value) if "_score_tone" in globals() else None
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
        raw_a = state.get("reference_a")
        raw_b = state.get("reference_b")
        reference_a = ("MYC" if raw_a is None else str(raw_a)).strip().upper()
        reference_b = ("PTEN" if raw_b is None else str(raw_b)).strip().upper()
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
        with st.popover(ui('Export figure', '导出图形'), icon=':material/image:', width='stretch'):
            render_download_buttons(fig, prefix, prefix, height=int(fig.layout.height or 600))
    with data_col:
        st.download_button(ui('Download results', '下载结果'),
                           selected.drop(columns=['gene_upper'], errors='ignore').to_csv(index=False),
                           f'{prefix}_selected_results.csv', 'text/csv', key=f'{prefix}_csv',
                           icon=':material/download:',
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
        'reference_genes': [g for g in (reference_a, reference_b) if g],
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
        st.caption(ui('ZIP includes selected results, complete cohort ranking, plot data, lineage preference CSV (when available), input matching, column diagnostics and versioned analysis parameters. Download figures separately.',
                      'ZIP 包含所选结果、当前范围完整排名、绘图数据、癌种偏好 CSV（若有）、输入匹配情况、列校验与版本参数。图形请单独下载。'))
        preview_meta = dict(metadata)
        pref_rows = preview_meta.get('lineage_preference')
        if isinstance(pref_rows, list) and len(pref_rows) > 30:
            preview_meta['lineage_preference'] = f"<{len(pref_rows)} rows stored in lineage_preference.csv & metadata.json inside ZIP>"
        st.json(preview_meta, expanded=False)


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


def _preference_scale():
    theme_name, _ = _active_plot_state()
    return _plotting.preference_colorscale(theme=theme_name)


def create_preference_plot(table, highlight_lineages=None, top_k=None):
    theme_name, lang = _active_plot_state()
    return _plotting.create_preference_plot(
        table, highlight_lineages=highlight_lineages, top_k=top_k, theme=theme_name, lang=lang,
    )


def create_preference_heatmap(table, highlight_lineages=None):
    theme_name, lang = _active_plot_state()
    return _plotting.create_preference_heatmap(
        table, highlight_lineages=highlight_lineages, theme=theme_name, lang=lang,
    )


def create_correlation_heatmap(corr: pd.DataFrame, counts: pd.DataFrame, method: str = "pearson"):
    theme_name, lang = _active_plot_state()
    return _plotting.create_correlation_heatmap(
        corr, counts, method=method, theme=theme_name, lang=lang,
    )


def render_lineage_preference(result, prefix="pref"):
    """Plot cancer-type preference as a per-gene lollipop or compact heatmap, with export."""
    table = result.get("preference")
    if table is None or len(table) == 0:
        return
    shown_all = table[table["n_valid"] >= PREFERENCE_MIN_N]
    genes_all = list(dict.fromkeys(shown_all["gene"])) if not shown_all.empty else []
    active_cohort = [x for x in result["metadata"].get("cohort_lineages", []) if x != "ALL"]
    scope_note = (
        ui(" Active sidebar cancer types are outlined/starred; baseline always uses the full dataset.",
           " 侧栏所选癌种已描边/标星提示；基准均值始终基于全队列计算。")
        if active_cohort else ""
    )
    caption_text = ui(
        "Each gene is compared against its own all-lineage mean. Left of zero (blue) means stronger dependency in that cancer type. "
        f"Types with fewer than {PREFERENCE_MIN_N} cell lines are omitted. This describes the means. It is not a test of selectivity or significance.{scope_note}",
        "每个基因与其全癌种总均值比较：位于 0 左侧（蓝色）表示该癌种依赖更强。"
        f"细胞系少于 {PREFERENCE_MIN_N} 的癌种不显示。这只是均值描述，不是选择性或显著性检验。{scope_note}",
    )
    if shown_all.empty:
        st.markdown(f"### {ui('Lineage preference', '癌种偏好')}")
        st.caption(caption_text)
        st.caption(ui("No cancer type has enough cell lines to plot.", "没有癌种达到可绘图的细胞系数。"))
        return

    default_style = "heatmap" if len(genes_all) > PREFERENCE_MAX_LOLLIPOP_GENES else "lollipop"
    ctrl1, ctrl2, ctrl3 = st.columns([1.4, 1.4, 1.0], vertical_alignment="bottom")
    with ctrl1:
        pref_style = st.segmented_control(
            ui("Display style", "展示形式"),
            ["lollipop", "heatmap"],
            default=default_style,
            format_func={"lollipop": ui("Lollipop", "棒棒糖图"), "heatmap": ui("Heatmap", "热图")}.get,
            key=f"{prefix}_pref_style",
            selection_mode="single",
            label_visibility="collapsed",
            width="stretch",
        ) or default_style
    top_k = None
    plot_table = table
    with ctrl2:
        if pref_style == "lollipop":
            max_lineages = int(shown_all.groupby("gene").size().max())
            if max_lineages > 12:
                scope_choice = st.segmented_control(
                    ui("Lineages shown", "显示癌种数"),
                    ["top", "all"],
                    default="top",
                    format_func={"top": ui("Top & bottom 6", "仅看首尾各 6 个"), "all": ui("All types", "全部癌种")}.get,
                    key=f"{prefix}_pref_topk",
                    selection_mode="single",
                    label_visibility="collapsed",
                    width="stretch",
                ) or "top"
                top_k = 6 if scope_choice == "top" else None
    if pref_style == "lollipop" and len(genes_all) > PREFERENCE_MAX_LOLLIPOP_GENES:
        pref_genes_key = f"{prefix}_pref_genes"
        if any(g not in genes_all for g in st.session_state.get(pref_genes_key, [])):
            st.session_state.pop(pref_genes_key, None)
        picked = st.multiselect(
            ui(f"Select up to {PREFERENCE_MAX_LOLLIPOP_GENES} genes for lollipop view (or switch to Heatmap for all {len(genes_all)} genes)",
               f"棒棒糖图最多选 {PREFERENCE_MAX_LOLLIPOP_GENES} 个基因（查看全部 {len(genes_all)} 个基因请切换至热图）"),
            genes_all,
            default=genes_all[:PREFERENCE_MAX_LOLLIPOP_GENES],
            max_selections=PREFERENCE_MAX_LOLLIPOP_GENES,
            key=pref_genes_key,
        )
        plot_table = table[table["gene"].isin(picked)] if picked else table.iloc[:0]

    figure = (
        create_preference_heatmap(plot_table, highlight_lineages=active_cohort)
        if pref_style == "heatmap"
        else create_preference_plot(plot_table, highlight_lineages=active_cohort, top_k=top_k)
    )
    with ctrl3:
        if figure is not None:
            with st.popover(ui("Export figure", "导出图形"), icon=":material/image:", width="stretch"):
                render_download_buttons(
                    figure, f"{prefix}_lineage_preference", f"{prefix}_pref",
                    height=int(figure.layout.height or 600),
                )
    st.caption(caption_text)
    if figure is None:
        st.caption(ui("Select at least one gene with sufficient cell lines.", "请至少选择一个具备足够细胞系数的基因。"))
        return
    st.plotly_chart(figure, config=PLOT_CONFIG, width="stretch", theme=None)


def render_top_cell_lines(result, prefix):
    """Wet-lab oriented lookup of most dependent and least dependent cell lines."""
    plot_data = result.get("plot_data")
    if plot_data is None or plot_data.empty:
        st.caption(ui("No finite cell-line observations available in the current scope.", "当前范围内没有有效细胞系分数。"))
        return
    genes = [g for g in result["genes"] if g in set(plot_data["gene"])]
    if not genes:
        st.caption(ui("No finite cell-line observations available in the current scope.", "当前范围内没有有效细胞系分数。"))
        return
    top_gene_key = f"{prefix}_top_gene"
    if st.session_state.get(top_gene_key) not in genes:
        st.session_state.pop(top_gene_key, None)
    c1, c2, c3 = st.columns([2, 1, 1], vertical_alignment="bottom")
    with c1:
        gene = st.selectbox(ui("Gene", "基因"), genes, key=top_gene_key)
    with c2:
        top_n = st.selectbox(ui("Rows per side", "每侧展示数"), [5, 10, 20], index=1, key=f"{prefix}_top_n")
    summary = top_cell_lines(plot_data, [gene], top_n=top_n)
    with c3:
        st.download_button(
            ui("Download CSV", "下载 CSV"),
            lambda frame=summary: frame.to_csv(index=False).encode("utf-8-sig"),
            f"{prefix}_{gene}_top_cell_lines.csv", "text/csv",
            key=f"{prefix}_top_csv", icon=":material/download:", on_click="ignore", width="stretch",
        )
    col_cfg = {
        "rank_in_group": st.column_config.NumberColumn("#", format="%d"),
        "cell_line_name": st.column_config.TextColumn(ui("Cell line", "细胞系")),
        "cell_line": st.column_config.TextColumn("ModelID"),
        "lineage": st.column_config.TextColumn(ui("Cancer type", "癌种")),
        "crispr_score": st.column_config.NumberColumn(ui("Score", "分数"), format="%.3f"),
    }
    cols = ["rank_in_group", "cell_line_name", "cell_line", "lineage", "crispr_score"]
    left, right = st.columns(2)
    with left:
        st.markdown(f"**{ui('Most dependent (lowest scores)', '最敏感（分数最低）')}**")
        show_grid(summary[summary["direction"] == "most_dependent"][cols], col_cfg)
    with right:
        st.markdown(f"**{ui('Least dependent (highest scores)', '最耐受（分数最高）')}**")
        show_grid(summary[summary["direction"] == "least_dependent"][cols], col_cfg)
    st.caption(ui(
        "Useful for selecting positive (sensitive) and negative (resistant) cell-line models within the active analysis scope.",
        "便于在当前分析范围内挑选阳性敏感株与阴性对照株用于湿实验验证。",
    ))


def render_deep_dive(result, prefix):
    """Show lineage preference and top sensitive/resistant cell lines below primary results."""
    has_pref = result.get("preference") is not None and len(result["preference"]) > 0
    labels = []
    keys = []
    if has_pref:
        labels.append(ui("Lineage preference", "癌种偏好"))
        keys.append("pref")
    labels.append(ui("Top cell lines", "敏感 / 耐受细胞系"))
    keys.append("cells")
    tabs = st.tabs(labels)
    for tab, key in zip(tabs, keys):
        with tab:
            if key == "pref":
                render_lineage_preference(result, prefix=prefix)
            elif key == "cells":
                render_top_cell_lines(result, prefix=prefix)


def render_pair_scatter_panel(gene_x: str, gene_y: str, rankings: pd.DataFrame, key_prefix: str):
    """Render an interactive cell-line scatter plot comparing two genes, with exports."""
    scatter_df, stats = gene_pair_scatter_frame(df_scope, gene_x, gene_y, rankings, lineage_col)
    if scatter_df.empty:
        st.caption(ui("No shared finite cell-line observations for this gene pair.", "这两个基因在当前范围内没有共有的有效细胞系分数。"))
        return
    fig = create_gene_pair_scatter(
        scatter_df, stats, gene_x, gene_y,
        th=get_theme(), lang=st.session_state.get("lang", "en"),
    )
    h1, h2, h3 = st.columns([2.2, 1.0, 1.0], vertical_alignment="center")
    with h1:
        st.markdown(f"#### {ui('Cell-line scatter', '细胞系相关性散点图')} · `{gene_x}` vs `{gene_y}`")
    with h2:
        with st.popover(ui("Export figure", "导出图形"), icon=":material/image:", width="stretch"):
            render_download_buttons(
                fig, f"{key_prefix}_{gene_x}_vs_{gene_y}", f"{key_prefix}_scatter",
                height=int(fig.layout.height or 520),
            )
    with h3:
        st.download_button(
            ui("Scatter CSV", "散点数据 CSV"),
            lambda frame=scatter_df: frame.to_csv(index=False).encode("utf-8-sig"),
            f"{key_prefix}_{gene_x}_vs_{gene_y}.csv", "text/csv",
            key=f"{key_prefix}_scatter_csv", icon=":material/download:", on_click="ignore", width="stretch",
        )
    centered_plot(fig)


@st.fragment
def render_codependency_result():
    """Render either genome-wide single-gene co-dependency scan or multi-gene correlation matrix."""
    result = st.session_state.get("codep_result")
    if not result:
        return
    mode = result["mode"]
    method = result["method"]
    rankings = result["rankings"]
    metric_label = "Pearson r" if method == "pearson" else "Spearman ρ"

    if mode == "scan":
        query = result["query_gene"]
        scan_df = result["scan_df"]
        if scan_df.empty:
            st.warning(ui(
                "Not enough shared finite cell lines in the current scope to compute co-dependencies.",
                "当前分析范围内共有有效细胞系不足，无法计算共依赖相关性。",
            ))
            return
        st.session_state.setdefault("codep_scan_topn", 10)
        c1, c2, c3 = st.columns([1.4, 1.0, 1.1], vertical_alignment="bottom")
        with c1:
            top_n = st.selectbox(
                ui("Top genes per direction", "正/负相关各展示数"),
                [10, 15, 20], key="codep_scan_topn",
            )
        fig = create_codependency_scan_plot(
            scan_df, query, method=method, top_n=top_n,
            th=get_theme(), lang=st.session_state.get("lang", "en"),
        )
        with c2:
            if fig is not None:
                with st.popover(ui("Export figure", "导出图形"), icon=":material/image:", width="stretch"):
                    render_download_buttons(
                        fig, f"codep_{query}_{method}_top{top_n}", "codep_scan",
                        height=int(fig.layout.height or 600),
                    )
        with c3:
            st.download_button(
                ui("Download full scan CSV", "下载全基因组 CSV"),
                lambda frame=scan_df: frame.to_csv(index=False).encode("utf-8-sig"),
                f"codep_{query}_{method}_genomewide.csv", "text/csv",
                key="codep_scan_csv", icon=":material/download:", on_click="ignore", width="stretch",
            )
        st.caption(ui(
            f"Genome-wide {metric_label} scan against {query} across cell lines in the active sidebar scope. "
            "Correlations are descriptive and may reflect shared pathways, lineage specificity, or common-essentiality structure.",
            f"在当前侧栏分析范围内的细胞系中，将 {query} 与其余全部基因逐一计算 {metric_label}。"
            "相关性为描述性指标，可能受共同通路、癌种背景或普遍必需基因结构影响。",
        ))
        if fig is not None:
            st.plotly_chart(fig, config=PLOT_CONFIG, width="stretch", theme=None)

        pos_df = scan_df[scan_df["correlation"] > 0].head(50)
        neg_df = scan_df[scan_df["correlation"] < 0].tail(50).iloc[::-1]
        partner_options = list(dict.fromkeys(pos_df["gene"].head(20).tolist() + neg_df["gene"].head(20).tolist()))
        if not partner_options:
            partner_options = scan_df["gene"].head(20).tolist()
        if st.session_state.get("codep_scan_partner") not in partner_options:
            st.session_state["codep_scan_partner"] = partner_options[0]

        st.markdown('<div class="section-gap"></div>', unsafe_allow_html=True)
        p_col, m_col = st.columns([2.2, 1.1], vertical_alignment="bottom")
        with p_col:
            partner = st.selectbox(
                ui(f"Inspect cell-line scatter: {query} vs partner gene", f"查看细胞系散点图：{query} vs 共依赖伙伴基因"),
                partner_options, key="codep_scan_partner",
            )
        with m_col:
            top_partners = [query] + (pos_df["gene"].head(5).tolist() if not pos_df.empty else scan_df["gene"].head(5).tolist())
            if len(top_partners) >= 2 and st.button(
                ui(f"Matrix of top {len(top_partners)}", f"用 Top {len(top_partners)} 生成矩阵"),
                key="codep_scan_to_matrix", icon=":material/grid_on:", width="stretch",
            ):
                st.session_state["codep_mode"] = "matrix"
                st.session_state["codep_genes"] = ", ".join(top_partners)
                st.session_state["_autorun_view"] = "codep"
                st.rerun()
        render_pair_scatter_panel(query, partner, rankings, "codep_scan")

        st.markdown('<div class="section-gap"></div>', unsafe_allow_html=True)
        col_cfg = {
            "gene": st.column_config.TextColumn(ui("Gene", "基因")),
            "correlation": st.column_config.NumberColumn(metric_label, format="%.3f"),
            "n_shared": st.column_config.NumberColumn(ui("Shared cell lines", "共有细胞系"), format="%d"),
            "mean_score": st.column_config.NumberColumn(ui("Mean score", "平均分"), format="%.3f"),
            "rank": st.column_config.NumberColumn(ui("Cohort rank", "范围排名"), format="%d"),
        }
        cols = ["gene", "correlation", "n_shared", "mean_score", "rank"]
        t_pos, t_neg = st.tabs([
            ui("Top positive (co-dependent)", "正相关 Top 50（共依赖）"),
            ui("Top negative (anti-correlated)", "负相关 Top 50（反向依赖）"),
        ])
        with t_pos:
            if pos_df.empty:
                st.caption(ui("No positively correlated genes found in the active scope.", "当前范围内无正相关基因。"))
            else:
                show_grid(pos_df[cols], col_cfg)
        with t_neg:
            if neg_df.empty:
                st.caption(ui("No negatively correlated genes found in the active scope.", "当前范围内无负相关基因。"))
            else:
                show_grid(neg_df[cols], col_cfg)
        return

    # Matrix mode (2-25 genes)
    genes = result["genes"]
    corr = result["corr"]
    counts = result["counts"]
    fig = create_correlation_heatmap(corr, counts, method=method)
    h1, h2, h3 = st.columns([2.0, 1.0, 1.1], vertical_alignment="center")
    with h1:
        st.markdown(f"### {ui('Pairwise co-dependency matrix', '多基因共依赖相关性矩阵')} · {metric_label}")
    with h2:
        with st.popover(ui("Export figure", "导出图形"), icon=":material/image:", width="stretch"):
            render_download_buttons(
                fig, f"codep_matrix_{method}", "codep_matrix",
                height=int(fig.layout.height or 500),
            )
    with h3:
        st.download_button(
            ui("Matrix CSV", "相关性矩阵 CSV"),
            lambda frame=corr: frame.to_csv().encode("utf-8-sig"),
            f"codep_matrix_{method}.csv", "text/csv",
            key="codep_matrix_csv", icon=":material/download:", on_click="ignore", width="stretch",
        )
    st.caption(ui(
        "Pairwise dependency correlation across cell lines in the active sidebar scope. High positive correlation often reflects shared complexes or pathways; correlations are descriptive.",
        "基于当前侧栏分析范围内的细胞系计算两两基因依赖相关性；高正相关常提示共复合体或同通路功能（描述性统计）。",
    ))
    centered_plot(fig)

    # Pick strongest non-diagonal pair as default for the scatter plot.
    default_x, default_y = genes[0], genes[1]
    best_abs = -1.0
    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            val = corr.iloc[i, j]
            if pd.notna(val) and abs(float(val)) > best_abs:
                best_abs = abs(float(val))
                default_x, default_y = genes[i], genes[j]
    if st.session_state.get("codep_pair_x") not in genes:
        st.session_state["codep_pair_x"] = default_x
    if st.session_state.get("codep_pair_y") not in genes or st.session_state.get("codep_pair_y") == st.session_state.get("codep_pair_x"):
        st.session_state["codep_pair_y"] = default_y if default_y != st.session_state["codep_pair_x"] else next(
            g for g in genes if g != st.session_state["codep_pair_x"]
        )

    st.markdown('<div class="section-gap"></div>', unsafe_allow_html=True)
    px1, px2 = st.columns(2)
    with px1:
        gene_x = st.selectbox(ui("Gene X (horizontal axis)", "基因 X（横轴）"), genes, key="codep_pair_x")
    y_options = [g for g in genes if g != gene_x]
    if st.session_state.get("codep_pair_y") not in y_options:
        st.session_state["codep_pair_y"] = y_options[0]
    with px2:
        gene_y = st.selectbox(ui("Gene Y (vertical axis)", "基因 Y（纵轴）"), y_options, key="codep_pair_y")
    render_pair_scatter_panel(gene_x, gene_y, rankings, "codep_matrix")


def attach_preference(result, genes):
    """Compare every cancer type with the all-lineage mean, independent of the sidebar filter."""
    if not lineage_col or not genes:
        return result
    observations = get_lineage_frame(df, genes, gene_rankings, lineage_col)
    table = lineage_preference(observations, threshold=ESSENTIALITY_THRESHOLD)
    if table.empty:
        return result
    result["preference"] = table
    result["metadata"]["lineage_preference"] = table.to_dict(orient="records")
    result["metadata"]["lineage_preference_definition"] = (
        "Per cancer type, mean and median score minus that gene's all-lineage average, plus fraction below threshold. "
        "Negative means a lower score in that type. Descriptive only; not a test of selectivity."
    )
    return result


def _resolve_reference_symbol(raw: str, rankings: pd.DataFrame | None):
    token = str(raw or "").strip().upper()
    if not token or rankings is None or rankings.empty:
        return token, None, False
    lookup = set(rankings["gene_upper"])
    if token in lookup:
        return token, None, False
    alias_target = GENE_ALIASES.get(token)
    if alias_target and alias_target in lookup:
        return alias_target, f"{token} → {alias_target}", False
    return token, None, True


def figure_controls(rankings=None):
    """Chart-only settings. Kept beside the plot so changing them does not rebuild the page."""
    st.session_state.setdefault("reference_a", "MYC")
    st.session_state.setdefault("reference_b", "PTEN")
    st.session_state.setdefault("show_labels", True)
    st.session_state.setdefault("point_size", 4)
    with st.expander(ui("Figure settings", "图形设置"), expanded=False, icon=":material/tune:"):
        c1, c2 = st.columns(2)
        with c1:
            raw_a = st.text_input(t("essential"), key="reference_a").strip().upper()
        with c2:
            raw_b = st.text_input(t("nonessential"), key="reference_b").strip().upper()
        essential, alias_a, missing_a = _resolve_reference_symbol(raw_a, rankings)
        nonessential, alias_b, missing_b = _resolve_reference_symbol(raw_b, rankings)
        notes = []
        if alias_a:
            notes.append(alias_a)
        if alias_b:
            notes.append(alias_b)
        if notes:
            st.caption(ui("Reference alias resolved: ", "参考基因别名解析：") + ", ".join(notes))
        missing_refs = [ref for ref, miss in ((raw_a, missing_a), (raw_b, missing_b)) if miss]
        if missing_refs:
            st.caption(ui("Reference gene not found in current scope: ", "当前范围内未找到参考基因：") + ", ".join(missing_refs))
        o1, o2, o3 = st.columns([1.2, 1.4, 1.2], vertical_alignment="center")
        with o1:
            labels = st.checkbox(t("show_labels"), key="show_labels")
        with o2:
            size = st.slider(t("point_size"), 2, 8, key="point_size")
        with o3:
            layout_mode = st.selectbox(
                ui("Canvas", "画布比例"),
                ["compact", "wide"],
                key="rank_layout",
                format_func={"compact": ui("Compact 4:3", "紧凑 4:3"), "wide": ui("Wide 16:10", "宽屏 16:10")}.get,
            )
    width, height = (840, 500) if layout_mode == "wide" else (640, 480)
    return essential, nonessential, labels, size, width, height


@st.fragment
def render_curve_result(prefix):
    """Redraw labels and highlight markers without rerunning dataset loading."""
    result = st.session_state.get(f"{prefix}_result")
    if not result:
        return
    essential, nonessential, labels, size, width, height = figure_controls(result["rankings"])
    n_cell_lines = result["metadata"].get("cohort_rows", 0)
    if prefix == "rank":
        fig = create_rank_plot(result["rankings"], result["genes"], essential, nonessential,
                               n_cell_lines, labels, size, width=width, height=height)
    else:
        cfg = result["config"]
        fig = create_multilayer_rank_plot(
            result["rankings"], cfg["background"], cfg["highlight"],
            cfg["background_color"], cfg["highlight_color"],
            essential, nonessential, n_cell_lines, labels, size, width=width, height=height)
    centered_plot(fig)
    result_downloads(result, fig, prefix)
    render_deep_dive(result, prefix)


def input_draft(text='', file=None, **options):
    return {'text': text, 'file_name': file.name if file else None,
            'file_sha256': hashlib.sha256(file.getvalue()).hexdigest() if file else None,
            **options}


def input_signature(draft):
    return hashlib.sha256(json.dumps(draft, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def apply_input_preset(prefix, clear=False):
    if prefix == 'codep':
        mode = st.session_state.get('codep_mode', 'scan')
        if mode == 'scan':
            st.session_state['codep_query'] = '' if clear else 'MDM2'
        else:
            st.session_state['codep_genes'] = '' if clear else 'TP53, CHEK2, MDM2, MDM4, PPM1D, USP7'
        return
    presets = {
        'rank': {'genes': 'KRAS, SOX10, BRAF, ESR1, TP53'},
        'box': {'genes': 'SOX10, ESR1'},
        'multi': {'bg': 'CDK1, CDK2, CDK4, CDK6, CCNB1, CCND1, CCNE1', 'hl': 'PLK1, AURKA'},
    }
    for key, value in presets[prefix].items():
        st.session_state[f'{prefix}_{key}'] = '' if clear else value
    if prefix in {'rank', 'box'}:
        st.session_state[f'{prefix}_upload_mode'] = False
    if prefix == 'box':
        st.session_state.pop('box_candidates', None)


def apply_curated_preset(prefix, genes_text):
    st.session_state[f'{prefix}_genes'] = genes_text
    if prefix in {'rank', 'box'}:
        st.session_state[f'{prefix}_upload_mode'] = False
    if prefix == 'box':
        st.session_state.pop('box_candidates', None)
    st.session_state['_autorun_view'] = prefix


def apply_multi_preset(bg_text, hl_text):
    st.session_state['multi_bg'] = bg_text
    st.session_state['multi_hl'] = hl_text
    st.session_state['_autorun_view'] = 'multi'


def apply_codep_preset(mode, text):
    st.session_state['codep_mode'] = mode
    st.session_state['codep_query' if mode == 'scan' else 'codep_genes'] = text
    st.session_state['_autorun_view'] = 'codep'


def send_genes_to_view(target_view, genes):
    chosen = list(dict.fromkeys(genes))
    if target_view == 'box':
        st.session_state['box_genes'] = ', '.join(chosen[:8])
        st.session_state['box_upload_mode'] = False
        st.session_state.pop('box_candidates', None)
    elif target_view == 'codep':
        if len(chosen) == 1:
            st.session_state['codep_mode'] = 'scan'
            st.session_state['codep_query'] = chosen[0]
        else:
            st.session_state['codep_mode'] = 'matrix'
            st.session_state['codep_genes'] = ', '.join(chosen[:25])
    else:
        st.session_state['rank_genes'] = ', '.join(chosen)
        st.session_state['rank_upload_mode'] = False
    st.session_state['analysis_tab'] = target_view
    st.session_state['_autorun_view'] = target_view


def sync_share_url(view_name, genes, active_cohort, extra=None):
    params = {'view': view_name, 'genes': ','.join(genes[:30])}
    if active_cohort:
        params['cohort'] = ','.join(active_cohort)
    if extra:
        for key, value in extra.items():
            if value:
                params[key] = str(value)
    st.query_params.clear()
    for key, value in params.items():
        st.query_params[key] = value
    st.session_state['_share_copied_view'] = view_name


def render_primary_action(prefix, is_dirty=False):
    """Single dominant Apple-blue primary CTA that quietly morphs when edits are pending."""
    if is_dirty:
        label = ui('Update chart', '更新图表')
        icon = ':material/refresh:'
    else:
        label = ui('Run analysis', '生成图表')
        icon = None
    return st.button(label, type='primary', key=f'{prefix}_run', icon=icon, width='stretch')


def _reset_view_substate(prefix):
    """Clear per-result deep-dive & export-dimension widget state on a fresh Run."""
    for suffix in ('_pref_style', '_pref_genes', '_top_gene', '_corr_genes'):
        st.session_state.pop(f'{prefix}{suffix}', None)
    if prefix == 'box':
        st.session_state.pop('box_height', None)
    if prefix == 'codep':
        for key in ('codep_scan_partner', 'codep_pair_x', 'codep_pair_y', 'codep_scan_topn'):
            st.session_state.pop(key, None)
    for key in [k for k in st.session_state if k.startswith(f'_export_{prefix}_')]:
        del st.session_state[key]


def render_examples_popover(prefix, mode=None):
    """Single sleek Examples menu replacing redundant Load example + inline preset rows."""
    available = set(gene_rankings['gene_upper'])
    if prefix in {'rank', 'box'}:
        lead = (
            (
                'rank_preset_iconic',
                ui('Iconic cancer panel', '经典癌症依赖谱'),
                'KRAS, SOX10, BRAF, ESR1, TP53',
                {'KRAS', 'SOX10', 'BRAF', 'ESR1', 'TP53'},
                'KRAS · SOX10 · BRAF · ESR1 · TP53',
            )
            if prefix == 'rank'
            else (
                'box_preset_contrast',
                ui('Melanoma vs Breast', '黑色素瘤 vs 乳腺癌'),
                'SOX10, ESR1',
                {'SOX10', 'ESR1'},
                'SOX10 · ESR1',
            )
        )
        items = [
            lead,
            (
                f'{prefix}_preset_lineage',
                ui('Lineage drivers', '谱系驱动基因'),
                'SOX10, IRF4, ESR1, FOXA1',
                {'SOX10', 'IRF4', 'ESR1', 'FOXA1'},
                'SOX10 · IRF4 · ESR1 · FOXA1',
            ),
            (
                f'{prefix}_preset_oncogene',
                ui('Oncogene addiction', '癌基因成瘾'),
                'KRAS, BRAF, EGFR, ERBB2',
                {'KRAS', 'BRAF', 'EGFR', 'ERBB2'},
                'KRAS · BRAF · EGFR · ERBB2',
            ),
            (
                f'{prefix}_preset_repair',
                ui('DNA repair', 'DNA 修复通路'),
                'BRCA1, BRCA2, ATR, PARG',
                {'BRCA1', 'BRCA2', 'ATR', 'PARG'},
                'BRCA1 · BRCA2 · ATR · PARG',
            ),
        ]
        visible = [item for item in items if item[3] <= available]
        if not visible:
            return
        with st.popover(ui('Examples', '示例'), icon=':material/auto_awesome:'):
            st.caption(ui('Select an example to load and plot immediately.', '点击任一示例即可载入并立即绘图。'))
            for key, label, genes_text, _, tip in visible:
                st.button(label, key=key, width='stretch', help=tip,
                          on_click=apply_curated_preset, args=(prefix, genes_text))
        return

    if prefix == 'multi':
        items = [
            (
                'multi_preset_cycle',
                ui('Cell cycle · Mitotic', '细胞周期 · 有丝分裂'),
                'CDK1, CDK2, CDK4, CDK6, CCNB1, CCND1, CCNE1',
                'PLK1, AURKA',
                {'CDK1', 'CDK2', 'CDK4', 'CDK6', 'CCNB1', 'CCND1', 'CCNE1', 'PLK1', 'AURKA'},
                'Background: CDKs & Cyclins · Highlight: PLK1, AURKA',
            ),
            (
                'multi_preset_mapk',
                ui('RTK / MAPK · BRAF axis', 'RTK / MAPK · BRAF 轴'),
                'EGFR, ERBB2, KRAS, NRAS, RAF1, MAP2K1, MAPK1',
                'BRAF, DUSP4',
                {'EGFR', 'ERBB2', 'KRAS', 'NRAS', 'RAF1', 'MAP2K1', 'MAPK1', 'BRAF', 'DUSP4'},
                'Background: RTK-RAS-MAPK · Highlight: BRAF, DUSP4',
            ),
        ]
        visible = [item for item in items if item[4] <= available]
        if not visible:
            return
        with st.popover(ui('Examples', '示例'), icon=':material/auto_awesome:'):
            st.caption(ui('Select a pathway preset to load and plot immediately.', '点击任一通路示例即可载入并立即绘图。'))
            for key, label, bg_t, hl_t, _, tip in visible:
                st.button(label, key=key, width='stretch', help=tip,
                          on_click=apply_multi_preset, args=(bg_t, hl_t))
        return

    if prefix == 'codep':
        if mode == 'scan':
            items = [
                ('codep_preset_mdm2', ui('MDM2 · p53 axis', 'MDM2 · p53 拮抗轴'), 'MDM2', {'MDM2'},
                 ui('Symmetric positive (PPM1D/MDM4/USP7) & negative (TP53/CDKN1A) partners', '正相关（PPM1D/MDM4）与负相关（TP53/CDKN1A）高度对称')),
                ('codep_preset_braf', ui('BRAF · MAPK', 'BRAF · MAPK 通路'), 'BRAF', {'BRAF'},
                 ui('MAPK signaling & melanoma lineage partners (MAPK1, SOX10, MITF)', 'MAPK 信号与黑色素瘤谱系共依赖（MAPK1、SOX10、MITF）')),
                ('codep_preset_brca1', ui('BRCA1 · HR repair', 'BRCA1 · 同源重组'), 'BRCA1', {'BRCA1'},
                 ui('Homologous recombination & DNA damage response partners', '同源重组与 DNA 损伤修复共依赖网络')),
            ]
        else:
            items = [
                (
                    'codep_preset_p53',
                    ui('p53 · MDM2 axis', 'p53 · MDM2 调控轴'),
                    'TP53, CHEK2, MDM2, MDM4, PPM1D, USP7',
                    {'TP53', 'CHEK2', 'MDM2', 'MDM4', 'PPM1D', 'USP7'},
                    'TP53 · CHEK2 vs MDM2 · MDM4 · PPM1D · USP7',
                ),
                (
                    'codep_preset_mapk',
                    ui('MAPK · Melanoma', 'MAPK · 黑色素瘤网络'),
                    'BRAF, MAPK1, MAP2K1, SOX10, MITF',
                    {'BRAF', 'MAPK1', 'MAP2K1', 'SOX10', 'MITF'},
                    'BRAF · MAPK1 · MAP2K1 · SOX10 · MITF',
                ),
                (
                    'codep_preset_repair_m',
                    ui('DNA repair', 'DNA 修复复合体'),
                    'BRCA1, BRCA2, ATR, PARG',
                    {'BRCA1', 'BRCA2', 'ATR', 'PARG'},
                    'BRCA1 · BRCA2 · ATR · PARG',
                ),
            ]
        visible = [item for item in items if item[3] <= available]
        if not visible:
            return
        with st.popover(ui('Examples', '示例'), icon=':material/auto_awesome:'):
            st.caption(ui('Select an example to load and run immediately.', '点击任一示例即可载入并立即计算。'))
            for key, label, text, _, tip in visible:
                st.button(label, key=key, width='stretch', help=tip,
                          on_click=apply_codep_preset, args=(mode, text))


def render_run_status(prefix, draft, result):
    if result is None:
        st.markdown(
            f'<div class="ready-note">{escape(ui("Enter gene symbols or pick from Examples above, then run.", "输入基因符号或从上方“示例”中选择，点击生成图表。"))}</div>',
            unsafe_allow_html=True,
        )
        return
    is_dirty = result.get('input_signature') != input_signature(draft)
    meta = result['metadata']
    scope = ', '.join(meta['cohort_lineages'])
    if scope == 'ALL':
        scope = ui('All cancer types', '全部癌种')
    timestamp = datetime.fromisoformat(meta['created_at_utc']).strftime('%H:%M UTC')
    share_extra = None
    if prefix == 'codep':
        metric_label = 'Pearson r' if result['method'] == 'pearson' else 'Spearman ρ'
        if result['method'] != 'pearson':
            share_extra = {'method': result['method']}
        if result['mode'] == 'scan':
            summary = ui(
                f"{result['query_gene']} · {len(result['scan_df']):,} genes · {metric_label} · {meta['cohort_rows']:,} lines · {scope}",
                f"{result['query_gene']} · {len(result['scan_df']):,} 基因 · {metric_label} · {meta['cohort_rows']:,} 细胞系 · {scope}",
            )
            transfer_genes = [result['query_gene']] + result['scan_df']['gene'].head(5).tolist()
            share_genes = [result['query_gene']]
        else:
            summary = ui(
                f"{len(result['genes'])} genes · {metric_label} · {meta['cohort_rows']:,} lines · {scope}",
                f"{len(result['genes'])} 基因 · {metric_label} · {meta['cohort_rows']:,} 细胞系 · {scope}",
            )
            transfer_genes = result['genes']
            share_genes = result['genes']
        actions = [
            (ui('Ranking', '基因排名'), 'codep_send_rank', 'rank', transfer_genes),
        ]
        if lineage_col:
            actions.append((ui('Lineages', '癌种比较'), 'codep_send_box', 'box', transfer_genes))
    else:
        summary = ui(f"{len(result['genes'])} genes · {meta['cohort_rows']:,} cell lines · {scope}",
                     f"{len(result['genes'])} 个基因 · {meta['cohort_rows']:,} 个细胞系 · {scope}")
        share_genes = result['genes']
        actions = []
        if prefix in {'rank', 'multi'}:
            target_genes = result['config'].get('highlight') or result['genes'] if prefix == 'multi' else result['genes']
            if lineage_col:
                actions.append((ui('Lineages', '癌种比较'), f'{prefix}_send_box', 'box', target_genes))
            actions.append((ui('Co-dependency', '共依赖'), f'{prefix}_send_codep', 'codep', target_genes))
        elif prefix == 'box':
            actions.append((ui('Ranking', '基因排名'), 'box_send_rank', 'rank', result['genes']))
            actions.append((ui('Co-dependency', '共依赖'), 'box_send_codep', 'codep', result['genes']))

    dot_class = 'status-dot pending' if is_dirty else 'status-dot'
    pending_html = (
        f'<span class="status-pending">· {escape(ui("Unsaved edit", "待更新"))}</span>'
        if is_dirty else ''
    )
    col_specs = [2.65] + [0.95] * len(actions) + [0.85]
    cols = st.columns(col_specs, vertical_alignment='center')
    with cols[0]:
        st.markdown(
            f'<div class="status-strip"><span class="{dot_class}"></span>'
            f'<span><strong>{escape(summary)}</strong></span>{pending_html}'
            f'<span class="status-time">{escape(timestamp)}</span></div>',
            unsafe_allow_html=True,
        )
    for col, (label, btn_key, target_view, target_genes) in zip(cols[1:-1], actions):
        with col:
            st.button(
                label, key=btn_key, icon=':material/arrow_forward:', width='stretch',
                on_click=send_genes_to_view, args=(target_view, target_genes),
            )
    with cols[-1]:
        active_cohort = [x for x in meta.get('cohort_lineages', []) if x != 'ALL']
        st.button(
            ui('Share', '分享'),
            key=f'{prefix}_share_url', icon=':material/share:', width='stretch',
            on_click=sync_share_url, args=(prefix, share_genes, active_cohort, share_extra),
            help=ui('Write current view, genes, and cohort into the browser URL for sharing.',
                    '将当前视图、基因列表和癌种范围写入浏览器地址栏以便分享。'),
        )
    if st.session_state.pop('_share_copied_view', None) == prefix:
        st.caption(ui('URL parameters updated — copy the address bar to share this view.',
                      '已将当前分析参数写入地址栏，直接复制浏览器链接即可分享。'))


def update_preferences(preference):
    options = TRANSLATIONS if preference == 'lang' else THEMES
    value = st.session_state.get(f'{preference}_select')
    if value in options:
        st.session_state[preference] = value


def render_preferences():
    """Always-visible compact language and theme controls at the top of the sidebar."""
    c1, c2 = st.sidebar.columns(2)
    with c1:
        st.selectbox('Language / 语言', ['en', 'zh'], key='lang_select',
                     index=['en', 'zh'].index(st.session_state.lang),
                     format_func={'en': 'English', 'zh': '中文'}.get,
                     on_change=update_preferences, args=('lang',))
    with c2:
        st.selectbox(t('theme'), ['light', 'dark'], key='theme_select',
                     index=['light', 'dark'].index(st.session_state.theme),
                     format_func={x: t(x) for x in ['light', 'dark']}.get,
                     on_change=update_preferences, args=('theme',))


# Apply URL query parameters once per session if explicitly provided.
if not st.session_state.get('_qp_initialized'):
    st.session_state['_qp_initialized'] = True
    qp = st.query_params
    qp_lang = qp.get('lang')
    if qp_lang in TRANSLATIONS:
        st.session_state.lang = qp_lang
        st.session_state.lang_select = qp_lang
    qp_theme = qp.get('theme')
    if qp_theme in THEMES:
        st.session_state.theme = qp_theme
        st.session_state.theme_select = qp_theme
    qp_view = qp.get('view')
    if qp_view in {'rank', 'box', 'multi', 'codep'}:
        st.session_state['analysis_tab'] = qp_view
    qp_method = qp.get('method')
    if qp_method in {'pearson', 'spearman'}:
        st.session_state['codep_method'] = qp_method
    qp_genes = qp.get('genes')
    if qp_genes:
        parsed_qp = parse_gene_text(qp_genes)['genes']
        if parsed_qp:
            if qp_view == 'codep':
                if len(parsed_qp) == 1:
                    st.session_state['codep_mode'] = 'scan'
                    st.session_state['codep_query'] = parsed_qp[0]
                else:
                    st.session_state['codep_mode'] = 'matrix'
                    st.session_state['codep_genes'] = ', '.join(parsed_qp[:25])
                st.session_state['_autorun_view'] = 'codep'
            else:
                target_view = qp_view if qp_view in {'rank', 'box'} else 'rank'
                st.session_state[f'{target_view}_genes'] = ', '.join(parsed_qp)
                st.session_state[f'{target_view}_upload_mode'] = False
                st.session_state['_autorun_view'] = target_view
    qp_cohort = qp.get('cohort')
    if qp_cohort:
        st.session_state['_qp_cohort'] = [item.strip() for item in str(qp_cohort).split(',') if item.strip()]

for draft_key in ['rank_genes', 'box_genes', 'multi_bg', 'multi_hl', 'rank_upload_mode', 'box_upload_mode',
                  'multi_bg_color', 'multi_hl_color', 'box_lineages', 'box_order', 'box_points',
                  'box_cols', 'box_compare', 'reference_a', 'reference_b', 'show_labels',
                  'point_size', 'rank_layout', 'codep_mode', 'codep_method', 'codep_query', 'codep_genes']:
    if draft_key in st.session_state:
        st.session_state[draft_key] = st.session_state[draft_key]

for preference, options in [('lang', TRANSLATIONS), ('theme', THEMES)]:
    if st.session_state.get(f'{preference}_select') not in options:
        st.session_state[f'{preference}_select'] = st.session_state[preference]
inject_css()
render_preferences()
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
        if key in {'rank_result', 'box_result', 'multi_result', 'codep_result', 'codep_scan_partner',
                   'codep_pair_x', 'codep_pair_y', 'codep_scan_topn', 'box_candidates', 'box_display_genes',
                   'gene_columns', 'custom_mapping', 'box_lineages', 'box_height', 'disease_filter'} or key.startswith('_export_'):
            del st.session_state[key]
    st.session_state['_source_id'] = source_id
if source_mode == 'custom' and uploaded_file is None:
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
        st.error(ui('No data rows. Upload a nonempty score matrix.', '数据没有有效行，请上传非空分数矩阵。'))
        st.stop()
except (ValueError, UnicodeError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
    st.error(ui('Cannot read this CSV: ', '无法读取 CSV：') + str(exc))
    st.stop()

df = crispr_data
source_sha = source_sha or df.attrs.get('sha256')
provenance = df.attrs.get('provenance') or {}
lineage_col = find_lineage_column(df)
disease_col = find_disease_column(df)
gene_columns = None
selected_diseases = []
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
        qp_cohort = st.session_state.pop('_qp_cohort', None)
        if qp_cohort:
            valid_qp = [x for x in qp_cohort if x in all_lineages]
            if valid_qp:
                st.session_state.cohort = valid_qp
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
    if disease_col:
        disease_labels = df_scope[disease_col].fillna('Unknown').astype(str).str.strip().replace('', 'Unknown')
        all_diseases = sorted(disease_labels.unique().tolist())
        if any(x not in all_diseases for x in st.session_state.get('disease_filter', [])):
            st.session_state.disease_filter = []
        selected_diseases = st.multiselect(
            ui('Disease / subtype (optional)', '细分病种 / 亚型（可选）'),
            all_diseases, key='disease_filter',
            placeholder=ui('All subtypes in scope', '当前范围内全部亚型'),
        )
        if selected_diseases:
            df_scope = df_scope.loc[disease_labels.isin(selected_diseases)]

context = json.dumps({'sha256': source_sha, 'cohort': sorted(cohort), 'diseases': sorted(selected_diseases),
                      'gene_columns': gene_columns, 'source': uploaded_file.name if uploaded_file else matrix_path.name}, sort_keys=True)
context_hash = hashlib.sha256(context.encode()).hexdigest()
if st.session_state.get('_analysis_context') != context_hash:
    for key in list(st.session_state):
        if key in {'rank_result', 'box_result', 'multi_result', 'codep_result', 'codep_scan_partner',
                   'codep_pair_x', 'codep_pair_y', 'codep_scan_topn', 'box_candidates', 'box_display_genes',
                   'box_lineages', 'box_height'} or key.startswith('_export_'):
            del st.session_state[key]
    st.session_state['_analysis_context'] = context_hash

gene_rankings, diagnostics = compute_gene_rankings(context_hash, df_scope, gene_columns)
all_gene_options = sorted(gene_rankings['gene'].tolist()) if not gene_rankings.empty else []
n_cell_lines = len(df_scope)
hero_version = t('custom_dataset') if uploaded_file is not None else provenance.get('release', DATA_VERSION)
scope_name = ', '.join(cohort) if cohort else ui('All cancer types', '全部癌种')
if selected_diseases:
    scope_name += ' · ' + ', '.join(selected_diseases)
stats_pills = (
    f'<span class="hero-metric"><strong>{n_cell_lines:,}</strong> {t("cell_lines")}</span>'
    f'<span class="hero-sep">/</span>'
    f'<span class="hero-metric"><strong>{len(gene_rankings):,}</strong> {t("gene_count")}</span>'
    f'<span class="hero-version">{escape(hero_version)}</span>'
) if not gene_rankings.empty else f'<span class="hero-version">{escape(hero_version)}</span>'
st.markdown(f'<section class="hero-shell"><div class="hero-copy">'
            f'<div class="hero-kicker">GENE DEPENDENCY EXPLORER</div>'
            f'<h1 class="main-header">CRISPR Score Analyzer</h1>'
            f'<p class="sub-header">{t("app_subtitle")}</p></div>'
            f'<div class="hero-stats">{stats_pills}</div></section>', unsafe_allow_html=True)
if gene_rankings.empty:
    st.warning(ui('No analyzable gene columns. For custom CSVs, select score columns in the sidebar. Gene symbols must be unique, with at least one finite numeric score per gene.',
                  '没有可分析的基因列。自定义 CSV 请在侧栏选择分数列；基因名不能重复，且每个基因至少需要一个有效数值分数。'))
    st.download_button(ui('Download validation report', '下载数据校验报告'),
                       lambda frame=diagnostics: frame.to_csv(index=False).encode('utf-8-sig'),
                       'column_diagnostics.csv', 'text/csv', key='validation_report', on_click='ignore')
    st.stop()

if cohort or selected_diseases:
    st.markdown(f'<div class="scope-note">{ui("Active scope", "当前分析范围")} · {escape(scope_name)}</div>',
                unsafe_allow_html=True)
with st.sidebar.expander(ui('Data overview', '数据概况')):
    st.metric(t('essential_genes').format(threshold=ESSENTIALITY_THRESHOLD),
              f"{(gene_rankings.mean_score < ESSENTIALITY_THRESHOLD).sum():,}")
    st.caption(ui('Mean-score range: ', '平均分范围：') +
               f'{gene_rankings.mean_score.min():.2f} – {gene_rankings.mean_score.max():.2f}')
    st.caption(ui('Lower score means stronger dependency. Rankings use means within the analysis scope; they do not establish cancer-type selectivity. The −0.5 line is a screening reference.',
                  '分数越低，依赖越强。排名使用分析范围内的平均分，本身不能证明癌种特异性。−0.5 是筛选参考线。'))
source_files = (provenance.get('files') or {}) if uploaded_file is None else {}
dataset_metadata = {
    'data_source': uploaded_file.name if uploaded_file else 'bundled matrix',
    'data_file': uploaded_file.name if uploaded_file else matrix_path.name,
    'data_release': 'custom' if uploaded_file else provenance.get('release', DATA_VERSION),
    'data_revision': None if uploaded_file else (source_files.get('CRISPRGeneEffect.csv') or {}).get('sha256'),
    'data_sha256': source_sha, 'source_files': source_files or None,
    'cohort_lineages': cohort or ['ALL'],
    'cohort_diseases': selected_diseases or None,
    'lineage_column': lineage_col, 'gene_column_mapping': gene_columns,
    'cohort_rows': n_cell_lines, 'ranked_genes': len(gene_rankings),
    'mean_score_cutoff': ESSENTIALITY_THRESHOLD,
    'tool_version': TOOL_VERSION,
    'tool_released': TOOL_RELEASED,
    'zenodo_version': ZENODO_VERSION,
    'zenodo_version_doi': ZENODO_VERSION_DOI,
    'zenodo_concept_doi': ZENODO_CONCEPT_DOI,
    'software': get_software_versions(),
    'analysis_schema': 4,
}

st.session_state.setdefault('analysis_tab', 'rank')
with st.container(border=False, key='analysis_nav'):
    view = st.segmented_control(ui('Analysis', '分析'), ['rank', 'box', 'multi', 'codep'],
                                format_func={'rank': t('tab1'), 'box': t('tab2'), 'multi': t('tab3'), 'codep': t('tab4')}.get,
                                key='analysis_tab', selection_mode='single', label_visibility='collapsed', width='stretch') or 'rank'
autorun_view = st.session_state.pop('_autorun_view', None)
allow_initial_autorun = not st.session_state.get('_initial_autorun_done')
if allow_initial_autorun:
    st.session_state['_initial_autorun_done'] = True

flash_info = st.session_state.get('_flash_info')
if flash_info and flash_info[0] == view:
    st.session_state.pop('_flash_info', None)
    st.info(flash_info[1])

if view == 'rank':
    st.session_state.setdefault('rank_upload_mode', False)
    with st.container(border=False, key='rank_input'):
        with st.container(horizontal=True, horizontal_alignment='distribute', vertical_alignment='center'):
            with st.container(horizontal=True, vertical_alignment='center', width='content'):
                st.markdown(f'<span class="input-title">{t("gene_list")}</span>', unsafe_allow_html=True)
                upload_mode = st.toggle(ui('Upload file', '上传文件'), key='rank_upload_mode')
            method = 'file' if upload_mode else 'text'
            with st.container(horizontal=True, vertical_alignment='center', width='content'):
                render_examples_popover('rank')
                if method == 'text':
                    render_gene_quick_add('rank', 'rank_genes')
                st.button(ui('Clear', '清空'), key='rank_clear', on_click=apply_input_preset, args=('rank', True))
        gene_text, gene_file = read_gene_input('rank', 'KRAS, SOX10, BRAF, ESR1, TP53', method)
        draft = input_draft(gene_text, gene_file, method=method)
        prev_result = st.session_state.get('rank_result')
        is_dirty = prev_result is not None and prev_result.get('input_signature') != input_signature(draft)
        initial_autorun = False
        if (allow_initial_autorun and uploaded_file is None and not cohort and not selected_diseases
                and prev_result is None and method == 'text' and gene_text.strip()):
            _, init_matched, _, _ = classify_input(gene_text, None)
            if init_matched:
                initial_autorun = True
        submitted = render_primary_action('rank', is_dirty) or (autorun_view == 'rank') or initial_autorun
    if submitted:
        try:
            parsed, matched, missing, resolved = classify_input(gene_text, gene_file)
            if matched:
                _reset_view_substate('rank')
                st.session_state.rank_result = attach_preference(build_result(
                    matched, {'targets': {'parsed': parsed, 'matched': matched, 'missing': missing, 'resolved': resolved}},
                    draft=draft), matched)
                if is_dirty:
                    st.rerun()
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
    if not lineage_col:
        st.info(t('lineage_missing'))
    else:
        st.session_state.setdefault('box_upload_mode', False)
        with st.container(border=False, key='box_input'):
            with st.container(horizontal=True, horizontal_alignment='distribute', vertical_alignment='center'):
                with st.container(horizontal=True, vertical_alignment='center', width='content'):
                    st.markdown(f'<span class="input-title">{t("gene_list")}</span>', unsafe_allow_html=True)
                    upload_mode = st.toggle(ui('Upload file', '上传文件'), key='box_upload_mode')
                method = 'file' if upload_mode else 'text'
                with st.container(horizontal=True, vertical_alignment='center', width='content'):
                    render_examples_popover('box')
                    if method == 'text':
                        render_gene_quick_add('box', 'box_genes')
                    st.button(ui('Clear', '清空'), key='box_clear', on_click=apply_input_preset, args=('box', True))
            gene_text, gene_file = read_gene_input(
                'box', 'SOX10, ESR1', method,
                placeholder=ui('Enter up to 8 genes to compare across cancer types (e.g. SOX10, ESR1)',
                               '输入最多 8 个基因比较各癌种依赖分布（如 SOX10、ESR1）'),
            )
            with st.expander(ui('Comparison options', '比较选项'), icon=':material/tune:'):
                box_lineages = st.multiselect(ui('Cancer types shown in the chart', '图中展示的癌种'),
                                              cohort or all_lineages, key='box_lineages',
                                              placeholder=ui('All types in the analysis scope', '分析范围内的全部癌种'),
                                              help=ui('This only limits the chart. Gene ranks still use the analysis scope in the sidebar.',
                                                      '此处仅限制图中分组；基因排名仍按侧栏的分析范围计算。'))
                opt1, opt2, opt3 = st.columns(3)
                with opt1:
                    sort_mode = st.selectbox(ui('Order', '排序'), ['median', 'alphabetical'], key='box_order',
                                             format_func={'median': ui('Median score', '中位数'), 'alphabetical': ui('Alphabetical', '字母顺序')}.get)
                with opt2:
                    box_compare = st.selectbox(ui('View mode', '比较视角'), ['by_lineage', 'by_gene'], key='box_compare',
                                               format_func={'by_lineage': ui('By cancer type', '按癌种拆分'),
                                                            'by_gene': ui('Compare genes together', '同轴比较基因')}.get)
                with opt3:
                    box_cols = st.selectbox(ui('Columns', '排布列数'), [1, 2], key='box_cols',
                                            format_func={1: ui('1 column', '单列堆叠'), 2: ui('2 columns', '双列并排')}.get)
                all_points = st.checkbox(ui('Show all cell lines', '显示全部细胞系散点'), key='box_points')
            draft = input_draft(gene_text, gene_file, method=method, lineages=box_lineages,
                                sort=sort_mode, all_points=all_points, columns=box_cols, compare_mode=box_compare)
            prev_result = st.session_state.get('box_result')
            is_dirty = prev_result is not None and prev_result.get('input_signature') != input_signature(draft)
            submitted = render_primary_action('box', is_dirty) or (autorun_view == 'box')
        candidate = st.session_state.get('box_candidates')
        if candidate and candidate['input_signature'] != input_signature(draft):
            st.session_state.pop('box_candidates', None)
        if submitted:
            try:
                parsed, matched, missing, resolved = classify_input(gene_text, gene_file)
                match = {'targets': {'parsed': parsed, 'matched': matched, 'missing': missing, 'resolved': resolved}}
                config = {'lineages': box_lineages, 'sort': sort_mode, 'all_points': all_points,
                          'columns': box_cols, 'compare_mode': box_compare}
                st.session_state.pop('box_candidates', None)
                if len(matched) > 8:
                    st.session_state.box_candidates = {'genes': matched, 'matches': match, 'config': config,
                                                       'draft': copy.deepcopy(draft), 'input_signature': input_signature(draft)}
                    st.session_state.pop('box_display_genes', None)
                elif matched:
                    frame = get_lineage_frame(df_scope, matched, gene_rankings, lineage_col)
                    if box_lineages:
                        frame = frame[frame.lineage.isin(box_lineages)]
                    _reset_view_substate('box')
                    st.session_state.box_result = attach_preference(build_result(matched, match, config, frame, draft), matched)
                    if is_dirty:
                        st.rerun()
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
                _reset_view_substate('box')
                st.session_state.box_result = attach_preference(
                    build_result(chosen, candidate['matches'], candidate['config'], frame, candidate['draft']), chosen)
                del st.session_state.box_candidates
                st.rerun()
        result = st.session_state.get('box_result')
        render_run_status('box', draft, result)
        render_draft_feedback(gene_text, gene_file, draft, result)
        if result:
            match = result['matches']['targets']
            show_matches(match['parsed'], match['matched'], match['missing'], resolved=match.get('resolved'))
            if len(result['genes']) < len(match['matched']):
                st.caption(ui('Displayed subset: ', '当前展示子集：') + ', '.join(result['genes']))
            valid_data = result['plot_data'].dropna(subset=['crispr_score'])
            if valid_data.empty:
                st.warning(ui('No finite scores in the selected groups.', '所选分组没有有效分数。'))
            else:
                absent = [g for g in result['genes'] if g not in set(valid_data.gene)]
                if absent:
                    st.info(ui('No valid scores in these groups: ', '这些基因在所选分组没有有效分数：') + ', '.join(absent))
                cfg = result['config']
                fig = create_lineage_boxplot(
                    valid_data, result['genes'], cfg['sort'], cfg['all_points'],
                    columns=cfg.get('columns', 1), compare_mode=cfg.get('compare_mode', 'by_lineage'),
                )
                st.plotly_chart(fig, config=PLOT_CONFIG, width='stretch', theme=None)
                result_downloads(result, fig, 'box')
                render_deep_dive(result, 'box')

elif view == 'multi':
    st.session_state.setdefault('multi_bg', 'CDK1, CDK2, CDK4, CDK6, CCNB1, CCND1, CCNE1')
    st.session_state.setdefault('multi_hl', 'PLK1, AURKA')
    st.session_state.setdefault('multi_bg_color', '#8E99A2')
    st.session_state.setdefault('multi_hl_color', '#3C5488')
    with st.container(border=False, key='multi_input'):
        with st.container(horizontal=True, horizontal_alignment='distribute', vertical_alignment='center'):
            st.markdown(f'<span class="input-title">{ui("Gene-set layers", "基因集分层")}</span>', unsafe_allow_html=True)
            with st.container(horizontal=True, vertical_alignment='center', width='content'):
                render_examples_popover('multi')
                st.button(ui('Clear', '清空'), key='multi_clear', on_click=apply_input_preset, args=('multi', True))
        c1, c2 = st.columns(2)
        with c1:
            with st.container(horizontal=True, horizontal_alignment='distribute', vertical_alignment='center'):
                st.caption(t('bg_gene_set'))
                render_gene_quick_add('multi_bg', 'multi_bg', ui('Search', '搜索'))
            bg_text = st.text_area(t('bg_gene_set'), key='multi_bg', height=68, label_visibility='collapsed')
        with c2:
            with st.container(horizontal=True, horizontal_alignment='distribute', vertical_alignment='center'):
                st.caption(t('hl_gene_set'))
                render_gene_quick_add('multi_hl', 'multi_hl', ui('Search', '搜索'))
            hl_text = st.text_area(t('hl_gene_set'), key='multi_hl', height=68, label_visibility='collapsed')
        with st.expander(ui('Annotation colors', '标注颜色'), icon=':material/palette:'):
            color1, color2 = st.columns(2)
            with color1:
                bg_color = st.color_picker(t('bg_color'), key='multi_bg_color')
            with color2:
                hl_color = st.color_picker(t('hl_color'), key='multi_hl_color')
        draft = {'background': bg_text, 'highlight': hl_text, 'background_color': bg_color, 'highlight_color': hl_color}
        prev_result = st.session_state.get('multi_result')
        is_dirty = prev_result is not None and prev_result.get('input_signature') != input_signature(draft)
        submitted = render_primary_action('multi', is_dirty) or (autorun_view == 'multi')
    if submitted:
        bg_parsed, hl_parsed = parse_gene_text(bg_text), parse_gene_text(hl_text)
        bg_matched, bg_missing, bg_resolved = resolve_gene_list(gene_rankings, bg_parsed['genes'], GENE_ALIASES)
        hl_matched, hl_missing, hl_resolved = resolve_gene_list(gene_rankings, hl_parsed['genes'], GENE_ALIASES)
        matches = {'background': {'parsed': bg_parsed, 'matched': bg_matched, 'missing': bg_missing, 'resolved': bg_resolved},
                   'highlight': {'parsed': hl_parsed, 'matched': hl_matched, 'missing': hl_missing, 'resolved': hl_resolved}}
        if bg_matched or hl_matched:
            genes = list(dict.fromkeys(bg_matched + hl_matched))
            pref_genes = hl_matched if hl_matched else bg_matched[:PREFERENCE_MAX_LOLLIPOP_GENES]
            _reset_view_substate('multi')
            st.session_state.multi_result = attach_preference(build_result(
                genes, matches, {'background': bg_matched, 'highlight': hl_matched,
                                 'background_color': bg_color, 'highlight_color': hl_color}, draft=draft), pref_genes)
            if is_dirty:
                st.rerun()
        else:
            st.warning(ui('No matched genes. The previous result is unchanged.', '没有匹配的基因，保留上次运行结果。'))
    result = st.session_state.get('multi_result')
    render_run_status('multi', draft, result)
    if result is None or result.get('input_signature') != input_signature(draft):
        bg_p, bg_m, bg_miss, bg_res = classify_input(bg_text, None) if bg_text.strip() else ({}, [], [], [])
        hl_p, hl_m, hl_miss, hl_res = classify_input(hl_text, None) if hl_text.strip() else ({}, [], [], [])
        if _has_match_issues(bg_p, bg_miss, bg_res) or _has_match_issues(hl_p, hl_miss, hl_res):
            fb1, fb2 = st.columns(2)
            with fb1:
                show_matches(bg_p, bg_m, bg_miss, t('bg_gene_set'), bg_res)
            with fb2:
                show_matches(hl_p, hl_m, hl_miss, t('hl_gene_set'), hl_res)
    if result:
        bg_match = result['matches']['background']
        hl_match = result['matches']['highlight']
        if (_has_match_issues(bg_match['parsed'], bg_match['missing'], bg_match.get('resolved'))
                or _has_match_issues(hl_match['parsed'], hl_match['missing'], hl_match.get('resolved'))):
            m1, m2 = st.columns(2)
            with m1:
                show_matches(bg_match['parsed'], bg_match['matched'], bg_match['missing'],
                             t('bg_gene_set'), bg_match.get('resolved'))
            with m2:
                show_matches(hl_match['parsed'], hl_match['matched'], hl_match['missing'],
                             t('hl_gene_set'), hl_match.get('resolved'))
        render_curve_result('multi')

elif view == 'codep':
    st.session_state.setdefault('codep_mode', 'scan')
    st.session_state.setdefault('codep_method', 'pearson')
    st.session_state.setdefault('codep_query', 'MDM2')
    st.session_state.setdefault('codep_genes', 'TP53, CHEK2, MDM2, MDM4, PPM1D, USP7')
    with st.container(border=False, key='codep_input'):
        m_col, meth_col = st.columns([2.0, 1.0], vertical_alignment='bottom')
        with m_col:
            codep_mode = st.segmented_control(
                ui('Co-dependency mode', '共依赖分析模式'),
                ['scan', 'matrix'],
                format_func={
                    'scan': ui('Genome-wide scan (1 gene)', '单基因全基因组扫描'),
                    'matrix': ui('Multi-gene matrix (2–25 genes)', '多基因相关性矩阵'),
                }.get,
                key='codep_mode', selection_mode='single', width='stretch',
            ) or 'scan'
        with meth_col:
            codep_method = st.segmented_control(
                ui('Correlation metric', '相关性指标'),
                ['pearson', 'spearman'],
                format_func={'pearson': 'Pearson r', 'spearman': 'Spearman ρ'}.get,
                key='codep_method', selection_mode='single', width='stretch',
            ) or 'pearson'

        with st.container(horizontal=True, horizontal_alignment='distribute', vertical_alignment='center'):
            title_text = ui('Query gene', '查询基因') if codep_mode == 'scan' else t('gene_list')
            st.markdown(f'<span class="input-title">{title_text}</span>', unsafe_allow_html=True)
            with st.container(horizontal=True, vertical_alignment='center', width='content'):
                render_examples_popover('codep', mode=codep_mode)
                if codep_mode == 'scan':
                    render_gene_quick_add('codep', 'codep_query', single_select=True)
                else:
                    render_gene_quick_add('codep', 'codep_genes')
                st.button(ui('Clear', '清空'), key='codep_clear', on_click=apply_input_preset, args=('codep', True))

        if codep_mode == 'scan':
            if st.session_state.get('codep_query') is None:
                st.session_state['codep_query'] = ''
            selected_query = st.text_input(
                ui('Query gene (symbol or alias)', '查询基因（支持官方符号或别名）'),
                key='codep_query',
                autocomplete='off',
                label_visibility='collapsed',
                placeholder=ui(
                    'Enter a gene symbol or alias (e.g. MDM2, TEAD1, BRAF, HER2)…',
                    '输入基因符号或别名（如 MDM2、TEAD1、BRAF、HER2）…',
                ),
                help=ui(
                    'Correlates this gene against every other gene in the active sidebar scope '
                    '(including single-lineage cohorts). Accepts official symbols, HGNC aliases, '
                    'or DepMap headers like "TEAD1 (7003)". Use Search to browse matching symbols.',
                    '在当前侧栏癌种范围内，将该基因与其余所有基因逐一计算相关性（支持单癌种内共依赖扫描）；'
                    '兼容官方符号、HGNC 别名及 "TEAD1 (7003)" 格式，也可点击右上角“搜索”按前缀查找。',
                ),
            )
            active_text = str(selected_query or '')
            draft = input_draft(active_text, None, mode='scan', method=codep_method)
        else:
            active_text = st.text_area(
                t('gene_list'), height=68, key='codep_genes',
                label_visibility='collapsed',
                placeholder=ui(
                    'Enter 2–25 genes to compute a pairwise correlation matrix (e.g. TP53, CHEK2, MDM2, MDM4)',
                    '输入 2–25 个基因计算两两共依赖矩阵（如 TP53、CHEK2、MDM2、MDM4）',
                ),
                help=ui('Separate with newlines, spaces, tabs, commas or semicolons. DepMap headers like "MYC (4609)" are also accepted.',
                        '支持换行、空格、制表符、中英文逗号和分号，也兼容直接粘贴 "MYC (4609)" 格式。'),
            )
            draft = input_draft(active_text, None, mode='matrix', method=codep_method)

        prev_result = st.session_state.get('codep_result')
        is_dirty = prev_result is not None and prev_result.get('input_signature') != input_signature(draft)
        submitted = render_primary_action('codep', is_dirty) or (autorun_view == 'codep')

    if submitted:
        try:
            parsed = parse_gene_text(active_text)
            matched, missing, resolved = resolve_gene_list(gene_rankings, parsed['genes'], GENE_ALIASES)
            if codep_mode == 'scan':
                if not matched:
                    st.warning(ui('No matched genes. The previous result is unchanged.', '没有匹配的基因，保留上次运行结果。'))
                elif len(gene_rankings) < 2:
                    st.warning(ui('At least 2 ranked genes are needed for a co-dependency scan.', '至少需要 2 个基因才能进行共依赖扫描。'))
                else:
                    query_gene = matched[0]
                    info_note = None
                    if len(matched) > 1:
                        info_note = ui(
                            f'Scanning {query_gene} (first matched gene). Switch to Multi-gene matrix to correlate multiple genes together.',
                            f'正在扫描首个匹配基因 {query_gene}；如需同时比较多个基因的两两相关性，请切换到“多基因相关性矩阵”。',
                        )
                    with st.spinner(ui(
                        f'Scanning {len(gene_rankings) - 1:,} genes against {query_gene}…',
                        f'正在将 {query_gene} 与 {len(gene_rankings) - 1:,} 个基因计算全基因组共依赖…',
                    )):
                        scan_df = compute_codependency_scan(context_hash, df_scope, query_gene, gene_rankings, codep_method)
                    _reset_view_substate('codep')
                    match_info = {'targets': {'parsed': parsed, 'matched': [query_gene], 'missing': missing, 'resolved': resolved}}
                    st.session_state.codep_result = {
                        'mode': 'scan',
                        'method': codep_method,
                        'query_gene': query_gene,
                        'genes': [query_gene],
                        'scan_df': scan_df,
                        'rankings': gene_rankings,
                        'matches': match_info,
                        'draft': copy.deepcopy(draft),
                        'input_signature': input_signature(draft),
                        'metadata': {
                            **dataset_metadata,
                            'created_at_utc': datetime.now(timezone.utc).isoformat(),
                            'selected_genes': [query_gene],
                            'parameters': {'mode': 'scan', 'method': codep_method, 'scanned_genes': len(scan_df)},
                            'input_matching': match_info,
                        },
                    }
                    if is_dirty:
                        if info_note:
                            st.session_state['_flash_info'] = ('codep', info_note)
                        st.rerun()
                    elif info_note:
                        st.info(info_note)
            else:
                if len(matched) < 2:
                    st.warning(ui(
                        'Enter at least 2 matched genes for a correlation matrix (or switch to Genome-wide scan for 1 gene).',
                        '相关性矩阵至少需要 2 个匹配基因（单个基因请切换到“单基因全基因组扫描”）。',
                    ))
                else:
                    chosen = matched[:25]
                    info_note = None
                    if len(matched) > 25:
                        info_note = ui('Showing the first 25 matched genes in the correlation matrix.',
                                       '相关性矩阵最多展示前 25 个匹配基因。')
                    corr, counts = gene_correlation_matrix(df_scope, chosen, gene_rankings, method=codep_method)
                    _reset_view_substate('codep')
                    match_info = {'targets': {'parsed': parsed, 'matched': chosen, 'missing': missing, 'resolved': resolved}}
                    st.session_state.codep_result = {
                        'mode': 'matrix',
                        'method': codep_method,
                        'genes': chosen,
                        'corr': corr,
                        'counts': counts,
                        'rankings': gene_rankings,
                        'matches': match_info,
                        'draft': copy.deepcopy(draft),
                        'input_signature': input_signature(draft),
                        'metadata': {
                            **dataset_metadata,
                            'created_at_utc': datetime.now(timezone.utc).isoformat(),
                            'selected_genes': chosen,
                            'parameters': {'mode': 'matrix', 'method': codep_method},
                            'input_matching': match_info,
                        },
                    }
                    if is_dirty:
                        if info_note:
                            st.session_state['_flash_info'] = ('codep', info_note)
                        st.rerun()
                    elif info_note:
                        st.info(info_note)
        except (ValueError, UnicodeError) as exc:
            st.error(str(exc))

    result = st.session_state.get('codep_result')
    render_run_status('codep', draft, result)
    render_draft_feedback(active_text, None, draft, result)
    if result:
        match = result['matches']['targets']
        show_matches(match['parsed'], match['matched'], match['missing'], resolved=match.get('resolved'))
        render_codependency_result()


# =============================================================================
# 页脚：致谢 + 引用
# =============================================================================
st.markdown('<div class="section-gap"></div>', unsafe_allow_html=True)
with st.expander(t('resources'), expanded=False, icon=':material/menu_book:'):
    st.markdown(f"""
    <div class="footer-card">
        <h4>{t('acknowledgements')}</h4>
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
    f'font-size:0.79rem; padding:1.1rem 0 0.5rem;">'
    f'CRISPR Score Analyzer {TOOL_VERSION}{footer_release} · Changhao Kan · '
    f'<a href="{GITHUB_URL}" target="_blank" style="color:{get_theme()["accent"]}; text-decoration:none;">GitHub</a>'
    f'</div>',
    unsafe_allow_html=True
)
