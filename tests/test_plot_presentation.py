"""Presentation checks against the app's real plotting helpers, without network I/O."""

import ast
from html import escape
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pytest

from figure_export import prepare_export_figure


@pytest.fixture
def plot_helpers():
    names = {
        'get_theme', 'ui', 'apply_theme_to_fig', '_rank_label_annotations',
        'build_rank_figure', 'create_rank_plot', 'create_multilayer_rank_plot',
        'create_lineage_boxplot',
    }
    constants = {'THEMES', 'PLOT_COLORS', 'FONT_FAMILY', 'ESSENTIALITY_THRESHOLD'}
    tree = ast.parse((Path(__file__).resolve().parents[1] / 'app.py').read_text())
    nodes = [node for node in tree.body
             if (isinstance(node, ast.FunctionDef) and node.name in names)
             or (isinstance(node, ast.Assign)
                 and any(isinstance(target, ast.Name) and target.id in constants
                         for target in node.targets))]
    namespace = {'go': go, 'np': np, 'escape': escape, 'make_subplots': make_subplots,
                 'st': SimpleNamespace(session_state={'lang': 'en', 'theme': 'light'})}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), 'app.py', 'exec'), namespace)
    return namespace


@pytest.fixture
def rankings():
    count = 1000
    frame = pd.DataFrame({
        'gene': [f'GENE{i}' for i in range(count)],
        'rank': np.arange(1, count + 1),
        'mean_score': np.linspace(-3, 1, count),
        'n_valid': np.full(count, 12),
        'percentile': np.arange(1, count + 1) / count * 100,
    })
    frame.loc[497:500, 'gene'] = ['MYC', 'PTEN', 'E2F1', 'E2F3']
    frame['gene_upper'] = frame.gene.str.upper()
    return frame


@pytest.mark.parametrize('theme', ['light', 'dark'])
def test_close_labels_are_separated_and_arrows_preserve_true_positions(plot_helpers, rankings, theme):
    plot_helpers['st'].session_state['theme'] = theme
    make_plot = plot_helpers['create_rank_plot']
    fig = make_plot(rankings, ['E2F1', 'E2F3'], n_cell_lines=12)
    annotations = fig.layout.annotations
    assert {annotation.text for annotation in annotations} == {'MYC', 'PTEN', 'E2F1', 'E2F3'}
    y_span = fig.layout.yaxis.range[1] - fig.layout.yaxis.range[0]
    x_span = fig.layout.xaxis.range[1] - fig.layout.xaxis.range[0]
    for i, annotation in enumerate(annotations):
        original = rankings[rankings.gene == annotation.text].iloc[0]
        assert annotation.x == original['rank']
        assert annotation.y == original.mean_score
        assert annotation.showarrow
        assert annotation.axref == 'x' and annotation.ayref == 'y'
        assert annotation.font.size >= 13
        assert annotation.bgcolor == 'rgba(0,0,0,0)'
        assert annotation.borderwidth == 0
        for other in annotations[i + 1:]:
            # At the default plotting area, these short gene
            # labels need either a full text width or a line height of space.
            horizontal = abs(annotation.ax - other.ax) / x_span * (
                fig.layout.width - fig.layout.margin.l - fig.layout.margin.r)
            vertical = abs(annotation.ay - other.ay) / y_span * (
                fig.layout.height - fig.layout.margin.t - fig.layout.margin.b)
            assert horizontal >= 44 or vertical >= 23
    repeated = make_plot(rankings, ['E2F3', 'E2F1'], n_cell_lines=12)
    assert [a.to_plotly_json() for a in annotations] == [a.to_plotly_json() for a in repeated.layout.annotations]
    exported = prepare_export_figure(fig)
    assert all(trace.type == 'scatter' for trace in exported.data)
    assert [(a.x, a.y, a.ax, a.ay) for a in exported.layout.annotations] == [
        (a.x, a.y, a.ax, a.ay) for a in annotations]


@pytest.mark.parametrize('language, all_genes, reference, valid_count', [
    ('en', 'All genes', 'Reference · MYC', 'Valid cell lines'),
    ('zh', '全部基因', '参考基因 · MYC', '有效细胞系数'),
])
def test_localized_legends_and_labels_off_retain_gene_hover(
        plot_helpers, rankings, language, all_genes, reference, valid_count):
    plot_helpers['st'].session_state['lang'] = language
    fig = plot_helpers['create_rank_plot'](rankings, ['E2F1', 'E2F3'], show_labels=False)
    assert len(fig.layout.annotations) == 0
    assert fig.data[0].name == all_genes
    assert fig.data[1].name == reference
    assert fig.data[1].marker.symbol != fig.data[2].marker.symbol
    assert fig.data[2].marker.symbol != fig.data[3].marker.symbol
    assert all(trace.mode == 'markers' for trace in fig.data)
    assert all('%{text}' in trace.hovertemplate and valid_count in trace.hovertemplate for trace in fig.data)
    assert set(fig.data[-1].text) == {'E2F1', 'E2F3'}
    assert fig.layout.legend.font.size >= 12
    assert fig.layout.xaxis.tickfont.size >= 12


def test_large_background_stays_one_webgl_trace_without_duplicate_genes(plot_helpers):
    count = 18435
    frame = pd.DataFrame({
        'gene': ['MYC', 'PTEN', 'E2F1'] + [f'G{i}' for i in range(count - 3)],
        'rank': np.arange(1, count + 1), 'mean_score': np.linspace(-4, 0.5, count),
        'n_valid': np.full(count, 1186), 'percentile': np.arange(1, count + 1) / count * 100,
    })
    frame['gene_upper'] = frame.gene
    fig = plot_helpers['create_rank_plot'](frame, ['MYC', 'E2F1'])
    assert fig.data[0].type == 'scattergl'
    assert len(fig.data[0].x) == count - 3
    assert len(fig.data) == 3
    displayed = [gene for trace in fig.data for gene in trace.text]
    assert len(displayed) == len(set(displayed)) == count
    assert len(fig.layout.annotations) == 3


def test_lineage_plot_shows_readable_counts_and_localized_point_hover(plot_helpers):
    plot_helpers['st'].session_state['lang'] = 'zh'
    scores = pd.DataFrame({
        'gene': ['E2F1'] * 4, 'lineage': ['Liver'] * 2 + ['Lung'] * 2,
        'crispr_score': [-1.8, -1.2, -0.4, 0.1],
        'cell_line': ['ACH-1', 'ACH-2', 'ACH-3', 'ACH-4'],
        'cell_line_name': ['Model A', 'Model B', 'Model C', 'Model D'],
    })
    fig = plot_helpers['create_lineage_boxplot'](scores, ['E2F1'], show_points=True)
    assert fig.layout.yaxis.categoryarray == ('Liver (n=2)', 'Lung (n=2)')
    assert fig.layout.yaxis.tickfont.size >= 12
    assert fig.layout.annotations[0].font.size >= 15
    assert all(trace.boxpoints == 'all' and '细胞系' in trace.hovertemplate for trace in fig.data)
    assert fig.layout.xaxis.title.text == 'CRISPR 分数'
