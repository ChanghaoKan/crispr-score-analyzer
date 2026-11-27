"""Presentation checks against the core plotting library, without Streamlit or network I/O."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from crispr_score_analyzer import plotting
from crispr_score_analyzer.figure_export import prepare_export_figure


@pytest.fixture
def plot_helpers():
    st_stub = SimpleNamespace(session_state={'lang': 'en', 'theme': 'light'})

    def _state():
        return {'theme': st_stub.session_state['theme'], 'lang': st_stub.session_state['lang']}

    return {
        'st': st_stub,
        'create_rank_plot': lambda *a, **kw: plotting.create_rank_plot(*a, **{**_state(), **kw}),
        'create_multilayer_rank_plot': lambda *a, **kw: plotting.create_multilayer_rank_plot(*a, **{**_state(), **kw}),
        'create_lineage_boxplot': lambda *a, **kw: plotting.create_lineage_boxplot(*a, **{**_state(), **kw}),
        'create_preference_plot': lambda *a, **kw: plotting.create_preference_plot(*a, **{**_state(), **kw}),
        'create_preference_heatmap': lambda *a, **kw: plotting.create_preference_heatmap(*a, **{**_state(), **kw}),
        'create_correlation_heatmap': lambda *a, **kw: plotting.create_correlation_heatmap(*a, **{**_state(), **kw}),
    }


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


def test_point_size_reuses_the_background_curve(plot_helpers, rankings):
    small = plot_helpers['create_rank_plot'](rankings, ['E2F1'], point_size=3)
    large = plot_helpers['create_rank_plot'](rankings, ['E2F1'], point_size=8)
    assert list(small.data[0].x) == list(large.data[0].x)
    assert small.data[-1].marker.size != large.data[-1].marker.size
    added = plot_helpers['create_rank_plot'](rankings, ['E2F1', 'GENE1'], point_size=8)
    assert len(added.data[0].x) == len(small.data[0].x) - 1


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


def test_preference_plot_is_a_labeled_lollipop_per_gene(plot_helpers):
    plot_helpers['st'].session_state['theme'] = 'dark'
    table = pd.DataFrame({
        'gene': ['MYC', 'MYC', 'MYC', 'PTEN', 'PTEN', 'PTEN'],
        'lineage': ['Lung', 'Liver', 'Bone', 'Lung', 'Liver', 'Bone'],
        'n_valid': [12, 8, 3, 12, 8, 9],
        'mean': [-2.2, -1.4, -1.0, 0.1, 0.4, 0.2],
        'panel_mean': [-1.8, -1.8, -1.8, 0.3, 0.3, 0.3],
        'difference': [-0.4, 0.4, 0.8, -0.2, 0.1, -0.1],
    })
    fig = plot_helpers['create_preference_plot'](table)
    assert fig is not None
    assert all(trace.type == 'scatter' for trace in fig.data)
    assert [annotation.text for annotation in fig.layout.annotations] == [
        'MYC · all-type mean -1.80', 'PTEN · all-type mean 0.30']
    # Bone has only 3 lines, below the plotting minimum, so MYC shows Lung then Liver.
    assert list(fig.layout.yaxis.categoryarray) == ['Liver', 'Lung']
    assert list(fig.layout.yaxis2.categoryarray) == ['Liver', 'Bone', 'Lung']
    assert fig.layout.xaxis.range[0] < 0 < fig.layout.xaxis.range[1]
    assert abs(fig.layout.xaxis.range[1]) > abs(fig.layout.xaxis2.range[1])
    marker_text = [trace.text for trace in fig.data if trace.mode == 'markers+text']
    assert marker_text[0] == ('+0.40', '-0.40')
    assert 'more dependent' in fig.layout.xaxis.title.text


def test_preference_top_k_and_heatmap_highlight_active_cohort(plot_helpers):
    rows = []
    for idx in range(16):
        diff = -0.8 + idx * 0.1
        rows.append({
            'gene': 'MYC', 'lineage': f'Type{idx:02d}', 'n_valid': 10,
            'mean': -1.0 + diff, 'panel_mean': -1.0, 'difference': diff,
            'median_difference': diff * 0.9, 'dependent_fraction': 0.8 if diff < 0 else 0.1,
        })
    table = pd.DataFrame(rows)
    lollipop = plot_helpers['create_preference_plot'](table, highlight_lineages=['Type07'], top_k=3)
    assert lollipop is not None
    # 3 lowest + 3 highest + 1 highlighted middle lineage = 7 categories.
    assert len(lollipop.layout.yaxis.categoryarray) == 7
    assert 'Type07' in lollipop.layout.yaxis.categoryarray
    marker_trace = [t for t in lollipop.data if t.mode == 'markers+text'][0]
    assert 'Median difference' in marker_trace.hovertemplate

    heatmap = plot_helpers['create_preference_heatmap'](table, highlight_lineages=['Type07'])
    assert heatmap is not None
    assert heatmap.data[0].type == 'heatmap'
    assert '★ Type07' in list(heatmap.data[0].y)


def test_correlation_heatmap_and_boxplot_modes(plot_helpers):
    corr = pd.DataFrame([[1.0, 0.82], [0.82, 1.0]], index=['MYC', 'MAX'], columns=['MYC', 'MAX'])
    counts = pd.DataFrame([[100, 98], [98, 100]], index=['MYC', 'MAX'], columns=['MYC', 'MAX'])
    hfig = plot_helpers['create_correlation_heatmap'](corr, counts, method='spearman')
    assert hfig.data[0].type == 'heatmap'
    assert 'Spearman' in hfig.data[0].hovertemplate

    scores = pd.DataFrame({
        'gene': ['E2F1'] * 4 + ['E2F2'] * 4,
        'lineage': ['Liver', 'Liver', 'Lung', 'Lung'] * 2,
        'crispr_score': [-1.6, -1.2, -0.3, 0.1, -0.8, -0.6, 0.2, 0.4],
        'cell_line': ['ACH-1', 'ACH-2', 'ACH-3', 'ACH-4'] * 2,
        'cell_line_name': ['A', 'B', 'C', 'D'] * 2,
    })
    grid_fig = plot_helpers['create_lineage_boxplot'](scores, ['E2F1', 'E2F2'], columns=2)
    assert grid_fig.layout.xaxis2 is not None
    by_gene = plot_helpers['create_lineage_boxplot'](scores, ['E2F1', 'E2F2'], compare_mode='by_gene')
    assert len(by_gene.data) == 2
    assert 'E2F1 (n=4)' in by_gene.layout.yaxis.categoryarray

