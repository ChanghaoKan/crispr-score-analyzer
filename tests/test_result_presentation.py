"""Browser formatting must preserve scientific values and downloadable provenance."""

import ast
import hashlib
import io
import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
import streamlit as st


class SessionState(dict):
    def __getattr__(self, name):
        return self[name]


class Block:
    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass


class ResultsUI:
    column_config = st.column_config

    def __init__(self):
        self.session_state = SessionState(lang='en', theme='light')
        self.frames = []
        self.downloads = {}
        self.clicked = set()

    def columns(self, widths, **kwargs):
        return [Block() for _ in widths]

    def popover(self, *args, **kwargs):
        return Block()

    def expander(self, *args, **kwargs):
        return Block()

    def spinner(self, *args, **kwargs):
        return Block()

    def markdown(self, *args, **kwargs):
        pass

    def caption(self, *args, **kwargs):
        pass

    def json(self, *args, **kwargs):
        pass

    def dataframe(self, frame, **kwargs):
        self.frames.append((frame.copy(deep=True), kwargs))

    def download_button(self, label, data, filename, mime, **kwargs):
        assert kwargs['on_click'] == 'ignore'
        self.downloads[kwargs['key']] = data

    def button(self, label, key, **kwargs):
        return key in self.clicked


@pytest.fixture
def presentation():
    ui = ResultsUI()
    bundles = []

    def bundle(rankings, genes, plot_data, metadata, diagnostics):
        bundles.append(metadata)
        return json.dumps(metadata, sort_keys=True).encode()

    namespace = {
        'st': ui,
        'ui': lambda en, zh: zh if ui.session_state.lang == 'zh' else en,
        'render_download_buttons': lambda *args, **kwargs: None,
        'make_analysis_bundle': bundle,
        'diagnostics': pd.DataFrame(),
        'essential_gene': 'MYC', 'nonessential_gene': 'PTEN',
        'show_labels': True, 'point_size': 4,
        'hashlib': hashlib, 'json': json,
    }
    source = (Path(__file__).resolve().parents[1] / 'app.py').read_text()
    function = next(node for node in ast.parse(source).body
                    if isinstance(node, ast.FunctionDef) and node.name == 'result_downloads')
    exec(compile(ast.Module(body=[function], type_ignores=[]), 'app.py', 'exec'), namespace)
    ranking = pd.DataFrame({
        'gene': ['MYC', 'PTEN'], 'gene_upper': ['MYC', 'PTEN'],
        'mean_score': [-1.23456789, .12345678], 'rank': [1, 2],
        'percentile': [50., 100.], 'n_valid': [3, 4], 'n_missing': [1, 0],
        'missing_fraction': [.25, 0.],
    })
    result = {'rankings': ranking, 'genes': ['PTEN', 'MYC'],
              'plot_data': pd.DataFrame({'gene': ['MYC'], 'crispr_score': [-1.23456789]}),
              'metadata': {'created_at_utc': '2026-09-25T01:02:03Z', 'selected_genes': ['PTEN', 'MYC']}}
    return ui, namespace['result_downloads'], result, SimpleNamespace(layout=SimpleNamespace(width=640, height=480)), bundles


def test_readable_localized_table_preserves_csv_values(presentation):
    ui, render, result, fig, bundles = presentation
    original = result['rankings'].copy(deep=True)
    for lang, gene_label, score_label in [('en', 'Gene', 'Mean score'), ('zh', '基因', '平均分')]:
        ui.session_state.lang = lang
        ui.frames.clear()
        render(result, fig, 'rank')
        frame, options = ui.frames[0]
        assert list(frame.columns) == ['gene', 'mean_score', 'rank', 'percentile', 'n_valid']
        assert frame['gene'].tolist() == ['PTEN', 'MYC']
        assert options['column_config']['gene']['label'] == gene_label
        assert options['column_config']['mean_score']['label'] == score_label
        assert options['column_config']['mean_score']['type_config']['format'] == '%.3f'
        assert options['column_config']['percentile']['type_config']['format'] == '%.2f%%'
        csv = pd.read_csv(io.StringIO(ui.downloads['rank_csv']))
        assert csv['gene'].tolist() == ['PTEN', 'MYC']
        assert csv['mean_score'].tolist() == [.12345678, -1.23456789]
        assert csv['missing_fraction'].tolist() == [0., .25]
        assert 'gene_upper' not in csv.columns
        assert result['rankings'].equals(original)
    assert bundles == []


def test_zip_is_prepared_on_demand_and_invalidated_for_new_analysis_or_dimensions(presentation):
    ui, render, result, fig, bundles = presentation
    render(result, fig, 'rank')
    assert bundles == []
    assert 'rank_zip' not in ui.downloads
    ui.clicked = {'rank_zip_btn'}
    render(result, fig, 'rank')
    assert len(bundles) == 1
    assert bundles[0]['export']['width'] == 640
    assert bundles[0]['export']['height'] == 480
    first = ui.downloads['rank_zip']
    ui.clicked.clear()
    render(result, fig, 'rank')
    assert len(bundles) == 1
    assert ui.downloads['rank_zip'] == first

    for mutation in [
        lambda: result['metadata'].update(created_at_utc='2026-09-25T02:03:04Z'),
        lambda: ui.session_state.update(rank_width=1500),
    ]:
        mutation()
        ui.downloads.clear()
        render(result, fig, 'rank')
        assert 'rank_zip' not in ui.downloads
        assert '_export_rank_zip' not in ui.session_state
        ui.clicked = {'rank_zip_btn'}
        render(result, fig, 'rank')
        assert ui.downloads['rank_zip'] != first
        first = ui.downloads['rank_zip']
        ui.clicked.clear()
    assert len(bundles) == 3
