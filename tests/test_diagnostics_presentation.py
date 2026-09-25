"""Diagnostics stay lazy while retaining actionable checks and complete downloads."""

import ast
import io
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from analysis_core import analyze_dataset


class Panel:
    def __init__(self, opened):
        self.open = opened

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass


class DiagnosticsUI:
    def __init__(self, *, opened=True, view=None, query=''):
        self.opened = opened
        self.view = view
        self.query = query
        self.frames = []
        self.downloads = {}
        self.captions = []
        self.messages = []

    def expander(self, *args, **kwargs):
        assert kwargs['on_change'] == 'rerun'
        return Panel(self.opened)

    def caption(self, text):
        self.captions.append(text)

    def selectbox(self, label, options, **kwargs):
        return self.view or options[0]

    def text_input(self, *args, **kwargs):
        return self.query

    def info(self, text):
        self.messages.append(text)

    def dataframe(self, frame, **kwargs):
        self.frames.append(frame.copy(deep=True))

    def download_button(self, label, data, filename, mime, **kwargs):
        assert kwargs['on_click'] == 'ignore'
        self.downloads[kwargs['key']] = data


def render(ui, diagnostics):
    source = (Path(__file__).resolve().parents[1] / 'app.py').read_text()
    function = next(node for node in ast.parse(source).body
                    if isinstance(node, ast.FunctionDef) and node.name == 'render_diagnostics')
    function.decorator_list = []
    namespace = {'st': ui, 'ui': lambda en, zh: en, 'diagnostics': diagnostics}
    exec(compile(ast.Module(body=[function], type_ignores=[]), 'app.py', 'exec'), namespace)
    namespace['render_diagnostics']()


@pytest.fixture
def diagnostics():
    frame = pd.DataFrame([
        ['ACH-1', 'A', -1., None, 'bad', np.inf, None, -.2, -.3],
        ['ACH-2', 'B', 0., 0., '1.25', -.5, None, -.4, -.6],
    ], columns=['depmap_id', 'cell_line_display_name', 'CLEAN (1)', 'MISSING (2)',
                'NONNUMERIC (3)', 'INFINITE (4)', 'EMPTY (5)', 'DUP (6)', 'DUP (7)'])
    return analyze_dataset(frame)[1]


def test_closed_panel_does_not_access_data_or_prepare_download(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail('A closed diagnostics panel must not serialize CSV data')

    monkeypatch.setattr(pd.DataFrame, 'to_csv', forbidden)
    ui = DiagnosticsUI(opened=False)
    render(ui, object())  # Any attempted indexing would fail before rendering.
    assert ui.frames == []
    assert ui.downloads == {}
    assert ui.captions == []


def test_default_view_keeps_all_problem_gene_types_and_hides_metadata(diagnostics, monkeypatch):
    original = diagnostics.copy(deep=True)

    def forbidden(*args, **kwargs):
        pytest.fail('Showing diagnostics must not eagerly serialize its download')

    monkeypatch.setattr(pd.DataFrame, 'to_csv', forbidden)
    ui = DiagnosticsUI()
    render(ui, diagnostics)
    assert len(ui.frames) == 1
    preview = ui.frames[0]
    assert preview['gene'].tolist() == ['MISSING', 'NONNUMERIC', 'INFINITE', 'EMPTY', 'DUP', 'DUP']
    assert preview.loc[preview.gene.eq('EMPTY'), 'reason'].item() == 'No valid scores'
    assert set(preview.loc[preview.gene.eq('DUP'), 'reason']) == {'Duplicate gene symbol'}
    assert 'position' not in preview.columns
    assert callable(ui.downloads['diagnostics_csv'])
    pd.testing.assert_frame_equal(diagnostics, original)


def test_other_columns_are_explained_but_download_retains_every_original_row(diagnostics, monkeypatch):
    serialized = []
    original_to_csv = pd.DataFrame.to_csv

    def tracked_to_csv(frame, *args, **kwargs):
        serialized.append(frame.copy(deep=True))
        return original_to_csv(frame, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, 'to_csv', tracked_to_csv)
    ui = DiagnosticsUI(view='other', query='DISPLAY_NAME')
    render(ui, diagnostics)
    assert ui.frames[0]['column'].tolist() == ['cell_line_display_name']
    assert ui.frames[0]['reason'].item() == 'Not recognized as a gene score'
    assert serialized == []

    content = ui.downloads['diagnostics_csv']()
    assert isinstance(content, bytes)
    downloaded = pd.read_csv(io.BytesIO(content), keep_default_na=False)
    assert downloaded['column'].tolist() == diagnostics['column'].tolist()
    assert downloaded['reason'].tolist() == diagnostics['reason'].tolist()
    assert downloaded.columns.tolist() == diagnostics.columns.tolist()
    assert len(serialized) == 1
    pd.testing.assert_frame_equal(serialized[0], diagnostics)


def test_preview_is_bounded_and_search_can_reach_rows_beyond_first_page():
    frame = pd.DataFrame([{f'GENE{i:03d} ({i + 1})': None for i in range(250)}])
    diagnostics = analyze_dataset(frame)[1]
    ui = DiagnosticsUI()
    render(ui, diagnostics)
    assert len(ui.frames[0]) == 200
    assert any('Showing 200 of 250 columns' in caption for caption in ui.captions)

    ui.query = 'gene249 ('  # Literal, case-insensitive source-column search.
    ui.frames.clear()
    render(ui, diagnostics)
    assert ui.frames[0]['gene'].tolist() == ['GENE249']
    downloaded = pd.read_csv(io.BytesIO(ui.downloads['diagnostics_csv']()))
    assert len(downloaded) == 250


def test_healthy_genes_remain_available_in_all_genes_view(diagnostics):
    ui = DiagnosticsUI(view='genes', query='clean')
    render(ui, diagnostics)
    assert ui.frames[0]['gene'].tolist() == ['CLEAN']
    assert ui.frames[0]['reason'].item() == 'Included in ranking'

    ui.query = 'not-a-column'
    ui.frames.clear()
    render(ui, diagnostics)
    assert ui.frames == []
    assert ui.messages == ['No columns to show for this selection.']
    assert callable(ui.downloads['diagnostics_csv'])
