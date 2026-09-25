import ast
import io
from pathlib import Path
import xml.etree.ElementTree as ET

import plotly.graph_objects as go
import pytest
from PIL import Image

from figure_export import export_fingerprint, prepare_export_figure, render_figure_bytes, MAX_PNG_PIXELS


@pytest.fixture
def figure():
    fig = go.Figure([
        go.Scattergl(
            x=[1, 2, 3], y=[-1.5, -0.7, 0.1],
            text=["MYC", "KIF18A", "PTEN"], name="All genes",
            mode="markers", marker=dict(color="#aaaaaa", size=4),
        ),
        go.Scatter(
            x=[2], y=[-0.7], text=["KIF18A"], name="Target",
            mode="markers+text", marker=dict(color="#0072B2", size=9),
            textfont=dict(color="white", family="Arial", size=13),
        ),
    ])
    fig.update_layout(
        template="plotly_dark", paper_bgcolor="#111111", plot_bgcolor="#111111",
        font=dict(color="white", family="Arial", size=12),
        title=dict(text="Dependency", font=dict(color="white", size=18)),
        legend=dict(font=dict(color="white")),
        xaxis=dict(title="Rank", title_font=dict(color="white")),
        yaxis=dict(title="Gene effect", tickfont=dict(color="white")),
    )
    fig.add_annotation(
        text="<i>KIF18A</i>", x=2, y=-0.7, showarrow=True,
        font=dict(color="white", size=15, family="Arial"), bgcolor="#182430",
    )
    fig.add_hline(y=-0.5, line_dash="dash", line_color="white")
    return fig


def test_fingerprint_changes_for_plot_and_export_options(figure):
    baseline = export_fingerprint(figure, "png")
    assert len(baseline) == 64
    assert baseline == export_fingerprint(go.Figure(figure), "PNG")
    variants = []
    for mutate in (
        lambda fig: setattr(fig.data[1], "text", ["MYC"]),
        lambda fig: setattr(fig.data[1], "y", [-1.1]),
        lambda fig: setattr(fig.data[1].marker, "color", "#D55E00"),
        lambda fig: setattr(fig.data[1].marker, "size", 12),
    ):
        changed = go.Figure(figure)
        mutate(changed)
        variants.append(export_fingerprint(changed, "png"))
    variants.extend([
        export_fingerprint(figure, "png", height=700),
        export_fingerprint(figure, "png", width=1200),
        export_fingerprint(figure, "png", scale=2),
        export_fingerprint(figure, "svg"),
    ])
    assert baseline not in variants
    assert len(set(variants)) == len(variants)


def test_export_copy_is_vector_readable_and_preserves_formatting(figure):
    original = figure.to_json()
    exported = prepare_export_figure(figure)
    assert [trace.type for trace in exported.data] == ["scatter", "scatter"]
    assert figure.to_json() == original
    assert figure.data[0].type == "scattergl"
    assert exported.layout.paper_bgcolor == "white"
    assert exported.layout.plot_bgcolor == "white"
    assert exported.layout.font.family == "Arial"
    assert exported.layout.font.color == "#1a1a1a"
    assert exported.layout.title.font.color == "#1a1a1a"
    assert exported.layout.legend.font.color == "#1a1a1a"
    assert exported.layout.xaxis.title.font.color == "#1a1a1a"
    assert exported.layout.yaxis.tickfont.color == "#1a1a1a"
    annotation = exported.layout.annotations[0]
    assert annotation.text == "<i>KIF18A</i>"
    assert annotation.font.family == "Arial"
    assert annotation.font.size == 15
    assert annotation.font.color == "#1a1a1a"
    assert annotation.bgcolor == "rgba(0,0,0,0)"
    assert annotation.borderwidth == 0
    assert figure.layout.annotations[0].bgcolor == "#182430"
    assert exported.data[1].textfont.size == 13
    assert exported.data[1].marker.color == "#0072B2"
    assert exported.layout.shapes[0].line.dash == "dash"
    assert prepare_export_figure(figure, full_vector=False).data[0].type == "scattergl"


@pytest.mark.parametrize("name", ["All genes", "全部基因"])
def test_background_points_export_consistently_in_both_languages(figure, name):
    figure.data[0].name = name
    figure.data[0].marker.color = "rgba(170,183,195,0.22)"
    exported = prepare_export_figure(figure)
    assert exported.data[0].name == name
    assert exported.data[0].marker.color == "rgba(150,150,150,0.55)"
    assert figure.data[0].marker.color == "rgba(170,183,195,0.22)"
    assert exported.data[1].marker.color == figure.data[1].marker.color


@pytest.mark.parametrize("options", [
    {"format": "jpeg"}, {"format": "png", "width": 0},
    {"format": "png", "height": -1}, {"format": "png", "scale": float("nan")},
    {"format": "png", "width": 100.5}, {"format": "png", "scale": True},
])
def test_reject_invalid_export_options(figure, options):
    with pytest.raises(ValueError):
        export_fingerprint(figure, **options)
    with pytest.raises(ValueError):
        render_figure_bytes(figure, **options)


def test_actual_svg_is_vector_pdf_is_valid_and_png_has_dpi(figure):
    original = figure.to_json()
    svg = render_figure_bytes(figure, "svg", width=360, height=240)
    svg_root = ET.fromstring(svg)
    assert svg_root.tag.endswith("svg")
    assert not any(element.tag.rsplit("}", 1)[-1] == "image" for element in svg_root.iter())
    assert b"KIF18A" in svg
    annotation_backgrounds = [element for group in svg_root.iter()
                              if group.get("class") == "annotation"
                              for element in group.iter()
                              if element.tag.endswith("rect") and element.get("class") == "bg"]
    assert annotation_backgrounds
    assert all("fill-opacity: 0;" in element.get("style", "")
               for element in annotation_backgrounds)
    pdf = render_figure_bytes(figure, "pdf", width=360, height=240)
    assert pdf.startswith(b"%PDF-")
    assert b"%%EOF" in pdf[-1024:]
    png = render_figure_bytes(figure, "png", width=360, height=240, scale=3)
    with Image.open(io.BytesIO(png)) as image:
        assert image.size == (1080, 720)
        assert image.info["dpi"] == pytest.approx((300, 300), abs=0.02)
        assert image.convert("RGB").getpixel((0, 0)) == (255, 255, 255)
    assert figure.to_json() == original


def test_oversized_png_rejected_before_rendering_but_vectors_allowed(figure, monkeypatch):
    calls = []
    def render(self, **options):
        calls.append(options)
        return b"vector-output"
    monkeypatch.setattr(go.Figure, "to_image", render)
    with pytest.raises(ValueError, match="24 million pixels"):
        render_figure_bytes(figure, "png", width=3000, height=12000)
    assert calls == []
    assert render_figure_bytes(figure, "svg", width=3000, height=12000) == b"vector-output"
    assert len(calls) == 1


def test_app_download_cache_replaces_bytes_when_plot_or_dimensions_change(figure):
    """Exercise the actual app helper without importing or launching Streamlit."""
    class RerunSignal(BaseException):
        pass

    class Block:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def number_input(self, label, minimum, maximum, value, step, key):
            return state.dimensions.get(key, value)

    class FakeStreamlit:
        session_state = {"_export_box_pdf": {"signature": "other-view", "bytes": b"other-view"}}
        dimensions = {}
        clicked = set()
        downloads = []

        def columns(self, count):
            return [Block() for _ in range(count)]

        def caption(self, *args):
            pass

        def button(self, label, key, **kwargs):
            return key in self.clicked

        def spinner(self, label):
            return Block()

        def expander(self, label, **kwargs):
            return Block()

        def download_button(self, label, data, filename, mime, **kwargs):
            assert kwargs["on_click"] == "ignore"
            self.downloads.append((filename, data))

        def rerun(self):
            raise RerunSignal

        def error(self, message):
            pytest.fail(message)

    state = FakeStreamlit()
    namespace = {
        "st": state,
        "ui": lambda en, zh: en,
        "export_fingerprint": export_fingerprint,
        "MAX_PNG_PIXELS": MAX_PNG_PIXELS,
        "render_figure_bytes": lambda fig, fmt, width, height: export_fingerprint(
            fig, fmt, width, height
        ).encode(),
    }
    source = (Path(__file__).resolve().parents[1] / "app.py").read_text()
    function = next(node for node in ast.parse(source).body
                    if isinstance(node, ast.FunctionDef) and node.name == "render_download_buttons")
    exec(compile(ast.Module(body=[function], type_ignores=[]), "app.py", "exec"), namespace)
    render_buttons = namespace["render_download_buttons"]

    state.clicked = {"rank_pdf_btn"}
    with pytest.raises(RerunSignal):
        render_buttons(figure, "rank", "rank")
    first_bytes = state.session_state["_export_rank_pdf"]["bytes"]
    state.clicked.clear()
    render_buttons(figure, "rank", "rank")
    assert state.downloads == [("rank.pdf", first_bytes)]

    changed = go.Figure(figure)
    changed.data[1].text = ["MYC"]
    state.downloads.clear()
    render_buttons(changed, "rank", "rank")
    assert state.downloads == []
    assert "_export_rank_pdf" not in state.session_state
    state.clicked = {"rank_pdf_btn"}
    with pytest.raises(RerunSignal):
        render_buttons(changed, "rank", "rank")
    second_bytes = state.session_state["_export_rank_pdf"]["bytes"]
    assert second_bytes != first_bytes

    state.clicked.clear()
    state.dimensions["rank_height"] = 800
    render_buttons(changed, "rank", "rank")
    assert state.downloads == []
    assert "_export_rank_pdf" not in state.session_state
    state.clicked = {"rank_pdf_btn"}
    with pytest.raises(RerunSignal):
        render_buttons(changed, "rank", "rank")
    assert state.session_state["_export_rank_pdf"]["bytes"] not in {first_bytes, second_bytes}
    assert set(state.session_state) == {"_export_rank_pdf", "_export_box_pdf"}
    assert state.session_state["_export_box_pdf"]["bytes"] == b"other-view"
