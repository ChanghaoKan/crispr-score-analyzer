"""Deterministic, publication-ready figure exports without UI state."""

import copy
import hashlib
import io
import json
import math

import plotly.graph_objects as go
from PIL import Image


_EXPORT_VERSION = 3
MAX_PNG_PIXELS = 24_000_000
_INK = "#1a1a1a"


def _export_options(format, width, height, scale):
    format = str(format).lower()
    if format not in {"pdf", "svg", "png"}:
        raise ValueError("Export format must be PDF, SVG, or PNG.")
    for name, value in (("width", width), ("height", height), ("scale", scale)):
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"Export {name} must be a positive number.")
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"Export {name} must be a positive number.")
    if int(width) != width or int(height) != height:
        raise ValueError("Export width and height must be whole pixels.")
    return format, int(width), int(height), float(scale)


def export_fingerprint(fig, format, width=1000, height=600, scale=3):
    """Identify the figure and every option that can change its export bytes."""
    format, width, height, scale = _export_options(format, width, height, scale)
    payload = {
        "export_version": _EXPORT_VERSION,
        "figure": json.loads(fig.to_json()),
        "format": format,
        "width": width,
        "height": height,
        "scale": scale,
        "png_dpi": 300,
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def prepare_export_figure(fig, full_vector=True):
    """Copy a figure onto white paper, optionally converting WebGL scatter traces."""
    payload = copy.deepcopy(fig.to_plotly_json())
    if full_vector:
        trace_groups = [payload.get("data", [])]
        trace_groups.extend(frame.get("data", []) for frame in payload.get("frames", []))
        for traces in trace_groups:
            for trace in traces:
                if trace.get("type") == "scattergl":
                    trace["type"] = "scatter"
    export_fig = go.Figure(payload)
    export_fig.update_layout(
        plot_bgcolor="white",
        paper_bgcolor="white",
        font=dict(color=_INK),
        title=dict(font=dict(color=_INK)),
        legend=dict(
            font=dict(color=_INK),
            title=dict(font=dict(color=_INK)),
            bgcolor="rgba(0,0,0,0)",
        ),
    )
    axis_style = dict(
        linecolor=_INK,
        tickcolor=_INK,
        tickfont=dict(color=_INK),
        title_font=dict(color=_INK),
        zerolinecolor="#667085",
        gridcolor="#e5e7eb",
    )
    export_fig.update_xaxes(**axis_style)
    export_fig.update_yaxes(**axis_style)
    export_fig.update_annotations(font=dict(color=_INK), arrowcolor=_INK,
                                  bgcolor="rgba(0,0,0,0)", borderwidth=0)
    for shape in export_fig.layout.shapes:
        shape.line.color = "#667085"
    for trace in export_fig.data:
        if trace.name in {"All genes", "全部基因"}:
            trace.marker.color = "rgba(150,150,150,0.55)"
        if trace.type == "box":
            trace.line.color = "#3f4b55"
        if hasattr(trace, "textfont"):
            trace.textfont.color = _INK
    return export_fig


def render_figure_bytes(fig, format, width=1000, height=600, scale=3):
    """Render vectors or a high-resolution PNG with explicit 300-DPI metadata."""
    format, width, height, scale = _export_options(format, width, height, scale)
    if format == "png" and width * height * scale ** 2 > MAX_PNG_PIXELS:
        raise ValueError("PNG exceeds 24 million pixels. Reduce width or height, or export PDF/SVG.")
    export_fig = prepare_export_figure(fig, full_vector=format in {"pdf", "svg"})
    kwargs = {"format": format, "width": width, "height": height}
    if format == "png":
        kwargs["scale"] = scale
    image_bytes = export_fig.to_image(**kwargs)
    if format != "png":
        return image_bytes
    with Image.open(io.BytesIO(image_bytes)) as image:
        output = io.BytesIO()
        image.save(output, format="PNG", dpi=(300, 300))
    return output.getvalue()
