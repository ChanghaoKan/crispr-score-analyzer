"""Backward-compatible shim re-exporting ``crispr_score_analyzer.analysis_core``."""

from crispr_score_analyzer import analysis_core as _impl
from crispr_score_analyzer.analysis_core import *  # noqa: F401, F403

for _name in dir(_impl):
    if not (_name.startswith("__") and _name.endswith("__")):
        globals()[_name] = getattr(_impl, _name)
del _name, _impl
