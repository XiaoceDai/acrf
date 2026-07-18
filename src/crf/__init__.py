"""
Chamber-aware Rewiring Framework (CRF) — portable method core.

Re-exports portable CRF definitions (see `crf.schema`) and auditable
implementations in `scripts/fig2c_chamber_decomposition.py`, the v3
integrative panels `scripts/v3/fig5b_v3_triangulation_heatmap.py` /
`scripts/v3/fig5c_v3_cross_axis_consistency.py`, and the legacy ten-program
bridge `scripts/fig5a_triangulation.py` (not used in the v3 oncoprint).
"""

from __future__ import annotations

from . import schema
from .acai import acai, acai_rel, cai, interpret_acai_rel, ncp_baseline_sign
from .rewiring import (
    BULK_CLASS_TO_STRENGTH,
    bulk_strength_from_labels,
    dur_direction,
    fig2c_trend_state,
    FIG2C_EPS_FLOOR,
    FIG2C_EPS_REVERSED,
    FIG2C_RATIO_AMPLIFIED,
    FIG2C_RATIO_ATTENUATED,
)
from .triangulation import TriangulationInput, triangulate

__all__ = [
    "schema",
    "acai",
    "acai_rel",
    "cai",
    "interpret_acai_rel",
    "ncp_baseline_sign",
    "BULK_CLASS_TO_STRENGTH",
    "bulk_strength_from_labels",
    "dur_direction",
    "fig2c_trend_state",
    "FIG2C_EPS_FLOOR",
    "FIG2C_EPS_REVERSED",
    "FIG2C_RATIO_AMPLIFIED",
    "FIG2C_RATIO_ATTENUATED",
    "TriangulationInput",
    "triangulate",
    "__version__",
]

__version__ = "0.1.1"
