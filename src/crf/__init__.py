"""Installable method core. Scoring lives in the figure scripts; this package repeats those rules."""

from __future__ import annotations

from . import schema
from .acai import acai, acai_rel, cai, interpret_acai_rel, ncp_baseline_sign
from .rewiring import (
    BULK_CLASS_TO_STRENGTH,
    FIG2C_EPS_FLOOR,
    FIG2C_EPS_REVERSED,
    FIG2C_RATIO_AMPLIFIED,
    FIG2C_RATIO_ATTENUATED,
    bulk_strength_from_labels,
    dur_direction,
    fig2c_trend_state,
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

__version__ = "0.2.0"
