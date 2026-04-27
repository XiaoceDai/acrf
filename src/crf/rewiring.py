"""
Rewiring labels and continuous→discrete rules (Fig 2C + Fig 5A bulk leg).

- Fig 2C state machine: `scripts/fig2c_chamber_decomposition.py` ``classify_state``.
- Fig 5A bulk strength: `scripts/fig5a_triangulation.py` ``bulk_score``.
"""

from __future__ import annotations

import math
from typing import Any

__all__ = [
    "FIG2C_EPS_REVERSED",
    "FIG2C_EPS_FLOOR",
    "FIG2C_RATIO_ATTENUATED",
    "FIG2C_RATIO_AMPLIFIED",
    "BULK_CLASS_TO_STRENGTH",
    "bulk_strength_from_labels",
    "dur_direction",
    "fig2c_trend_state",
]

# —— Fig 2C (trend / exploratory) — same constants as fig2c_chamber_decomposition.py
FIG2C_EPS_REVERSED = 0.02
FIG2C_EPS_FLOOR = 0.02
FIG2C_RATIO_ATTENUATED = 0.67
FIG2C_RATIO_AMPLIFIED = 1.33

# —— Fig 5A: consensus / trend string → signed strength in {-2,-1,0,+1,+2} ——
BULK_CLASS_TO_STRENGTH: dict[str, int] = {
    "Reversed": -2,
    "Reversed-weak": -1,
    "Attenuated": -1,
    "Attenuated-weak": -1,
    "Preserved": 0,
    "": 0,
    "Low-support": 0,
    "Amplified-weak": 1,
    "Amplified": 2,
}


def _strip(s: Any) -> str:
    if s is None or (isinstance(s, float) and math.isnan(s)):
        return ""
    t = str(s).strip().strip('"')
    return t


def bulk_strength_from_labels(
    consensus_class: Any,
    trend_class: Any,
) -> tuple[int, str]:
    """
    Map ``consensus_class`` (permutation-gated) with ``trend_class`` fallback
    to a signed ``bulk_strength`` and resolved label.
    """
    c = _strip(consensus_class)
    t = _strip(trend_class)
    label = c or t
    strength = BULK_CLASS_TO_STRENGTH.get(label, 0)
    return strength, label or "n/a"


def dur_direction(rho: float | None, thr: float = 0.30) -> int:
    """LAA primary duration trend: sign of Spearman ρ if |ρ| >= thr, else 0 (Fig 5A)."""
    if rho is None:
        return 0
    r = float(rho)
    if math.isnan(r):
        return 0
    if r <= -thr:
        return -1
    if r >= thr:
        return 1
    return 0


def fig2c_trend_state(af_aligned: float, ref_aligned: float) -> str:
    """
    Preserved / Attenuated / Reversed / Amplified from *aligned* AF vs SR means.

    ``af_aligned`` and ``ref_aligned`` must already be signed by program
    direction (see Fig 2C pipeline: ``aligned_AF``, ``reference_aligned_SR``).
    """
    if af_aligned <= -FIG2C_EPS_REVERSED:
        return "Reversed"
    ref = max(float(ref_aligned), FIG2C_EPS_FLOOR)
    ratio = float(af_aligned) / ref
    if ratio < FIG2C_RATIO_ATTENUATED:
        return "Attenuated"
    if ratio <= FIG2C_RATIO_AMPLIFIED:
        return "Preserved"
    return "Amplified"
