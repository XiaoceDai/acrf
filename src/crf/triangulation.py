"""
Figure 5A triangulation class (Convergent / Structure-led / Metabolism-led / Fragmented).

Rules are locked to match ``scripts/fig5a_triangulation.py`` `triangulate` + docstring.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

TriClass = Literal["Metabolism-led", "Structure-led", "Convergent", "Fragmented"]


@dataclass(frozen=True)
class TriangulationInput:
    """One program row: bulk strength, duration direction, NCP arm, localisation text."""

    bulk_strength: int
    """Signed strength in {-2,-1,0,+1,+2} from :func:`crf.rewiring.bulk_strength_from_labels`."""
    dur_dir: int
    """-1, 0, or +1 from :func:`crf.rewiring.dur_direction`."""
    baseline_direction: int
    """+1 = LA-biased (metabolic arm), -1 = RA-biased (structural arm)."""
    af_major_carrier_ct: str
    """e.g. fibroblast; may be "Unresolved"."""
    healthy_top_expressor_ct: str = ""
    """Optional Fig 1 expressor-based fallback for structure-led."""

    def is_metab(self) -> bool:
        return self.baseline_direction == 1

    def is_struct(self) -> bool:
        return self.baseline_direction == -1


def _loc_lower(s: str) -> str:
    return (s or "").lower()


def triangulate(row: TriangulationInput) -> TriClass:
    bulk = int(row.bulk_strength)
    dd = int(row.dur_dir)
    loc = _loc_lower(row.af_major_carrier_ct)
    fib_fallback = "fibroblast" in _loc_lower(row.healthy_top_expressor_ct)

    # Metabolism-led
    if row.is_metab() and dd == -1:
        return "Metabolism-led"

    # Structure-led
    if row.is_struct() and bulk <= -1:
        if "fibroblast" in loc:
            return "Structure-led"
        if fib_fallback:
            return "Structure-led"

    # Convergent
    if bulk != 0 and dd != 0:
        if (row.is_metab() and bulk <= -1 and dd == -1) or (row.is_struct() and bulk <= -1 and dd == 1):
            return "Convergent"

    return "Fragmented"
