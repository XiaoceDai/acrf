"""Fig 5A class. Same branches as scripts/fig5a_triangulation.py triangulate()."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class TriangulationInput:
    bulk_strength: int
    dur_dir: int
    baseline_direction: int
    af_major_carrier_ct: str
    healthy_top_expressor_ct: str = ""


def triangulate(row):
    bulk = int(row.bulk_strength)
    dd = int(row.dur_dir)
    is_metab = row.baseline_direction == 1
    is_struct = row.baseline_direction == -1
    loc = (row.af_major_carrier_ct or "").lower()

    # Metabolism-led
    if is_metab and dd == -1:
        return "Metabolism-led"

    # Structure-led
    if is_struct and bulk <= -1:
        if "fibroblast" in loc:
            return "Structure-led"
        # Fig 1 expressor, used when the AF carrier cell type is missing.
        if "fibroblast" in (row.healthy_top_expressor_ct or "").lower():
            return "Structure-led"

    # Convergent: bulk and duration agree on weakening/decline.
    # LA-high: AF side declines (bulk -1/-2 and dur -1).
    # RA-high in the LAA arm: AF score may rise (bulk -1/-2 and dur +1).
    if bulk != 0 and dd != 0:
        if (is_metab and bulk <= -1 and dd == -1) or (is_struct and bulk <= -1 and dd == 1):
            return "Convergent"

    return "Fragmented"
