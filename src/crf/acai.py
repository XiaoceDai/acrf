"""
Aligned chamber asymmetry: CAI, aCAI, aCAI_rel (see `crf.schema` and package docstring).
"""

from __future__ import annotations

import math
from typing import Literal

__all__ = [
    "cai",
    "acai",
    "acai_rel",
    "ncp_baseline_sign",
    "interpret_acai_rel",
]


def cai(la_score: float, ra_score: float) -> float:
    """Chamber asymmetry: score(LA) − score(RA) for one donor × program."""
    return float(la_score) - float(ra_score)


def ncp_baseline_sign(direction: Literal["LA-biased", "RA-biased"]) -> int:
    """
    NCP direction to alignment sign s_p: LA-biased → +1, RA-biased → -1
    (healthy bulk is the reference for each program’s stable bias).
    """
    if direction == "LA-biased":
        return 1
    if direction == "RA-biased":
        return -1
    raise ValueError(f"Unknown program direction: {direction!r}")


def acai(cai_value: float, baseline_direction: int) -> float:
    """
    aCAI = s_p · CAI with baseline_direction ∈ {+1, -1} from NCP.
    ``baseline_direction`` matches ``scripts`` / panel column ``baseline_direction``.
    """
    if baseline_direction not in (-1, 1):
        raise ValueError("baseline_direction must be -1 (RA) or +1 (LA)")
    return float(cai_value) * float(baseline_direction)


def acai_rel(
    acai_value: float,
    healthy_bulk_magnitude: float,
    *,
    min_denominator: float | None = None,
) -> float:
    """
    aCAI_rel = aCAI / |aCAI^healthy_bulk(p)|.

    R pipelines use the cohort healthy reference magnitude; optional floor avoids
    division blow-ups when the reference is near zero (set per analysis policy).
    """
    d = float(healthy_bulk_magnitude)
    if min_denominator is not None and abs(d) < min_denominator:
        d = math.copysign(float(min_denominator), d) if d != 0.0 else float(min_denominator)
    if d == 0.0:
        raise ZeroDivisionError("healthy_bulk_magnitude is zero; supply min_denominator or filter")
    return float(acai_value) / abs(d)


def interpret_acai_rel(value: float) -> str:
    """Framework interpretation table for aCAI_rel (main text, §2.2)."""
    v = float(value)
    if v > 1.0:
        return "amplified"
    if v >= 0.0:
        return "attenuated"
    return "reversed"
