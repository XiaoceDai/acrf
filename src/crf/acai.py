"""LA-RA contrast, then the healthy-direction sign and the healthy magnitude."""

from __future__ import annotations

import math


def cai(la_score, ra_score):
    return float(la_score) - float(ra_score)


def ncp_baseline_sign(direction):
    # Same sign the panel uses: LA-biased programs stay positive.
    if direction == "LA-biased":
        return 1
    if direction == "RA-biased":
        return -1
    raise ValueError(f"Unknown program direction: {direction!r}")


def acai(cai_value, baseline_direction):
    if baseline_direction not in (-1, 1):
        raise ValueError("baseline_direction must be -1 (RA) or +1 (LA)")
    return float(cai_value) * float(baseline_direction)


def acai_rel(acai_value, healthy_bulk_magnitude, *, min_denominator=None):
    d = float(healthy_bulk_magnitude)
    if min_denominator is not None and abs(d) < min_denominator:
        if d == 0.0:
            d = float(min_denominator)
        else:
            d = math.copysign(float(min_denominator), d)
    if d == 0.0:
        raise ZeroDivisionError(
            "healthy_bulk_magnitude is zero; supply min_denominator or filter"
        )
    return float(acai_value) / abs(d)


def interpret_acai_rel(value):
    v = float(value)
    if v > 1.0:
        return "amplified"
    if v >= 0.0:
        return "attenuated"
    return "reversed"
