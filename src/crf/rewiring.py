"""Fig 2C state labels and the Fig 5A bulk/duration scores.

Thresholds match scripts/ncp_v3_singscore_reanalysis/regenerate_fig2.py
and scripts/fig5a_triangulation.py. Do not retune them here.
"""

from __future__ import annotations

import math

# Same cutoffs as classify_state() in the fig2 redraw.
FIG2C_EPS_REVERSED = 0.02
FIG2C_EPS_FLOOR = 0.02
FIG2C_RATIO_ATTENUATED = 0.67
FIG2C_RATIO_AMPLIFIED = 1.33

BULK_CLASS_TO_STRENGTH = {
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


def _strip(s):
    if s is None or (isinstance(s, float) and math.isnan(s)):
        return ""
    return str(s).strip().strip('"')


def bulk_strength_from_labels(consensus_class, trend_class):
    # Permutation label if present, otherwise the trend label.
    label = _strip(consensus_class) or _strip(trend_class)
    strength = BULK_CLASS_TO_STRENGTH.get(label, 0)
    return strength, label or "n/a"


def dur_direction(rho, thr=0.30):
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


def fig2c_trend_state(af_aligned, ref_aligned):
    # af_aligned / ref_aligned are already signed by program direction.
    if af_aligned <= -FIG2C_EPS_REVERSED:
        return "Reversed"
    ref = max(float(ref_aligned), FIG2C_EPS_FLOOR)
    ratio = float(af_aligned) / ref
    if ratio < FIG2C_RATIO_ATTENUATED:
        return "Attenuated"
    if ratio <= FIG2C_RATIO_AMPLIFIED:
        return "Preserved"
    return "Amplified"
