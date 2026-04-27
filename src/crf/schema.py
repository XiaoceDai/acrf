"""
Locked column names and table roles for TSV / dataframe interchange.

Use these constants when writing exporters so downstream tools (including
paper figures) stay byte-compatible with the repo contract.
"""

from __future__ import annotations

# --- program-level panel (e.g. results/figure3_program_panel.tsv) ---

COL_PROGRAM_ID = "program_id"
COL_PROGRAM_ORDER = "program_order"
COL_BASELINE_DIRECTION = "baseline_direction"  # +1 = LA-biased, -1 = RA-biased (NCP)
COL_CONSENSUS_CLASS = "consensus_class"
# Fig 2 overlay / fallback label (cohort-suffixed keys vary by export)
COL_TREND_CLASS = "trend_class_GSE128188"

# --- donor x program long table (e.g. fig2_primary_acai_long) ---

COL_DONOR_ID = "donor_id"
COL_RHYTHM = "rhythm"  # "AF" | "SR" etc.
COL_ACAI = "aCAI"  # aligned CAI (signed by NCP direction upstream)
COL_COHORT = "cohort"
COL_FAMILY_ID = "family_id"
COL_DIRECTION = "direction"  # "LA-biased" | "RA-biased"

# --- Fig 2C / trend layer ---

COL_MEAN_ACAI_AF = "mean_aCAI_AF"
COL_MEAN_ACAI_SR = "mean_aCAI_SR"
COL_AF_OVER_REF = "af_over_ref"
COL_STATE = "state"  # Preserved | Attenuated | Reversed | Amplified

# --- Fig 5A triangulation inputs ---

COL_RHO = "rho"
COL_AF_MAJOR_CARRIER_CT = "AF_major_carrier_CT"
COL_HEALTHY_TOP_EXPRESSOR_CT = "healthy_top_expressor_CT"

# --- aCAI_rel wide / long (see h5_falsification.R, h7_methodology_audit.R) ---

COL_ACAI_REL = "aCAI_rel"
