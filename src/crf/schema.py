"""Column names shared with the TSV exports. Renaming one breaks the figure scripts."""

# program-level panel, e.g. results/figure3_program_panel.tsv
COL_PROGRAM_ID = "program_id"
COL_PROGRAM_ORDER = "program_order"
COL_BASELINE_DIRECTION = "baseline_direction"  # +1 LA-biased, -1 RA-biased
COL_CONSENSUS_CLASS = "consensus_class"
COL_TREND_CLASS = "trend_class_GSE128188"

# donor x program long table
COL_DONOR_ID = "donor_id"
COL_RHYTHM = "rhythm"
COL_ACAI = "aCAI"
COL_COHORT = "cohort"
COL_FAMILY_ID = "family_id"
COL_DIRECTION = "direction"

# Fig 2C
COL_MEAN_ACAI_AF = "mean_aCAI_AF"
COL_MEAN_ACAI_SR = "mean_aCAI_SR"
COL_AF_OVER_REF = "af_over_ref"
COL_STATE = "state"

# Fig 5A
COL_RHO = "rho"
COL_AF_MAJOR_CARRIER_CT = "AF_major_carrier_CT"
COL_HEALTHY_TOP_EXPRESSOR_CT = "healthy_top_expressor_CT"
COL_ACAI_REL = "aCAI_rel"
