# Innovation Workflow (Version-Locked)

This document defines the **frozen implementation path** for the innovation
module and its deliverables. It is intended for reproducibility, internal QA,
and reviewer-facing traceability.

## Scope

The innovation module extends the baseline pipeline with:
- direction-aware consensus pathway scoring
- diversity- and quota-aware panel selection (v3 panel)
- perturbation robustness profiling
- v2-vs-v3 transferability benchmark
- v3 Figure 2 panel pack (b/c/d/e)

The pipeline promotes `NCP_v3` to the **primary manuscript panel** while
retaining `NCP_v2` as a **frozen baseline comparator** for method validation.

## Manuscript version lock (frozen)

- **Primary panel for main figures and core claims**: `NCP_v3`
- **Comparator panel (SI / benchmark only)**: `NCP_v2`
- **Comparator evidence files**:
  - `results/innovation_v2_vs_v3_metrics.tsv`
  - `results/innovation_v2_vs_v3_decision_table.tsv`
  - `results/innovation_v2_vs_v3_axis_effects.tsv`
  - `figures/innovation_v2_vs_v3_benchmark.pdf`

Interpretation lock:
- `NCP_v3` is the deployed panel.
- `NCP_v2` is preserved to quantify upgrade gain, not to compete as an
  alternative primary analysis path.

## Frozen Scripts

- `scripts/innovation_ncp_consensus.py`
- `scripts/innovation_ncp_robustness.py`
- `scripts/innovation_recompute_acai_panel.py`
- `scripts/innovation_v2_vs_v3_benchmark.py`
- `scripts/innovation_build_v3_axis_families.py`
- `scripts/innovation_fig2e_replication_lollipop_v3.py`
- `scripts/innovation_generate_fig2_v3_pack.py`
- `scripts/innovation_release_bundle.py`

## One-Command Build

Run:

`python3 scripts/innovation_release_bundle.py`

This executes the full pipeline in a fixed order and writes release metadata.

## Primary Deliverables

Core innovation tables:
- `results/innovation_ncp_consensus_long.tsv`
- `results/innovation_ncp_v3_panel.tsv`
- `results/innovation_ncp_v3_summary.tsv`
- `results/innovation_ncp_bootstrap_frequency.tsv`
- `results/innovation_ncp_bootstrap_summary.tsv`

Benchmark tables:
- `results/innovation_v2_vs_v3_metrics.tsv`
- `results/innovation_v2_vs_v3_axis_effects.tsv`
- `results/innovation_v2_vs_v3_decision_table.tsv`

v3 transfer tables:
- `results/innovation_ncp_v2_acai_long.tsv`
- `results/innovation_ncp_v2_acai_summary.tsv`
- `results/innovation_ncp_v3_acai_long.tsv`
- `results/innovation_ncp_v3_acai_summary.tsv`

v3 Figure 2 support tables:
- `results/innovation_ncp_v3_axis_families.tsv`
- `results/fig2c_chamber_decomposition_v3.tsv`
- `results/fig2d_family_forest_v3.tsv`
- `results/fig2e_replication_lollipop_v3.tsv`

Figures:
- `figures/innovation_ncp_bootstrap_frequency.pdf`
- `figures/innovation_v2_vs_v3_benchmark.pdf`
- `figures/fig2b_donor_acai_heatmap_v3.pdf`
- `figures/fig2c_chamber_decomposition_v3.pdf`
- `figures/fig2d_family_forest_v3.pdf`
- `figures/fig2e_replication_lollipop_v3.pdf`

Release integrity:
- `results/innovation_release_manifest.tsv`
- `results/innovation_release_summary.tsv`

## Operational Rule

- Do not manually edit generated TSV/PDF outputs.
- If any upstream script changes, rerun `innovation_release_bundle.py`.
- Use `innovation_release_manifest.tsv` checksums as the canonical release fingerprint.
- Keep both `NCP_v2` and `NCP_v3` outputs on disk for traceability; never delete
  `v2` benchmark artifacts when updating primary figures.
