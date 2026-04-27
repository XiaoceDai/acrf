# Figure and pipeline scripts

## Layout

- **`scripts/v3/`** — NCP_v3 primary figure drivers (`*_v3*.py`), Fig 2 pack, v2→v3 benchmark, integrated atlas builder.
- **`scripts/`** (root) — legacy / cohort pipelines, pre-v3 figures, manuscript DOCX helpers, shared **`plot_style.py`** / **`_pathway_labels.py`**.

## v3 primary chain (sixteen families)

Primary v3 drivers target `figures/*_v3*.pdf` outputs under this repository layout.

| Area | Scripts under `scripts/v3/` (representative) |
|------|-----------------------------------------------|
| Fig 1 | (root) `fig1b_baseline_donor_aCAI.py`, `innovation_ncp_consensus.py`, `innovation_ncp_robustness.py`; **v3** `fig1e_v3_robustness_stability.py`, `fig1f_v3_HCA_baseline_attribution.py`, `fig1f_v3_HCA_figure.py` |
| Fig 2 | `innovation_build_v3_axis_families.py`, `innovation_generate_fig2_v3_pack.py`, `innovation_fig2e_replication_lollipop_v3.py` (+ root `fig2b_donor_acai_heatmap.py`, `fig2c_chamber_decomposition.py`, `fig2d_family_forest.py` with v3 TSV paths) |
| Fig 3 | `fig3a_v3_schematic.py` … `fig3f_v3_F4_BCAA_mini.py`, `fig3_v3_carrier_fidelity_class.py`, `fig3_v3_program_celltype_delta.py` |
| Fig 4 | `fig4a_v3_schematic.py` … `fig4def_v3_case_studies.py`, `fig4_v3_chronicity_class.py`, `fig4_v3_duration_trend.py` |
| Fig 5 | `fig_integrated_af_rewiring_atlas.py`, `fig5a_v3_schematic.py`, `fig5b_v3_triangulation_heatmap.py`, `fig5c_v3_cross_axis_consistency.py`, `fig5d_v3_geomx_roi.py`, `fig5e_v3_mechanism_cartoon.py` |
| Fig 6 | Locus / anchor tables under `results/figure6_mr_anchors_v3*.tsv`, `results/fig6a_v3_locus_zoom_summary.tsv`; R drivers (e.g. `h4_mr_coloc_anchors.R`) live under `scripts/`. |
| v2 → v3 benchmark | `innovation_v2_vs_v3_benchmark.py` |

Supplementary v3 matrices: `figS3_v3_*.py`, `figS4_v3_*.py`, `figS4b_v3_*.py`.

## Conventions

- **Repo root** for files in `scripts/v3/`: `Path(__file__).resolve().parents[2]`. For files in `scripts/`: `Path(__file__).resolve().parent.parent`.
- Add **`sys.path.insert(0, str(ROOT / "scripts"))`** when importing **`plot_style`** or **`_pathway_labels`** from `scripts/`.
- Prefer **`import crf`** for shared definitions (`pyproject.toml` distribution **`acrf`**).
- **Fig. 5B–5C v3** = `scripts/v3/fig5b_v3_triangulation_heatmap.py` and `scripts/v3/fig5c_v3_cross_axis_consistency.py`; legacy ten-program bridge remains `scripts/fig5a_triangulation.py`.
