# Figure and pipeline scripts

## Layout

- **`scripts/v3/`** — NCP primary figure drivers (`*_v3*.py`), Fig 2 pack, prototype-vs-NCP benchmark, integrated atlas builder.
- **`scripts/`** (root) — legacy / cohort pipelines, pre-v3 figures, manuscript DOCX helpers, shared **`plot_style.py`** / **`_pathway_labels.py`**.

**Supporting tables (submission):** `p1_demographics_table.py` → `results/cohort_demographics_table.tsv`; then `m_supp_table_s1_cohort_demographics.py` → `manuscript_drafts/Supplementary_Table_S1_cohort_demographics.docx` (export PDF in Word if required). **S6 anchor matrix:** `m10_supp_table_s6_anchor_matrix.py` → `manuscript_drafts/Supplementary_Table_S6_anchor_matrix.docx` from `results/figure6_mr_anchors_v3_pipeline_wide.tsv`. **Fig 6 prose-facing MR slice:** run `build_figure6_mr_v3_manuscript_tables.py` → `results/figure6_mr_anchors_v3_consolidated.tsv` + gate-rules / calibration / locus-window helpers.

## v3 primary chain (sixteen families)

Primary v3 drivers target `figures/*_v3*.pdf` outputs under this repository layout.

| Area | Scripts under `scripts/v3/` (representative) |
|------|-----------------------------------------------|
| Fig 1 | (root) `fig1b_baseline_donor_aCAI.py` (TSV only), `innovation_ncp_consensus.py`, `innovation_ncp_robustness.py`; **v3** `fig1e_v3_robustness_stability.py`, `fig1f_v3_HCA_baseline_attribution.py`, `fig1f_v3_HCA_figure.py` |
| Fig 2 | `innovation_build_v3_axis_families.py`, `innovation_generate_fig2_v3_pack.py`, `export_fig2b_donor_acairel_tables.py`, `innovation_fig2e_replication_lollipop_v3.py` (+ root `fig2b_donor_acai_heatmap.py`, `fig2c_chamber_decomposition.py`, `fig2d_family_forest.py` with v3 TSV paths) |
| Fig 3 | `fig3a_v3_schematic.py` … `fig3f_v3_F4_BCAA_mini.py`, `fig3_v3_carrier_fidelity_class.py`, `fig3_v3_program_celltype_delta.py` |
| Fig 4 | `fig4a_v3_schematic.py` … `fig4def_v3_case_studies.py`, `fig4_v3_chronicity_class.py`, `fig4_v3_duration_trend.py` |
| Fig 5 | `fig_integrated_af_rewiring_atlas.py`, `fig5a_v3_schematic.py`, `fig5b_v3_triangulation_heatmap.py`, `fig5c_v3_cross_axis_consistency.py`, `fig5d_v3_geomx_roi.py`, `fig5e_v3_mechanism_cartoon.py` |
| Fig 6 | Genetic traceability / coloc calibration. **Phase 1 (frozen):** `fig6a_anchor_coloc_calibration.py` → `figures/fig6a_anchor_coloc_calibration.pdf`; `fig6b_pp3_pp4_leaddistance.py` → `figures/fig6b_pp3_pp4_leaddistance.pdf`. Lead-variant offsets: `compute_anchor_lead_distance.py` → `results/fig6a_v3_locus_zoom_summary.tsv` (from `Project_P6/cache/` tabix slices). **Phase 2 (family overlap + G-tier):** `phase2_connectivity_probe.py` → `results/phase2_connectivity_probe.tsv`; `fetch_external_af_gene_sets.py` → `results/external_af_gene_sets.tsv`; `fig6c_family_af_overlap.py` → `figures/fig6c_family_af_overlap.pdf`; `build_figure6_g_tier_table.py` + `fig6d_traceability_matrix.py` → `figures/fig6d_traceability_matrix.pdf`. **Tables:** `build_figure6_mr_v3_manuscript_tables.py` → consolidated / gate-rules / calibration / locus-window helpers; G-tier → `results/figure6_g_tier_traceability.tsv`. **Frozen sources (do not overwrite without approval):** `results/figure6_mr_anchors_v3_pipeline_wide.tsv` + `results/fig6a_v3_locus_zoom_summary.tsv`. **R drivers (network tabix):** `h4_mr_coloc_anchors.R`, `h4c_mr_coloc_positive_controls.R`, `h4d_mr_coloc_v3_all_anchors.R` (diff-only vs frozen via `compare_coloc_frozen_vs_rerun.py`). **Universe prep:** `build_ncp_v3_family_gene_sets.py` → `results/ncp_v3_gene_universe.tsv` (N=14,635) + family leading-edge genes. |
| v2 → v3 benchmark | `innovation_v2_vs_v3_benchmark.py` |

Supplementary v3 matrices: `figS3_v3_*.py`, `figS4_v3_*.py`, `figS4b_v3_*.py`.

> **Fig 6 naming / method note.** Filenames and columns containing **`mr`**
> (e.g. `h4_mr_coloc_anchors.R`, `figure6_mr_*`) implement **cis-eQTL × AF-GWAS
> colocalisation (`coloc.abf`) anchor diagnostics** and **do not perform
> Mendelian randomization**. The AF GWAS is **Yuan et al. 2025 cross-population**
> (GWAS Catalog `GCST90624411`); Roselli 2025, Selewa 2023 and Choi 2025 are
> external risk-gene / regulatory references for the family-overlap layer, not
> the coloc GWAS. Anchor coloc numbers and lead-variant offsets are **frozen
> analysis artifacts**; Fig 6A/6B regenerate from them offline. Phase 2 adds
> family×external-set overlap (6C) and integrated G-tier readout (6D); external
> AF gene sets are reference overlays, not replacements for the Yuan coloc GWAS.
> End-to-end coloc recomputation uses `h4d_mr_coloc_v3_all_anchors.R` with
> diff-only reporting; frozen `pipeline_wide` is not overwritten unless approved.
> **Cursor / remote SSH preview:** each Fig 6 driver also writes a `.png` beside
> the `.pdf` (and copies to `figures/_audit_png/`) for in-IDE viewing when PDF
> preview fails.

## Conventions

- **Repo root** for files in `scripts/v3/`: `Path(__file__).resolve().parents[2]`. For files in `scripts/`: `Path(__file__).resolve().parent.parent`.
- Add **`sys.path.insert(0, str(ROOT / "scripts"))`** when importing **`plot_style`** or **`_pathway_labels`** from `scripts/`.
- Prefer **`import crf`** for shared definitions (`pyproject.toml` distribution **`acrf`**).
- **Fig. 5B–5C v3** = `scripts/v3/fig5b_v3_triangulation_heatmap.py` and `scripts/v3/fig5c_v3_cross_axis_consistency.py`; legacy ten-program bridge remains `scripts/fig5a_triangulation.py`.
