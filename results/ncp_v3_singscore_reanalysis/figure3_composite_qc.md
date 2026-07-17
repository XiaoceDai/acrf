# Final Figure 3 composite QC

Status: **PASS**

## Phase 6B0-R display correction
- scientific_data_changed = NO
- scientific_rerun = NO
- display_text_changed = YES
- change_reason = Figure 3 score-computability threshold clarification
- decision_code = PRODUCTION_3_ONLY_UPDATE_DISPLAY_NO_RERUN

## Physical geometry and vector assembly
- Final canvas: 178 × 225 mm.
- A: (0, 0), 178 × 48 mm; B: (0, 51.5), 178 × 86 mm.
- C: (0, 141), 87 × 84 mm; D: (91, 141), 87 × 84 mm.
- All panels nested at 100% physical scale with equal x/y transform; no stretching.
- PDF/SVG retain vector text, symbols, donor plots, and the B matrix; PNG is 600 dpi.
- Exactly one C/D y-axis title, one shared footer, and one badge per case block.

## Typography and containment
- Panel A: minimum visible font 6.0 pt; text bbox violations 0; numerical-annotation overlaps 0.
- Panel B: minimum visible font 6.0 pt; text bbox violations 0; numerical-annotation overlaps 0.
- Panel C: minimum visible font 6.0 pt; text bbox violations 0; numerical-annotation overlaps 0.
- Panel D: minimum visible font 6.0 pt; text bbox violations 0; numerical-annotation overlaps 0.
- Right-edge labels, class badges, module glyphs, and panel titles are contained.

## Scientific assertions
- A: 16-family axis; 34 donors; 9 compartments; 10,000 primary permutations.
- B: canonical 16-family order and 9-compartment order retained; carrier outlines and nominal-support dots retained.
- B class counts: Fidelity 1; Relay 2; Extension 1; Below-detection 12.
- F10 Fidelity: Fibroblast is the sole carrier and supported compartment.
- F11 Relay: Endothelial is carrier/not supported; Adipocyte is non-carrier/supported.
- F15 Relay: Fibroblast FALSE, Endothelial TRUE, Lymphocyte FALSE support; precision-sensitivity audit retained.
- F16 Extension: Fibroblast and aCM TRUE, Adipocyte FALSE support.
- C/D share y limits, ticks, point/mean/zero-line marks, block heights, separator, badge dimensions, and 6-pt annotation style.
- Renderers open frozen family-level sources only; no gene- or program-level result is shown.
- No visible `NCP_v3` text.

## Palette accessibility
- Semantic class palette matches Figure 3B in A–D.
- Minimum class RGB distance: 0.1811.
- Minimum class grayscale luminance distance: 0.0067.
- Heatmap grayscale sign distance: 0.1590.
- Heatmap minimum simulated CVD distance: 0.0929.

## SHA256
- `figures/ncp_v3_singscore_reanalysis/figure3_la_single_nucleus_localization.svg`: `293175c9e1c6726bc190391afa3d19d2ecbbf2b57db184589479de2858674efe`
- `figures/ncp_v3_singscore_reanalysis/figure3_la_single_nucleus_localization.pdf`: `03c8e626d075fe9d6b65432ed5699c522ffb3b1b77bbdd430a7d7456e38af392`
- `figures/ncp_v3_singscore_reanalysis/figure3_la_single_nucleus_localization.png`: `c23226e4d4dbb07c34411822c61bb972c55ccc913ab7a87f2484b5c4ab074547`
- `figures/ncp_v3_singscore_reanalysis/figure3_la_single_nucleus_localization_preview.png`: `45c04389071b2aef9d9143daa72bdaf65b6548af1001a445b7931452d091a413`
- `results/ncp_v3_singscore_reanalysis/figure3_composite_manifest.tsv`: `58325ce1e9c64060b167c2aa3ac146072f3da7fc6b042b498ba6c8134cd4244e`
- `results/ncp_v3_singscore_reanalysis/figure3cd_case_study_composite_geometry.tsv`: `a105e0f46f7a300e50dcd9be8da40e43dd2c1e7a9a5e5adbab8cd1e57992eed9`
- `results/ncp_v3_singscore_reanalysis/figure3cd_case_study_axis_contract.tsv`: `061566fa14fd27513c6731162300a26db32ddba5dbda35a7d1a5a89eddc83cc6`
