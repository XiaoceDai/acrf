# Figure 1 composite QC

Status: **PASS**

## Assembly geometry
- final size: 183.0 x 254.5 mm
- outer margins: 4.0 mm
- horizontal and vertical gutters: 4.5 mm
- columns: 68.2 / 102.3 mm (0.40 / 0.60)
- rows: 105.0 / 137.0 mm
- reading order: A -> B -> C -> D
- vector assembly: panel PDFs into composite PDF; self-contained panel SVGs into composite SVG
- raster export: 600 dpi
- intended manuscript width: 183.0 mm
- font family: DejaVu Sans
- minimum publication-facing text size: 6.5 pt
- PDF panel titles and forbidden visible-path text: checked
- SVG external file references: none

## Frozen-content assertions
- programs: 30
- program direction composition: 22 LA / 8 RA
- selection lineage: 80 -> 38 -> 25 strict + 5 fallback -> 30
- fixed topology: 30 programs -> 16 families
- Figure 1B and 1D canonical program order: identical
- canonical direction colors: LA #C97A74; RA #7D99B8
- panel source data or frozen identifiers modified by assembler: no

## Inputs
- panel_A: `figures/fig1a_donor_paired_estimand.pdf` — `13cf7997698383fb341823c7e9d347ee2cbafea117570c1d0020544a07b03af7`
- panel_B: `figures/ncp_v3_singscore_reanalysis/fig1b_healthy_reference_direction_magnitude.pdf` — `4ea85b9acc4a0acb892afda089d36eb4dae222851ca9ec7cb58d18f2829daccb`
- panel_C: `figures/fig1c_selection_provenance.pdf` — `2ea8a6fd87e8ccdcf345c5da4d70abd4d57d8a0803cafa8d27b6a6cd5ebdbd78`
- panel_D: `figures/fig1d_program_family_topology.pdf` — `41c97ed8f40dc164b5b1cbb7120307d438de0be40f66db813aeaf906ae383f28`
- frozen_source: `results/innovation_ncp_v3_panel.tsv` — `fc6da52a2bfe793ad836f49b6b3c5b38eaafd366408af135356de9d1d3e4a89c`
- frozen_source: `results/innovation_ncp_v3_axis_families.tsv` — `732921ee0be19ffd369f424961790748ce9bc7dfa49905355bf3926b1f51d3bb`
- frozen_source: `results/ncp_v3_canonical_display_order.tsv` — `7d065c73aa8dbe6ca7e8cf2847b2c31dbf1da53314c962489672e7f21cde5448`
- frozen_source: `results/ncp_v3_singscore_reanalysis/fig1b_canonical_display_order.tsv` — `c26b52bfa335c7add783e425bc11c2995ebba7c955fd01427c881f79a50965bb`

## Outputs
- `figures/figure1_composite.pdf` — `48f4be15a2b2268730df69c5936f089f7f0d6bb39f087a4134e87aa8dd24e827`
- `figures/figure1_composite.svg` — `cc07c00e7eafd9601d9f84a1c045642203e4f37b672bf99b7ec5a296a364d844`
- `figures/figure1_composite.png` — `03da64fcb4c81fe1f0fb8b921ba4b63aa827c1c8aa593a14ecc5f28fcfc3e559`
- `figures/figure1_composite_preview.png` — `af4fc6a5825701a184edbe31cfe6985edd674465ff439046bad1b682bdb63148`
