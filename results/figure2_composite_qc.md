# Figure 2 composite QC

Status: **PASS**

## Assembly geometry
- final size: 183.0 x 311.0 mm
- vertical margins: 4.0 mm
- horizontal gutters: 4.5 mm
- A/B width: 178.0 mm (centered with 2.5 mm side margins)
- C|D: 88.0 + 4.5 + 88.0 = 180.5 mm (centered with 1.25 mm side margins)
- row heights: A 55.0 / B 105.0 / C|D 134.0 mm
- reading order: A -> B -> C|D
- C and D use composite panel variants with shared row contract (row 5.60 mm; TABLE_TOP 27.0 mm; panel 134 mm)
- vector assembly: panel PDFs into composite PDF; self-contained panel SVGs into composite SVG
- raster export: 600 dpi
- intended manuscript width: 183.0 mm
- font family: DejaVu Sans
- PDF panel titles and forbidden visible-path text: checked
- SVG external file references: none

## Inputs
- panel_A: `figures/ncp_v3_singscore_reanalysis/fig2a_frozen_ncp_bulk_deployment.pdf` — `a426188abba343a525cf59e844c3d4a6e8466eef8423c151a9530a1ac96517c4`
- panel_B: `figures/fig2b_donor_acai_heatmap.pdf` — `0f9c118c6e9a6935f1860a89ac60802fc258a8ec1877116df58a75ef6d6d3afd`
- panel_C: `figures/ncp_v3_singscore_reanalysis/fig2c_family_cohort_state_map_composite.pdf` — `42ef7f73896918b2630d98dc88606ad6dd15f02138c780cbd97db9cdd770f4a0`
- panel_D: `figures/ncp_v3_singscore_reanalysis/fig2d_family_effects_consistency_composite.pdf` — `2e21f6745c4ec9e6a89e6c8a68d10f2af9909004d87bdcbde6b24027693af77b`
- composite_geometry: `results/ncp_v3_singscore_reanalysis/figure2_family_row_geometry_composite.tsv` — `73e203f76a57ce0f67e24c02b6dca63ad031787a5a7eb13bcc8b1ce491db937f`

## Outputs
- `figures/figure2_composite.pdf` — `42552821e90303fbffb08315d5a77683a5391c0c1876dd12b80f0a02b9de395d`
- `figures/figure2_composite.svg` — `3f6a1922273e74f92608d6d852f5bc11ec649bc1ce0a08ca397337241a13f9f2`
- `figures/figure2_composite.png` — `2b272d875604b346cca95f6c30383c6c4cbb0c5b0f8c60bcc90ace7406ea7682`
- `figures/figure2_composite_preview.png` — `d0afa8d40f92ec2b6a6296da239837cf9eaa71ab7cd68d4841abcd6e888b3a70`
