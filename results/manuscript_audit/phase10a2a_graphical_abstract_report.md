# Phase 10A2-A — Graphical Abstract Final Freeze

**Final verdict:** `GRAPHICAL_ABSTRACT_FROZEN`

Snapshot UTC: 2026-07-17T08:21:02.778645+00:00

## 1. Builder SHA256

| Item | SHA256 |
|------|--------|
| Pre-freeze builder | `594c8072eeb3d8a1a216be0246c4742283d72163ca86a86e8e7d7434925c938a` |
| Freeze builder (post-edit) | `f2fc2b0432c4a85f827277899b36425477ee557b2dbff76706d93189adf705b8` |

Allowed edit: output routing for Phase 10A2-A audit paths; `WRITE_BUILDER_SIDE_AUDITS=False`; three boundary-phrase scrubs in disclaimer lines only (scientific numbers unchanged).

## 2. Previous asset hashes (pre-freeze)

| Asset | SHA256 | Bytes |
|-------|--------|-------|
| PDF | `03093e5f9d5a9579beec53b1d4efd480bbeb1acbe3ea4b1c55be3378c17333ba` | 54279 |
| SVG | `dd5c0c9ebe0b6a6f7b685f5760b2d7d81e6132ad3cf48bab68dc19814aac6289` | 77550 |
| PNG | `b257abc2cab4d6748e7ea61e97392ccf9c8229701e705e9332475591398d30e1` | 1561314 |

Source: `phase10a2a_before_graphical_abstract_state.tsv`

## 3. Final asset hashes

| Asset | SHA256 | Bytes |
|-------|--------|-------|
| PDF | `a33fb1294fdfee8aee724d8b187fc4fcfc2fb68233c55c1d2bb569f2269dd6e0` | 54570 |
| SVG | `bc36e5ae6f52e551c4d158e5d0344a135c76280e584a8a2e8f912dfa6177b362` | 77545 |
| PNG | `d9435af9de41b49382f82eeadef1824a6d53b94a001b319c63784d2ff527fddc` | 1560867 |
| Preview | `d9435af9de41b49382f82eeadef1824a6d53b94a001b319c63784d2ff527fddc` | 1560867 |

## 4. Canvas dimensions

- Specified: **240 × 160 mm**
- Observed after `fig.canvas.draw()`: **240.0 × 160.0 mm** (9.449 × 6.299 in)

### Panel rectangles (figure fraction)

- `ax`: x0=0.02, y0=0.925, w=0.96, h=0.0625
- `column_1`: x0=0.015, y0=0.15, w=0.3, h=0.7562
- `column_2`: x0=0.325, y0=0.15, w=0.405, h=0.7562
- `column_3`: x0=0.74, y0=0.15, w=0.245, h=0.7562
- `bulk_panel`: x0=0.3371, y0=0.6406, w=0.3808, h=0.1875
- `snrna_panel`: x0=0.3371, y0=0.375, w=0.3808, h=0.25
- `duration_panel`: x0=0.3371, y0=0.1625, w=0.3808, h=0.1969
- `path_a_panel`: x0=0.7521, y0=0.5531, w=0.2208, h=0.275
- `geomx_panel`: x0=0.7521, y0=0.3625, w=0.2208, h=0.175
- `genetics_panel`: x0=0.7521, y0=0.1625, w=0.2208, h=0.1844
- `bottom_summary`: x0=0.015, y0=0.0438, w=0.97, h=0.0875
- `footer`: x0=0.015, y0=0.0125, w=0.97, h=0.0219

## 5. Scientific assertion audit

Status: **PASS**

Required values verified in searchable SVG text (Healthy reference funnel, bulk High/Moderate/Low, snRNA F10/F11/F15/F16, duration Stable/Weak/Strong, Path A F10/F15, GeoMx 24 ROIs / 4 patients, genetics 7/12 / SMAD9 / PP3 / PP4, Non-nested).

Missing: none

## 6. Stale-term scan

Status: **PASS**

Forbidden-term hits: none

## 7. Overclaim scan

Status: **PASS**

Hits: none

## 8. Geometry QC

| Metric | Value | Required | Status |
|--------|-------|----------|--------|
| text_overlap_count | 0 | 0 | PASS |
| clipped_text_count | 0 | 0 | PASS |
| panel_boundary_violation_count | 0 | 0 | PASS |

## 9. Font QC

- minimum_font_size: **5.00 pt** (required ≥ 4.8) → **PASS**

## 10. Vector-text QC

- PDF vector_text: **YES** (TrueType/Type0=3, Type3=0)
- SVG vector_text / searchable `<text>`: **YES**

## 11. Embedded-raster QC

- SVG `<image>` tags: **0**
- SVG `data:image`: **NO**
- embedded_raster (SVG): **NO**

## 12. Panel-specific QC

| Panel | Status | Notes |
|-------|--------|-------|
| 1 Healthy reference | PASS | NCP expanded; funnel; not AF discovery |
| 2A Bulk | PASS | High 6 / Moderate 8 / Low 2 |
| 2B snRNA cells | PASS | four schematics; no cell-connecting arrows |
| 2C Duration | PASS | classification only; Stable 8 / Weak 6 / Strong 2 |
| 3 Path A + context | PASS | three gates; F10/F15; GeoMx/genetics separate |

## 13. Export settings

- `bbox_inches=None`
- no `tight_layout()`
- no `constrained_layout=True`
- PNG dpi = **600**
- white facecolor

## 14. Files touched (Phase 10A2-A)

- `scripts/graphical_abstract/build_graphical_abstract.py`
- `figures/graphical_abstract/graphical_abstract_final.pdf`
- `figures/graphical_abstract/graphical_abstract_final.svg`
- `figures/graphical_abstract/graphical_abstract_final.png`
- `results/manuscript_audit/phase10a2a_before_graphical_abstract_state.tsv`
- `results/manuscript_audit/phase10a2a_graphical_abstract_manifest.tsv`
- `results/manuscript_audit/phase10a2a_graphical_abstract_qc.tsv`
- `results/manuscript_audit/phase10a2a_graphical_abstract_report.md`
- `results/manuscript_audit/phase10a2a_graphical_abstract_preview.png`

`results/graphical_abstract_*` sidecars unchanged: **YES**

## 15. Final freeze verdict

`GRAPHICAL_ABSTRACT_FROZEN`

STOP. Do not begin license reconciliation, release commit, Zenodo, routing, or submission-pack migration.
