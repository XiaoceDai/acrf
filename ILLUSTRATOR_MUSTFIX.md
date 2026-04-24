# Illustrator must-fix list — final NC-level polish

After the matplotlib polish pass, the following items still benefit from
**Adobe Illustrator** (or Inkscape) editing of the generated PDFs.  This is
the standard last-mile workflow for Cell / Nature / Nat. Commun.; it is
**not** a re-analysis step.  Each item below is local (panel-level) and
should not affect underlying analysis or numbers.

---

## Tier 1 — Highly recommended (visible in first-look review)

### Figure 5

| Panel | Edit | Why |
|---|---|---|
| **5A** triangulation heatmap | Adjust column-header line breaks; remove vertical jitter on the program names; align the right "Triangulation class" pills to a single x-coordinate | matplotlib draws header text per axis; manual alignment looks sharper |
| **5C** localization inset | Increase scatter point size by ~20%; nudge the AA-position box / Spearman-ρ annotation away from the y-axis | matplotlib auto-padding crowds these |
| **5D** dual-axis mechanism | Replace the cell-type chips with **BioRender** glyphs (cardiomyocyte rod, fibroblast spindle, ECM mesh) and group + tighten the right-hand axis cards.  Currently the chips are typographic placeholders — readable but flat | BioRender is the standard for mechanism cartoons in this journal tier |

### Figure 6

| Panel | Edit | Why |
|---|---|---|
| **6A** axis forest | Tighten left-side cohort labels and the right-side `g [CI]` columns into a 2-column table inside the existing card; align the diamond (RE pooled) markers to one shared x-grid | makes the forest look like a published meta-analysis figure |
| **6B** program pooled | Rotate axis-membership side strips into a left-margin colour ribbon (instead of inline strips); align the per-row I² annotations | reads more like a Cell forest plot |
| **6C** anchor evidence matrix | Convert the per-cell "lead p = …" text + glyph stack into a single-row pill (icon + p-value side by side) for tighter layout | currently uses a 2-row vertical stack per cell |
| **6D** synthesis | Already typographic; replace the status pills (`ROBUST` / `ANCHOR` / `NOT SUPP.`) with subtle outlined chips + a small icon (✓ / ● / ✗) drawn in Illustrator if the journal style prefers icons | matplotlib's pills are deliberately flat |

### Graphical abstract

| Edit | Why |
|---|---|
| Replace the typographic "atrial silhouette" + cell mini-icons with **BioRender** atrial silhouette + cardiomyocyte / fibroblast pictograms | NC graphical abstracts are almost always BioRender or Illustrator |
| Tighten the dashed environment → mediator → axis arrow chain | matplotlib dashed arrows look uniform; manual control of dash lengths and arrow heads improves the upstream-modifier feel |

---

## Tier 2 — Nice to have

### Spec-side
- **Pull all of Figs 1–6 + S5–S11 into a single multi-page PDF master**
  (one figure per page, in publication order) for reviewer convenience.
- **Generate 300-dpi TIFF / PNG mirror exports** for journal upload (most
  NC submissions accept PDF directly, but a 300-dpi raster fallback is
  expected for some panels).

### Per-panel
- Fig 4C representative scatter: consider tightening `wspace` between the
  4 sub-panels by hand (matplotlib's `subplots_adjust` is uniform).
- Fig S8 falsification: the 4-panel S8a–d block is tall (9.6 inch); in
  Illustrator, consider a 2-column re-layout of S8a + S8d with S8b and
  S8c stacked at the right.
- Fig S10 / S11 methodology audit: convert the "target ≥ 0.7" annotations
  into right-aligned axis ticks for cleaner lay-out.

### Typography
- Pick **one** numeric font (e.g. JetBrains Mono / Inconsolata) for all
  monospace `g [CI]` and `p_emp` annotations across figures.
- Replace `× / ↑ / ↓` arrows with the journal's preferred Unicode glyphs
  (some journals require ASCII fallbacks).

---

## Tier 3 — Probably not worth it

- Hand-redrawing the 3-cohort centroid arrows in Fig 5C as Bezier paths.
- Adding journal-house-style colour gradients to the rewiring class
  pills in Fig 5A.

---

## Workflow note

Open each PDF in Illustrator with **"Convert text to outlines = OFF"**
to keep the typography editable.  All panels were exported with
`pdf.fonttype = 42` (TrueType embedded), so text is fully editable in
Illustrator without rasterisation.  Layers are flat — re-group them
into "Background", "Data", "Annotation" to make manual edits cleaner.

---

## Done in matplotlib already (not for Illustrator)

- Unified `set_nature_style()` rcParams across all panels (Arial /
  Helvetica, axis 6.4pt, ticks 5.6pt, captions 4.8pt, spine 0.6pt, no
  top/right spines).
- Uniform colour palette (LA blue / RA red, matrix orange, metabolic
  purple, AF brown / SR teal where applicable).
- Vector PDF output for every figure.
- `pad_inches = 0.04` consistent margin.
- Caveats / boundary statements in spec-locked verbatim form
  (`MANUSCRIPT_FRAMEWORK.md` §5).
