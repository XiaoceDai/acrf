# Phase 11B0-BMC — Locked Limited Benchmark and Ablation Protocol

**Status:** PROTOCOL FREEZE (planning only — no results calculated)  
**Target use:** BMC Bioinformatics Research Article method evaluation  
**Working title context:** acrf chamber-aware reference framework  

Companion tables:

- `phase11b0_benchmark_input_lock.tsv`
- `phase11b0_benchmark_design.tsv`
- `phase11b0_metric_definitions.tsv`
- `phase11b0_comparator_estimand_map.tsv`
- `phase11b0_execution_plan.tsv`

---

## 1. Purpose

Evaluate **cross-cohort comparability**, **stability**, **contribution of acrf components**, and **robustness of the chamber-aware coordinate** on the locked paired-bulk layer.

Does **not** evaluate clinical prediction, diagnostic accuracy, causal correctness, therapeutic utility, or biological ground-truth accuracy.

**Forbidden metrics:** accuracy, AUROC, AUPRC, sensitivity, specificity (no prediction task).

---

## 2. Locked inputs

| Item | Lock |
|------|------|
| Cohorts | GSE128188, GSE79768, GSE115574 (3) |
| Complete paired donors | 51 (SR 26 / AF 25) |
| Programs / families | 30 / 16 frozen |
| Leading edges | Selection-time frozen LE (same for all methods) |
| Healthy reference | GSE112339; 8 healthy paired donors; frozen \(s_p\), \(M_p\) |
| Missing-gene rule | `MIN_MEASURED_GENES = 3` |
| Aggregation | Program-first equal-weight family mean |
| Layer | **Paired-bulk only** |

**Prohibited in all methods:** AF-outcome pathway reselection; family redefinition; unjustified donor-subset changes; tuning on benchmark results; use of snRNA, duration, GeoMx, or genetics.

Authorities and SHA256 hashes: `phase11b0_benchmark_input_lock.tsv`.  
**Do not use rendered figures as numerical authority.**

---

## 3. Reference method — A1_FULL_ACRF

Required components:

1. Donor-paired LA–RA contrast: \(\mathrm{CAI} = s_{\mathrm{LA}} - s_{\mathrm{RA}}\)  
2. Frozen healthy chamber direction \(s_p \in \{+1,-1\}\)  
3. Direction alignment: \(\mathrm{aCAI} = s_p \cdot \mathrm{CAI}\)  
4. Program-specific healthy-magnitude scaling: \(\mathrm{aCAI_{rel}} = \mathrm{aCAI}/|M_p|\)  
5. Frozen program definitions (exact singscore on frozen LE)  
6. Program-first family aggregation  
7. Identical cross-cohort synthesis: cohort Hedges’ \(g\) (AF vs SR) → FE-IVW → High/Moderate/Low  

Consistency class rule (locked):

- **High:** pooled CI excludes 0 **and** \(n_{\mathrm{sign\_match}} \ge 2\)  
- **Moderate:** \(n_{\mathrm{sign\_match}} \ge 2\) (pooled not supported)  
- **Low:** otherwise  

A1 is the **reference configuration**, not ground truth. Locked `fig2d_family_effects_consistency_source.tsv` may be reused if bit-reproduced.

---

## 4. Design ablations (exactly four)

### A2_NO_DONOR_PAIRING

**Question:** What is lost when chamber contrasts ignore within-donor pairing?

**Frozen fair estimand (STOP condition cleared):**

1. Use the same complete-pair donors and LA/RA samples.  
2. Within each rhythm \(r \in \{\mathrm{AF},\mathrm{SR}\}\), define the **unpaired chamber mean contrast**  
   \[
   C_r = \overline{s}_{\mathrm{LA},r} - \overline{s}_{\mathrm{RA},r}
   \]
   treating LA and RA as independent samples of size \(n_r\) (same \(n_r\) as donors).  
3. Align and scale: \(aC_{\mathrm{rel},r} = s_p \cdot C_r / |M_p|\).  
4. Program-first family aggregation on \(aC_{\mathrm{rel},r}\).  
5. AF−SR effect = \(aC_{\mathrm{rel},\mathrm{AF}} - aC_{\mathrm{rel},\mathrm{SR}}\).  
6. **Uncertainty (frozen before execution):** Welch SE for each unpaired chamber contrast ignoring within-donor covariance; propagate to AF−SR difference under independence.  

**Not allowed:** adding samples; treating technical replicates as donors; using pairing-permutation stress as a substitute for A2 (permutation is optional robustness R5).

### A3_NO_DIRECTION_ALIGNMENT

Retain pairing and \(M_p\) scaling; omit only \(s_p\).

\[
\mathrm{CAI_{rel}} = \frac{\mathrm{CAI}}{|M_p|}
\]

Program-first family mean of **unaligned** \(\mathrm{CAI_{rel}}\). Opposite-direction programs enter family means in native LA−RA orientation (possible cancellation) — must be documented in outputs.

### A4_NO_HEALTHY_MAGNITUDE_SCALING

Retain pairing and \(s_p\); use unscaled \(\mathrm{aCAI}\) (no `/|M_p|\); no alternative post hoc normalization.

### A5_CHAMBER_AGNOSTIC

**Recommended fair design (frozen):**

1. Score LA and RA with the same frozen singscore.  
2. Within each donor: \(\mathrm{donor\_score} = \mathrm{mean}(s_{\mathrm{LA}}, s_{\mathrm{RA}})\).  
3. Program-first family mean of donor scores.  
4. Hedges’ \(g\) AF vs SR on donor-level family scores.  

**Forbidden:** pooling LA and RA samples as independent observations (pseudoreplication STOP).

---

## 5. Scoring comparator — B1_PLAIN_SINGSCORE

**Purpose:** Compare acrf chamber calibration with plain frozen-gene-set scoring.

**Frozen design:**

- Same frozen LE, samples, cohorts, exact singscore.  
- **No** acrf direction alignment or healthy-magnitude scaling.  
- **Chamber-stratified estimands:** `B1_LA` and `B1_RA` separately (donor-level family scores within chamber → Hedges’ \(g\) AF vs SR).  

**Do not** claim singscore was designed to solve chamber calibration.  
**Do not** average chambers for B1 (that would collapse into A5).

### B1 vs A5 redundancy verdict

**NOT REDUNDANT.**

| | A5 | B1 |
|--|----|----|
| Chamber handling | Within-donor average | Stratified LA and RA arms |
| Calibration | None | None |
| Primary outputs | One effect vector | Two effect vectors (`B1_LA`, `B1_RA`) |

STOP condition for B1/A5 redundancy: **cleared**.

---

## 6. Primary metrics

### M1_CROSS_COHORT_EFFECT_CORRELATION

For each method/arm: Spearman correlation of the 16 family AF−SR effects for each of 3 cohort pairs.  
Report: three pairwise correlations, **median**, **range**.  
No P-value as primary summary.

### M2_CROSS_COHORT_SIGN_CONCORDANCE

Concordant family×cohort-pair sign matches / **48** (16×3).  
Report count, proportion, family-level breakdown.  
Sign convention: \(\mathrm{sign}(0)=0\); zero does not concord with nonzero.  
Not biological truth.

---

## 7. Secondary metrics

| ID | Summary |
|----|---------|
| M3 | Family-rank Spearman; median across pairs |
| M4 | Cochran Q / I² descriptive; n=3 caveat; not sole endpoint |
| M5 | High/Mod/Low retention vs A1; A1 not ground truth |
| M6 | Pooled-sign agreement with A1 (descriptive) |
| M7 | Estimand completeness / missingness |
| M8 | Runtime: 1 warm-up + ≥3 timed runs; median/range |
| M9 | Peak memory only if reproducible; else `NOT_MEASURED` |

---

## 8. Optional robustness (not in minimum)

| ID | Classification |
|----|----------------|
| R1 LOCO | REQUIRES_LIMITED_NEW_RUN or reaggregation; DEFER_TO_ADDITIONAL_FILE; NOT_REQUIRED |
| R2 Healthy-reference donor bootstrap | REQUIRES_LIMITED_NEW_RUN (selection bootstrap ≠ this); DEFER; NOT_REQUIRED |
| R3 Reference subsampling | REQUIRES_LIMITED_NEW_RUN; DEFER; NOT_REQUIRED |
| R4 Missing-gene sensitivity | DEFER / NOT_REQUIRED (legacy unverified) |
| R5 Pairing permutation stress | Distinct from A2; lineage-verify figS15 or re-run; DEFER; NOT_REQUIRED |

Do not expand the minimum benchmark without justification.

---

## 9. Statistical reporting

Benchmark is **primarily descriptive**. No multiple-hypothesis starring across comparator metrics.  
If resampling is added later: define unit, preserve donor structure, same design across comparable methods, **freeze seeds before execution**.

---

## 10. Planned output tables (do not generate now)

1. `benchmark_method_definitions.tsv`  
2. `benchmark_family_cohort_effects.tsv`  
3. `benchmark_pairwise_correlations.tsv`  
4. `benchmark_sign_concordance.tsv`  
5. `benchmark_rank_stability.tsv`  
6. `benchmark_heterogeneity.tsv`  
7. `benchmark_class_stability.tsv`  
8. `benchmark_runtime.tsv`  
9. `benchmark_summary.tsv`  

**Required columns on future rows:**  
`method_id`, `cohort`, `family`, `estimand`, `source_input_hash`, `software_version`, `seed` (where applicable).

---

## 11. Proposed BMC benchmark figure (not generated)

| Panel | Content |
|-------|---------|
| A | acrf and ablation definitions |
| B | Cross-cohort effect-correlation comparison (M1) |
| C | Sign concordance (M2) and rank stability (M3) |
| D | Consistency-class retention/reclassification (M5) |

Runtime → Additional file.

---

## 12. Claim boundaries

**Permitted if supported:** increased cross-cohort concordance; pairing/alignment/scaling contributed to stability/comparability; chamber-agnostic scoring differed or was less stable.

**Prohibited:** recovered true biology; universally superior; more accurate without ground truth; predicts clinical outcomes; proves causal remodeling.

---

## 13. STOP conditions (protocol gate)

| Condition | Status |
|-----------|--------|
| Ablation changes frozen program set | Forbidden by lock |
| Different donor eligibility without justification | Forbidden |
| A5 donor pseudoreplication | Forbidden; donor-mean design frozen |
| B1 indistinguishable from A5 | **Cleared — not redundant** |
| No-pairing estimator not definable | **Cleared — unpaired mean contrast + Welch SE frozen** |
| Comparator not mappable to fair estimand | Cleared via estimand map |
| Primary metrics selected after viewing results | Forbidden — M1/M2 predeclared |

---

## 14. Execution order

See `phase11b0_execution_plan.tsv` (S0–S12).  
**This phase ends at S0.** No benchmark computation.

---

## 15. Verdict

**BMC_BENCHMARK_PROTOCOL_LOCKED**
