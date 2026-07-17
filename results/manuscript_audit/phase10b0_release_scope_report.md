# Phase 10B0 — Release scope preparation (read-only)

**Final verdict:** `RELEASE_SCOPE_READY_FOR_CLEAN_COMMIT`

Snapshot UTC: 2026-07-17T08:35:02.762972+00:00  
Git porcelain lines: **455** · Tracked files: **1551**

## 1. Release philosophy

Public release ships a **reproducible final-artifact state**: installable `acrf`/`crf` code, locked manuscript sources, final main/SI figures and tables, graphical abstract, provenance manifests, environment/hash verification metadata, and MIT license identity.

Public release does **not** redistribute raw external omics/GWAS blobs, obsolete narrative branches (Fig 7 / 5E / 6E / G-tier / PriorityScore / dual-axis), temporary previews, or stale submission-pack duplicates.

> The release environment documents reproduction of the final artifact state; it is not a claim of complete historical dependency reconstruction (Phase 10A3).

## 2. Included assets (`PUBLIC_RELEASE_REQUIRED` / optional)

### Required (summary)
| Group | Examples |
|-------|----------|
| Manuscript sources | Abstract/Author Summary, Introduction, Results, Methods, Discussion, Figure captions, Data/Code Availability, Author metadata |
| Generators | `scripts/m0`–`m6` manuscript builders; GA builder |
| Main figures | Figure 1–6 PDF/SVG/PNG (canonical paths) |
| SI figures | S1–S8 PDF/SVG/PNG under `figures/supplementary/` |
| Graphical abstract | `figures/graphical_abstract/graphical_abstract_final.{pdf,svg,png}` |
| SI tables | S1–S6 xlsx + machine-readable TSVs |
| Code | `src/crf/`, `scripts/`, `pyproject.toml`, tests (optional) |
| Metadata | `LICENSE`, `README.md` |
| Provenance | claims registry, GWAS provenance package, gene-set manifest, NCP frozen tables, Phase 10A3 env/hash manifests |

**Included asset rows (required+optional classifications):** 105

### Optional
- `.zenodo.json`
- `CITATION.cff`
- `results/figure1_composite_qc.md`
- `results/figure2_composite_qc.md`
- `results/figure6_composite_qc.md`
- `results/ncp_v3_singscore_reanalysis/figure3_composite_qc.md`
- `tests/`
- `tools/ldsc/`


## 3. Excluded assets

See `phase10b0_release_exclusion_manifest.tsv` (81 rows).

### Explicit large / raw exclusions
| Path | Reason | Reference |
|------|--------|-----------|
| `data/gwas/.../GCST90624411.h.tsv.gz` | large external sumstats (~1.60 GiB) | accession + `results/provenance/gwas_gcst90624411/` |
| `data/gwas/` | GWAS download tree | provenance package |
| `global_raw.h5ad` | raw expression object | public GEO via S1 |
| `Project_P2/`…`Project_P6/`, `bulk/` | cohort working trees / RAW.tar | processed results |
| `results/figure6_ldsc/*munge*.tsv.gz` | large intermediate | provenance + accession |

### Legacy / obsolete scientific narratives
Excluded when present on disk: Figure 5E, Figure 7, Figure 6E-like, PriorityScore/G-tier/dual-axis remnants, legacy SI ≥ S9, obsolete composites (`_composite_fig1`, integrated atlas), GA preview/debug sidecars under `results/graphical_abstract_*`.

### Development / stale packs
| Path | Status |
|------|--------|
| `submission/PLOS_CB_submission_pack/` | `SUBMISSION_ONLY` (stale; rebuild Phase 10D) |
| `results/handoff/fulltext_audit_v1/` | exclude superseded |
| `results/handoff/fulltext_audit_v1_1/` | archive selective (keep claims/manifests) |
| `Manuscript_outputs/` | exclude duplicates |

**Excluded asset count (union):** 82

**Large-file exclusion verdict:** PASS  
**Legacy asset exclusion verdict:** PASS

## 4. Optional assets

- Package `tests/`
- Figure QC markdown companions
- `tools/ldsc/` (vendored; not required for acrf import)
- `CITATION.cff` / `.zenodo.json` (create at tag/Zenodo; not blocking scope)

## 5. Owner decisions

**Owner decision count:** 0

No blocking owner decisions remain for scope classification. License MIT + holder already reconciled (Phase 10A2-B). GWAS sumstats exclusion follows explicit policy in this phase.

## 6. Estimated release size

| Bucket | Estimate |
|--------|----------|
| Manuscript assets (+ generators listed separately in code) | 434.2 KiB |
| Figures (main + SI + GA) | 19.9 MiB |
| Tables (SI) | 197.5 KiB |
| Code (package + scripts + tests) | 3.5 MiB |
| Provenance + QC + metadata | 268.3 KiB |
| **Total (sum of classified INCLUDE paths)** | **24.3 MiB** |
| Excluded large/legacy (catalogued) | 12.63 GiB |

Note: `scripts/` umbrella size dominates code; overlaps with individually listed generators are possible in the sum. No archive created.

## 7. Git readiness

| Metric | Value |
|--------|-------|
| Porcelain lines | 455 |
| Tracked files | 1551 |
| Commit-candidate rows | 105 |
| Git readiness verdict | **READY_WITH_FORCE_ADD_OR_GITIGNORE_ADJUST** |

### Git notes before a clean release commit
- Worktree is **DIRTY** (expected): many deleted legacy figure panels, modified LICENSE/pyproject/README, untracked GA/scripts trees.
- `docs/` and `manuscript_drafts` symlink are **gitignored** — locked manuscript sources must be force-added or `.gitignore` adjusted at commit phase (scope allows; not executed here).
- `figures/graphical_abstract/` and `scripts/graphical_abstract/` currently untracked — commit candidates.
- Do **not** force-add `data/gwas/**/*.h.tsv.gz`, `global_raw.h5ad`, or `Project_P*`.

Git path warnings (ignored/missing commit candidates):
- manuscript_drafts/Abstract_AuthorSummary_plos_cb.md: IGNORED; via manuscript_drafts/Abstract_AuthorSummary_plos_cb.md
- manuscript_drafts/Introduction.docx: IGNORED; via manuscript_drafts/Introduction.docx
- manuscript_drafts/Results.docx: IGNORED; via manuscript_drafts/Results.docx
- manuscript_drafts/Methods.docx: IGNORED; via manuscript_drafts/Methods.docx
- manuscript_drafts/Discussion.docx: IGNORED; via manuscript_drafts/Discussion.docx
- manuscript_drafts/Figure_Captions.docx: IGNORED; via manuscript_drafts/Figure_Captions.docx
- manuscript_drafts/Data_Code_Availability.docx: IGNORED; via manuscript_drafts/Data_Code_Availability.docx
- docs/manuscript/drafts/Author_Metadata.md: IGNORED


## 8. Unresolved issues before clean release commit

1. Adjust ignore/force-add policy for `docs/manuscript/drafts/*` locked sources.
2. Stage only INCLUDE commit candidates; drop deleted legacy panels intentionally.
3. Rebuild combined manuscript + submission pack in later phases (not required for code-release scope definition).
4. CITATION.cff / Zenodo JSON at tag time.
5. Optional: emit pip/R lockfile from Phase 10A3 inventory.
6. Propagate Figure 3 hashes into stale handoff manifests (Phase 10D/E; not scope-blocking).

## Final verdict

`RELEASE_SCOPE_READY_FOR_CLEAN_COMMIT`

STOP. Do not begin commit/tag creation.
