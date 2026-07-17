# Phase 10A3 — Environment freeze and canonical output hash verification

**Final verdict:** `HISTORICAL_ENVIRONMENT_UNRESOLVED_BUT_RELEASE_FREEZABLE`

Snapshot UTC: 2026-07-17T08:29:19.853701+00:00

> The release environment documents a reproducible execution environment for the
> final artifact state. It should not be interpreted as a complete reconstruction
> of every historical production dependency version.

## 1. Environment inventory

| Component | Version / status |
|-----------|------------------|
| OS | Linux 6.17.0-22-generic (Linux-6.17.0-22-generic-x86_64-with-glibc2.39) |
| Python | 3.12.3 |
| R | 4.5.1 |
| pip | 24.0 |
| renv.lock | ABSENT |
| root conda/environment.yml | ABSENT |
| Dockerfile | ABSENT |
| CI workflows | ABSENT |
| pyproject / acrf | 0.1.0 (MIT) |
| crf (src) | 0.1.0 |

See `phase10a3_environment_manifest.tsv`.

## 2. Current package versions

### Critical analysis packages
| Package | Version | Notes |
|---------|---------|-------|
| singscore | 1.28.1 | release record |
| matrixStats | 1.5.0 | release record |
| GSEABase | 1.70.1 | release record |
| edgeR | 4.6.3 | release record |
| limma | 3.64.3 | release record |
| coloc | 5.2.3 | **current/release only** |
| acrf | 0.1.0 | distribution |
| crf | 0.1.0 | import / src |

### Python scientific stack (current)
| Package | Version |
|---------|---------|
| pandas | 2.3.3 |
| numpy | 2.4.4 |
| scipy | 1.17.1 |
| matplotlib | 3.10.8 |
| scikit-learn | 1.8.0 |
| openpyxl | 3.1.5 |
| python-docx | 1.2.0 |

## 3. Historical environment limitations

| Item | Status |
|------|--------|
| Exact historical coloc version | **NOT RECOVERABLE** |
| Exact historical Python environment | **NOT RECOVERABLE** |

Do **not** claim exact historical reproducibility, identical historical package environment, or that current coloc was used historically.

## 4. Release environment proposal

Record the **current** runtime as the release reproduction environment:

- Python 3.12.3
- R 4.5.1
- coloc 5.2.3 (current only)
- singscore/matrixStats/GSEABase/edgeR/limma as tabulated above
- acrf/crf 0.1.0

A machine-readable pip/R lockfile was **not** written in this phase (allowed outputs restricted to audit files). Future release packaging may emit lockfiles from this inventory without scientific rerun.

## 5. Canonical hash verification summary

| Scope | PASS | FAIL/other | Total |
|-------|------|------------|-------|
| All verified assets | 78 | 0 | 78 |
| Release-required assets | 74 | 0 | 74 |

**Required rule:** all canonical assets `match = PASS` → SATISFIED

## 6. Figure 3 verification (Phase 6B0-R)

| Asset | Expected | Match |
|-------|----------|-------|
| PDF | `03c8e626…af392` | PASS |
| SVG | `293175c9…674efe` | PASS |
| PNG | `c23226e4…074547` | PASS |

**Figure 3 hash verdict:** PASS

Stale historical handoff/manifest rows were **not** modified; mismatches there are reported only in prior audits.

## 7. Graphical abstract verification (Phase 10A2-A)

| Asset | Expected | Match |
|-------|----------|-------|
| PDF | `a33fb129…dd6e0` | PASS |
| SVG | `bc36e5ae…17b362` | PASS |
| PNG | `d9435af9…27fddc` | PASS |

**Graphical abstract hash verdict:** PASS

## 8. Scientific rerun boundary

| Item | Status |
|------|--------|
| Scientific analyses | LOCKED |
| Scientific outputs | LOCKED |
| Figures | HASH VERIFIED |
| Tables | HASH VERIFIED |
| Historical environment | NOT RECOVERABLE |
| Release environment | FREEZABLE |
| New analysis | NOT REQUIRED |

**Scientific rerun required:** **NO**

## 9. Unresolved release items

- Historical coloc / Python pins: NOT RECOVERABLE
- Release pip/R lockfile file not emitted in this phase (audit-only)
- CITATION.cff / Zenodo DOI / release commit+tag: pending later phases
- DCA final identifiers: pending
- Combined manuscript / routing / submission pack: not started

## 10. Final verdict

`HISTORICAL_ENVIRONMENT_UNRESOLVED_BUT_RELEASE_FREEZABLE`

### Blocking issues
- None for canonical hash verification
- Historical coloc/Python pins remain NOT RECOVERABLE (non-blocking for release freezability)

STOP. Do not begin Data Availability final identifiers, release commit, release tag, Zenodo, combined manuscript, or routing/submission-pack migration.
