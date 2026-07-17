# Phase 10A2-B — License reconciliation and release metadata preparation

**Final verdict:** `LICENSE_METADATA_SYNCHRONIZED`

## 1. Owner decision

| Field | Value |
|-------|-------|
| License | MIT |
| Copyright year | 2026 |
| Copyright holder | Xiao-Ce Dai |

No institutional ownership, additional copyright holders, author licenses, or contributor agreements inferred.

## 2. Files changed

| File | Action |
|------|--------|
| `LICENSE` | Updated copyright holder to `Xiao-Ce Dai` (MIT body retained) |
| `pyproject.toml` | `license = { text = "Proprietary" }` → `license = { text = "MIT" }` |
| `README.md` | Unchanged (no license statement present) |
| `CITATION.cff` | Not present; not created |
| `.zenodo.json` | Not present; not created |

## 3. License before / after

| File | Before | After |
|------|--------|-------|
| LICENSE | MIT; Copyright (c) 2026 Project P02 authors | MIT; Copyright (c) 2026 Xiao-Ce Dai |
| pyproject.toml | Proprietary | MIT |
| README.md | none explicit | none explicit (unchanged) |
| CITATION.cff | ABSENT | ABSENT (pending) |

## 4. Package identity

| Field | Expected | Observed | Status |
|-------|----------|----------|--------|
| distribution | acrf | acrf | PASS |
| import | crf | ['crf'] | PASS |
| version | 0.1.0 | 0.1.0 | PASS |

**Package identity verdict:** PASS

## 5. Metadata consistency

| Check | Status |
|-------|--------|
| LICENSE MIT + holder Xiao-Ce Dai | PASS |
| pyproject MIT | PASS |
| README (if stated) | N/A (silent) |
| CITATION (if present) | N/A (absent) |
| LICENSE ↔ pyproject conflict resolved | PASS |

**License consistency verdict:** PASS

## 6. Manuscript hash integrity

| Section | Status | SHA256 |
|---------|--------|--------|
| Abstract | PASS | `151e6a06613549716646bc578715a84a9df1e478f72dc21b1c4688cb0bbfe4bf` |
| Introduction | PASS | `c6f34a644f735f847b73a5e86be904c253309ba0979a4b5cb664e536c50e481c` |
| Results | PASS | `db2430c53382d1dd73eefc63d1f90658ccc06feda929e60e4609e0d1c68d6221` |
| Methods | PASS | `290926e278f916619d476c7902ae645b70c84be9bda0eb48c2a3d4a9c9e37ddc` |
| Discussion | PASS | `d73363f1c7f0e33a9b3aa909a9668b6d504dba8c063312fa3614c3ece8ec518d` |
| Figure_Captions | PASS | `b6521efcbff6fa6a0b4ed71eeab4b55b97bd24b769ba78dfacf5da8b02af7ed5` |
| Data_Code_Availability | PASS | `50a5b1c5f1badd55dbe4203934fd599baad6b680dfb01b2907b043bb846e5257` |

**Manuscript hash integrity verdict:** PASS

## 7. Scientific asset integrity

- Figure 3 PDF live hash match (phase10a0): PASS (`03c8e626d075fe9d6b65432ed5699c522ffb3b1b77bbdd430a7d7456e38af392`)
- Graphical abstract PDF frozen match (phase10a2a): PASS
- Claims registry: PASS (`results/handoff/fulltext_audit_v1_1/02_final_claims_registry.tsv`)
- Canon manifest live matches: 64; mismatches: 0 (stale handoff rows may pre-exist; not modified by this phase)

**Scientific asset integrity verdict:** PASS

Notes: No scientific assets modified in Phase 10A2-B.

## 8. Unresolved release items

- CITATION.cff absent (pending; not created)
- .zenodo.json absent (Zenodo out of scope)
- DCA license language still silent until later DCA finalize (correct)
- Release commit / tag / Zenodo DOI not created (out of scope)

## Audit artifacts

- `results/manuscript_audit/phase10a2b_before_license_audit.tsv`
- `results/manuscript_audit/phase10a2b_license_reconciliation.tsv`
- `results/manuscript_audit/phase10a2b_license_report.md`

## Final verdict

`LICENSE_METADATA_SYNCHRONIZED`

STOP. Do not begin environment freeze, release commit, Zenodo, combined manuscript, routing, or submission-pack migration.
