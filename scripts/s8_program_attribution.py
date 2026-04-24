#!/usr/bin/env python3
"""
S8 - Program-level attribution of bulk LA/RA programs to cell types.

Inputs
------
- normal_chamber_programs.tsv        (bulk: 10 stable programs)
- results/{ct}_chamber_fgsea.tsv     (sc fgsea for ct in {atrialCM, fibro, endo, immune})

Attribution logic
-----------------
For each bulk program P with direction D (LA-biased / RA-biased):
  For each cell type C, look up P's pathway id in C's sc fgsea and compute:
    - sc_NES, sc_padj, sc_direction (sc direction is "LA_up" / "RA_up")
    - leading_edge_overlap = |bulk.core_genes intersect sc.leading_edge| / |bulk.core_genes|
    - leading_edge_overlap_n = |bulk.core_genes intersect sc.leading_edge|
  C is a *carrier* of P if:
    - sc_direction matches D          (LA-biased <-> LA_up, RA-biased <-> RA_up)
    - sc_padj < PADJ_THRESH (0.25)
    - leading_edge_overlap_n >= OVERLAP_MIN_N (3)
    - leading_edge_overlap >= OVERLAP_MIN_JACCARD (0.10)  # overlap / |core_genes|

Labels (user spec)
------------------
  0 carriers  -> "unresolved"
  1 carrier   -> that cell type's display name
  >=2 carriers-> "mixed"

dominant_cell_type  = carrier with max |sc_NES| (or the label above)
support_cell_types  = other carriers (if mixed), comma-separated, or "none"
attribution_strength = |sc_NES| * leading_edge_overlap  of the strongest carrier
                       (0 if unresolved)

Output
------
results/normal_chamber_program_attribution.tsv
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

CTS = ["atrialCM", "fibro", "endo", "immune"]
CT_DISPLAY = {
    "atrialCM": "atrialCM",
    "fibro":    "fibroblast",
    "endo":     "endothelial_or_vascular",
    "immune":   "immune",
}

PADJ_THRESH        = 0.25
OVERLAP_MIN_N      = 3
OVERLAP_MIN_JACC   = 0.10   # intersection / |bulk.core_genes|


def strip_bulk_prefix(program_id: str) -> str:
    """`HALLMARK__HALLMARK_OXIDATIVE_PHOSPHORYLATION` -> `HALLMARK_OXIDATIVE_PHOSPHORYLATION`."""
    s = program_id.strip('"').strip()
    s = re.sub(r"^(HALLMARK|REACTOME)__", "", s)
    return s


def split_genes(cell: object) -> list[str]:
    if pd.isna(cell):
        return []
    s = str(cell).strip().strip('"')
    if not s or s.upper() == "NA":
        return []
    return [g.strip() for g in s.split(",") if g.strip()]


def bulk_dir_to_sc(d: str) -> str:
    d = str(d).strip().strip('"')
    if d.startswith("LA"):
        return "LA_up"
    if d.startswith("RA"):
        return "RA_up"
    return d


def main() -> None:
    bulk = pd.read_csv(ROOT / "normal_chamber_programs.tsv", sep="\t")
    bulk["program_id"]   = bulk["program_id"].astype(str).str.strip('"')
    bulk["program"]      = bulk["program"].astype(str).str.strip('"')
    bulk["direction"]    = bulk["direction"].astype(str).str.strip('"')
    bulk["pathway_key"]  = bulk["program_id"].map(strip_bulk_prefix)
    bulk["core_gene_list"] = bulk["core_genes"].map(split_genes)
    bulk["expected_sc_dir"] = bulk["direction"].map(bulk_dir_to_sc)

    sc = {}
    for ct in CTS:
        sc_df = pd.read_csv(ROOT / f"results/{ct}_chamber_fgsea.tsv", sep="\t")
        sc[ct] = sc_df.set_index("pathway")

    rows = []
    for _, brow in bulk.iterrows():
        pkey = brow["pathway_key"]
        bulk_dir_sc = brow["expected_sc_dir"]
        core_genes = brow["core_gene_list"]
        core_n = max(len(core_genes), 1)

        per_ct: dict[str, dict] = {}
        for ct in CTS:
            if pkey in sc[ct].index:
                scrow = sc[ct].loc[pkey]
                le_genes = split_genes(scrow.get("leading_edge_genes", ""))
                overlap = set(core_genes) & set(le_genes)
                overlap_n = len(overlap)
                jaccard = overlap_n / core_n
                per_ct[ct] = {
                    "nes":         float(scrow["NES"])  if pd.notna(scrow["NES"])  else float("nan"),
                    "padj":        float(scrow["padj"]) if pd.notna(scrow["padj"]) else float("nan"),
                    "direction":   str(scrow["direction"]),
                    "overlap_n":   overlap_n,
                    "overlap_jac": jaccard,
                }
            else:
                per_ct[ct] = {
                    "nes": float("nan"),
                    "padj": float("nan"),
                    "direction": "NA",
                    "overlap_n": 0,
                    "overlap_jac": 0.0,
                }

        carriers = []
        for ct in CTS:
            v = per_ct[ct]
            if (
                v["direction"] == bulk_dir_sc
                and pd.notna(v["padj"]) and v["padj"] < PADJ_THRESH
                and v["overlap_n"] >= OVERLAP_MIN_N
                and v["overlap_jac"] >= OVERLAP_MIN_JACC
            ):
                carriers.append(ct)

        if len(carriers) == 0:
            label = "unresolved"
            dominant = "unresolved"
            support = "none"
            strength = 0.0
        else:
            carriers_sorted = sorted(
                carriers, key=lambda c: abs(per_ct[c]["nes"]), reverse=True
            )
            best = carriers_sorted[0]
            strength = abs(per_ct[best]["nes"]) * per_ct[best]["overlap_jac"]
            if len(carriers) == 1:
                label = CT_DISPLAY[best]
                dominant = CT_DISPLAY[best]
                support = "none"
            else:
                label = "mixed"
                dominant = CT_DISPLAY[best]
                support = ",".join(CT_DISPLAY[c] for c in carriers_sorted[1:])

        # Compose row with rich context (for traceability)
        row = {
            "pathway":            brow["program"],
            "program_id":         brow["program_id"],
            "bulk_direction":     brow["direction"],
            "bulk_NES":           brow["NES"],
            "bulk_FDR":           brow["FDR"],
            "bulk_core_gene_n":   len(core_genes),
            "dominant_cell_type": dominant,
            "support_cell_types": support,
            "attribution_strength": round(strength, 4),
            "label":              label,
        }
        for ct in CTS:
            v = per_ct[ct]
            short = ct
            row[f"NES_{short}"]       = v["nes"]
            row[f"padj_{short}"]      = v["padj"]
            row[f"sc_dir_{short}"]    = v["direction"]
            row[f"LE_overlap_n_{short}"]   = v["overlap_n"]
            row[f"LE_overlap_jac_{short}"] = round(v["overlap_jac"], 4)
        rows.append(row)

    out = pd.DataFrame(rows)

    # Reorder: user-required core columns first, context after
    core_cols = [
        "pathway", "bulk_direction", "dominant_cell_type",
        "support_cell_types", "attribution_strength",
    ]
    other_cols = [c for c in out.columns if c not in core_cols]
    out = out[core_cols + other_cols]

    outpath = ROOT / "results/normal_chamber_program_attribution.tsv"
    out.to_csv(outpath, sep="\t", index=False, float_format="%.4g")
    print(f"[S8] Wrote {outpath.relative_to(ROOT)}  ({len(out)} rows)")

    print("\n[S8] Per-program attribution:")
    for _, r in out.iterrows():
        print(
            f"  {r['pathway']:<55s} | {r['bulk_direction']:<10s} | "
            f"dominant={r['dominant_cell_type']:<25s} "
            f"support={r['support_cell_types']:<30s} "
            f"strength={r['attribution_strength']:.3f}"
        )

    print("\n[S8] Label counts:")
    print(out.groupby(["bulk_direction", "label"]).size().to_string())


if __name__ == "__main__":
    main()
