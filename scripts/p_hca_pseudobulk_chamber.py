#!/usr/bin/env python3
"""
Build v3 chamber-resolved pseudobulk from the HCA Heart Cell Atlas
(global_raw.h5ad, Litvinukova et al. 2020 Nature).

Output one folder per (mapped) compartment under
    Project_P3/HCA/pseudobulk_chamber/<compartment>/
        expr.tsv      mean expression matrix (genes x donor_chamber sample_id)
        samples.tsv   sample metadata (sample_id, donor, region, n_cells, source, ...)

Aggregation:
    donor x region (LA / RA only) x cell_type
    weighted mean of raw counts -> sum(X) / n_cells (matches r2b_*.R style
    so downstream singscore / row-z behaves identically to the GSE255612
    pseudobulk_full pipeline).

Filters:
    - source == "Nuclei"            (Harvard-Nuclei + Sanger-Nuclei)
    - region in {LA, RA}
    - Used == "Yes"                 (HCA QC pass flag)
    - cell_type in HCA_TO_V3 mapping
    - n_cells per (compartment, donor, region) >= 30

Outputs
-------
    Project_P3/HCA/pseudobulk_chamber/<compartment>/expr.tsv
    Project_P3/HCA/pseudobulk_chamber/<compartment>/samples.tsv
    Project_P3/HCA/pseudobulk_chamber/summary.tsv
    Project_P3/HCA/pseudobulk_chamber/paired_donors.tsv
    Project_P3/HCA/pseudobulk_chamber/HCA_to_v3_compartment_map.tsv
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parent.parent
H5AD = ROOT / "global_raw.h5ad"
OUT = ROOT / "Project_P3" / "HCA" / "pseudobulk_chamber"

# HCA cell_type -> v3 compartment ontology (matches GSE255612 pseudobulk_full
# naming where possible; unmapped HCA labels are excluded).
HCA_TO_V3: dict[str, str] = {
    "Atrial_Cardiomyocyte": "atrialCM",
    "Fibroblast": "fibroblast",
    "Endothelial": "endothelial",
    "Pericytes": "pericyte",
    "Smooth_muscle_cells": "vascular_smooth_muscle",
    "Adipocytes": "adipocyte",
    "Mesothelial": "epicardial",
    "Neuronal": "neuronal",
    "Myeloid": "macrophage",
    "Lymphoid": "lymphocyte",
}
# Cell types intentionally dropped:
#   - Ventricular_Cardiomyocyte (not atrial)
#   - NotAssigned, doublets (QC flagged)

REGIONS = ("LA", "RA")
N_CELLS_MIN = 30
PAIRED_DONOR_MIN = 4


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main() -> None:
    if not H5AD.exists():
        sys.exit(f"missing input: {H5AD}")
    OUT.mkdir(parents=True, exist_ok=True)

    log(f"opening {H5AD} (backed='r')")
    a = ad.read_h5ad(H5AD, backed="r")
    log(f"  n_obs={a.n_obs:,}  n_vars={a.n_vars:,}")

    obs = a.obs
    used_col = obs["Used"].astype(str) if "Used" in obs.columns else None
    src_col = obs["source"].astype(str) if "source" in obs.columns else obs["cell_source"].astype(str)
    mask = (
        (src_col == "Nuclei")
        & obs["region"].isin(list(REGIONS))
        & obs["cell_type"].isin(list(HCA_TO_V3.keys()))
    )
    if used_col is not None:
        mask &= used_col == "Yes"
    keep_idx = np.where(mask.values)[0]
    log(f"after Nuclei+LA/RA+Used+ct filter: kept {keep_idx.size:,} cells")
    if keep_idx.size == 0:
        sys.exit("no cells passed filter; aborting")

    log("loading filtered subset into memory (sparse)")
    sub = a[keep_idx].to_memory()
    a.file.close()

    X = sub.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    else:
        X = X.tocsr()
    log(f"  X subset: shape={X.shape} nnz={X.nnz:,} dtype={X.dtype}")

    meta = sub.obs[["donor", "region", "cell_type", "source"]].copy()
    meta["compartment"] = meta["cell_type"].astype(str).map(HCA_TO_V3)
    meta = meta.dropna(subset=["compartment"])
    meta["sample_id"] = (
        meta["compartment"].astype(str)
        + "__"
        + meta["donor"].astype(str)
        + "_"
        + meta["region"].astype(str)
    )

    sample_ids, inv = np.unique(meta["sample_id"].values, return_inverse=True)
    n_samples = sample_ids.size
    log(f"raw (compartment,donor,region) groups: {n_samples}")

    rows = inv
    cols = np.arange(len(meta))
    data = np.ones_like(rows, dtype=np.int32)
    M_indicator = sparse.csr_matrix((data, (rows, cols)), shape=(n_samples, X.shape[0]))

    log("aggregating: pseudobulk = M @ X")
    pb_sum = (M_indicator @ X).astype(np.float64)
    n_cells = np.asarray(M_indicator.sum(axis=1)).ravel().astype(np.int64)
    pb_mean = pb_sum.multiply(1.0 / np.maximum(n_cells, 1)[:, None])
    pb_mean = np.asarray(pb_mean.todense())
    log(f"  pseudobulk mean: shape={pb_mean.shape}")

    agg = pd.DataFrame({"sample_id": sample_ids, "n_cells": n_cells})
    info = (
        meta.groupby("sample_id", observed=True)
        .agg(
            compartment=("compartment", "first"),
            donor=("donor", "first"),
            region=("region", "first"),
            source=("source", "first"),
        )
        .reset_index()
    )
    agg = agg.merge(info, on="sample_id", how="left")

    pass_n = agg["n_cells"] >= N_CELLS_MIN
    log(f"n_cells>={N_CELLS_MIN}: kept {pass_n.sum()} / dropped {(~pass_n).sum()}")
    agg["pass_ncells"] = pass_n.values

    pair_rows = []
    main_keep_ct: list[str] = []
    for comp, grp in agg[agg["pass_ncells"]].groupby("compartment", observed=True):
        donors_la = set(grp.loc[grp["region"] == "LA", "donor"])
        donors_ra = set(grp.loc[grp["region"] == "RA", "donor"])
        paired = sorted(donors_la & donors_ra)
        pair_rows.append({
            "compartment": comp,
            "n_paired": len(paired),
            "n_la_only": len(donors_la - donors_ra),
            "n_ra_only": len(donors_ra - donors_la),
            "paired_donors": ",".join(paired),
        })
        flag = "MAIN" if len(paired) >= PAIRED_DONOR_MIN else "EXPL_ONLY"
        log(f"  {comp:24s} paired_donors={len(paired):2d}  LA_only={len(donors_la - donors_ra):2d}  RA_only={len(donors_ra - donors_la):2d}  -> {flag}")
        if len(paired) >= PAIRED_DONOR_MIN:
            main_keep_ct.append(comp)

    pair_df = pd.DataFrame(pair_rows).sort_values("compartment")
    pair_df.to_csv(OUT / "paired_donors.tsv", sep="\t", index=False)
    log(f"wrote {OUT / 'paired_donors.tsv'}")

    var = sub.var
    gene_index = pd.Index(var.index.astype(str), name="gene")

    written = 0
    for comp in main_keep_ct:
        comp_dir = OUT / comp
        comp_dir.mkdir(parents=True, exist_ok=True)

        mask_ct = (agg["compartment"] == comp) & agg["pass_ncells"]
        ct_meta = agg.loc[mask_ct].reset_index(drop=True)
        order = np.array([np.where(sample_ids == sid)[0][0] for sid in ct_meta["sample_id"]])
        pb_ct = pb_mean[order, :]
        log(f"[{comp}] writing samples={len(ct_meta)}  matrix={pb_ct.shape}")

        expr = pd.DataFrame(pb_ct.T, index=gene_index, columns=ct_meta["sample_id"].values)
        nz = expr.sum(axis=1) > 0
        expr = expr.loc[nz]
        expr.to_csv(comp_dir / "expr.tsv", sep="\t", float_format="%.6g")
        ct_meta_out = ct_meta[["sample_id", "donor", "region", "n_cells", "source", "compartment"]].copy()
        ct_meta_out.to_csv(comp_dir / "samples.tsv", sep="\t", index=False)
        written += 1

    log(f"wrote {written} compartments")

    summary = agg.copy()
    summary["main_analysis"] = summary["compartment"].isin(main_keep_ct)
    summary.to_csv(OUT / "summary.tsv", sep="\t", index=False)

    map_df = pd.DataFrame(
        [{"hca_cell_type": k, "v3_compartment": v} for k, v in HCA_TO_V3.items()]
    )
    map_df.to_csv(OUT / "HCA_to_v3_compartment_map.tsv", sep="\t", index=False)
    log("done.")


if __name__ == "__main__":
    main()
