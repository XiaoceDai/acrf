"""
S2 — build donor x chamber x top_ct pseudobulk.

Filters (per locked pseudobulk spec + decision A):
  - source == "Nuclei"                     (Harvard-Nuclei + Sanger-Nuclei)
  - region in {LA, RA}
  - Used == "Yes"
  - cell_type in {Atrial_Cardiomyocyte, Fibroblast, Endothelial, Myeloid, Lymphoid}

Top cell type mapping:
  Atrial_Cardiomyocyte      -> atrial_cardiomyocyte
  Fibroblast                -> fibroblast
  Endothelial               -> endothelial_or_vascular
  Myeloid, Lymphoid         -> immune

Aggregation: sum raw counts per donor x region x top_ct.

Gating:
  - drop combos with n_cells < 30
  - cell type enters main analysis only if paired-donor count >= 4.

Outputs (one folder per cell type):
  pseudobulk/{ct}/counts.tsv    -> rows=gene symbols, cols=sample_id
  pseudobulk/{ct}/samples.tsv   -> sample metadata (sample_id, donor, chamber, n_cells)
  pseudobulk/summary.tsv        -> one row per (ct, donor, chamber), includes kept flag
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

H5AD = Path("/home/bio/Project_P02/global_raw.h5ad")
OUT = Path("/home/bio/Project_P02/pseudobulk")
OUT.mkdir(parents=True, exist_ok=True)

TOP_CT_MAP = {
    "Atrial_Cardiomyocyte": "atrial_cardiomyocyte",
    "Fibroblast": "fibroblast",
    "Endothelial": "endothelial_or_vascular",
    "Myeloid": "immune",
    "Lymphoid": "immune",
}
CT_SHORT = {
    "atrial_cardiomyocyte": "atrialCM",
    "fibroblast": "fibro",
    "endothelial_or_vascular": "endo",
    "immune": "immune",
}

N_CELLS_MIN = 30
PAIRED_DONOR_MIN = 4


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def main() -> None:
    log(f"opening {H5AD} in backed mode")
    a = ad.read_h5ad(H5AD, backed="r")
    log(f"  n_obs={a.n_obs:,}  n_vars={a.n_vars:,}")

    obs = a.obs
    mask = (
        (obs["source"].astype(str) == "Nuclei")
        & obs["region"].isin(["LA", "RA"])
        & (obs["Used"].astype(str) == "Yes")
        & obs["cell_type"].isin(list(TOP_CT_MAP.keys()))
    )
    idx = np.where(mask.values)[0]
    log(f"filter -> {idx.size:,} cells kept")

    # Load subset into memory (still sparse).
    log("loading filtered subset into memory (sparse)...")
    sub = a[idx].to_memory()
    a.file.close()

    X: sparse.csr_matrix = sub.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    log(f"  subset X: shape={X.shape} nnz={X.nnz:,} dtype={X.dtype}")

    meta = sub.obs[["donor", "region", "cell_type", "cell_source", "sample"]].copy()
    meta["top_ct"] = meta["cell_type"].astype(str).map(TOP_CT_MAP)
    meta["sample_id"] = (
        meta["top_ct"].astype(str)
        + "__"
        + meta["donor"].astype(str)
        + "_"
        + meta["region"].astype(str)
    )
    meta = meta.reset_index(drop=False).rename(columns={"index": "cell_id"})

    # Build indicator matrix M (n_samples, n_cells) with M[s,i]=1 if cell i -> sample s.
    sample_ids, inverse = np.unique(meta["sample_id"].values, return_inverse=True)
    n_samples = sample_ids.size
    log(f"distinct raw samples (before n_cells filter): {n_samples}")
    rows = inverse
    cols = np.arange(len(meta))
    data = np.ones_like(rows, dtype=np.int32)
    M = sparse.csr_matrix((data, (rows, cols)), shape=(n_samples, X.shape[0]))
    log(f"  indicator M: shape={M.shape} nnz={M.nnz:,}")

    # pseudobulk = M @ X  -> (n_samples, n_genes)  int32 counts
    log("computing pseudobulk = M @ X ...")
    pb = (M @ X).astype(np.int64)  # avoid overflow in huge CM pool
    log(f"  pseudobulk: shape={pb.shape} nnz={pb.nnz:,} dtype={pb.dtype}")

    # Collapse meta to one row per sample_id
    agg = (
        meta.groupby("sample_id", observed=True)
        .agg(
            donor=("donor", "first"),
            region=("region", "first"),
            top_ct=("top_ct", "first"),
            n_cells=("cell_id", "size"),
        )
        .reindex(sample_ids)
        .reset_index()
        .rename(columns={"index": "sample_id"})
    )
    agg["chamber"] = agg["region"]  # alias for SOP
    log(f"  aggregated samples: {len(agg)}")

    # Filter: n_cells >= 30
    kept = agg["n_cells"] >= N_CELLS_MIN
    log(f"n_cells>={N_CELLS_MIN}: kept {kept.sum()} / dropped {(~kept).sum()}")
    agg["pass_ncells"] = kept.values

    # Filter: paired-donor count per top_ct
    pair_stats = {}
    main_keep_ct: list[str] = []
    for ct, grp in agg[agg["pass_ncells"]].groupby("top_ct", observed=True):
        donors_la = set(grp.loc[grp["region"] == "LA", "donor"])
        donors_ra = set(grp.loc[grp["region"] == "RA", "donor"])
        paired = sorted(donors_la & donors_ra)
        pair_stats[ct] = paired
        log(f"  {ct}: paired_donors={len(paired)}  -> {'MAIN' if len(paired) >= PAIRED_DONOR_MIN else 'EXPLORATORY ONLY'}")
        if len(paired) >= PAIRED_DONOR_MIN:
            main_keep_ct.append(ct)

    # Gene symbols
    var = sub.var.copy()
    gene_ids = var.index.astype(str).values
    log(f"  gene space: {gene_ids.size:,}  (using var.index which is gene symbols)")

    # For each main cell type, write counts.tsv + samples.tsv
    # Only keep samples that (a) pass n_cells and (b) whose donor is in the paired set
    pb_csr = pb.tocsr()

    for ct in main_keep_ct:
        ct_dir = OUT / ct
        ct_dir.mkdir(parents=True, exist_ok=True)

        paired_donors = set(pair_stats[ct])
        # rows in `agg` we want
        mask_ct = (
            (agg["top_ct"] == ct)
            & agg["pass_ncells"]
            & agg["donor"].isin(paired_donors)
        )
        ct_meta = agg.loc[mask_ct, ["sample_id", "donor", "chamber", "n_cells"]].reset_index(drop=True)
        # Corresponding row indices in pb
        row_idx = np.array([np.where(sample_ids == sid)[0][0] for sid in ct_meta["sample_id"]])
        pb_ct = pb_csr[row_idx].toarray().astype(np.int64)
        log(f"  [{ct}] samples={len(ct_meta)}  pb shape={pb_ct.shape}")

        # counts.tsv: rows=genes, cols=sample_id
        counts_df = pd.DataFrame(
            pb_ct.T,
            index=pd.Index(gene_ids, name="gene"),
            columns=ct_meta["sample_id"].values,
        )
        # Drop genes that are 0 everywhere in this cell type (keeps file smaller;
        # downstream edgeR::filterByExpr will do a stricter pass anyway).
        nz_genes = (counts_df.sum(axis=1) > 0)
        counts_df = counts_df.loc[nz_genes]
        log(f"    counts nonzero-in-any-sample genes: {counts_df.shape[0]:,} / {gene_ids.size:,}")
        counts_df.to_csv(ct_dir / "counts.tsv", sep="\t")
        ct_meta.to_csv(ct_dir / "samples.tsv", sep="\t", index=False)
        log(f"    wrote {ct_dir}/counts.tsv  +  samples.tsv")

    # Global summary
    summary = agg.copy()
    summary["paired_donors_in_ct"] = summary["top_ct"].map(
        {ct: len(v) for ct, v in pair_stats.items()}
    ).fillna(0).astype(int)
    summary["main_analysis"] = summary["top_ct"].isin(main_keep_ct)
    summary.to_csv(OUT / "summary.tsv", sep="\t", index=False)
    log(f"wrote {OUT}/summary.tsv")

    # Pairing report
    pair_df = pd.DataFrame(
        [{"top_ct": ct, "n_paired": len(v), "paired_donors": ",".join(v)} for ct, v in pair_stats.items()]
    )
    pair_df.to_csv(OUT / "paired_donors.tsv", sep="\t", index=False)
    log(f"wrote {OUT}/paired_donors.tsv")

    log("done.")


if __name__ == "__main__":
    main()
