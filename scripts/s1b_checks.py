"""
S1b — three quick checks that gate S2 design decisions.

1. Is X integer-valued counts? (edgeR needs integers)
2. Are modalities (cell_source/source) mixed inside each donor x region?
3. What does `Used`==Yes look like across region/cell_type?
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

H5AD = Path("/home/bio/Project_P02/global_raw.h5ad")
REPORT = Path("/home/bio/Project_P02/reports/s1b_checks.txt")
REPORT.parent.mkdir(parents=True, exist_ok=True)


def main() -> None:
    a = ad.read_h5ad(H5AD, backed="r")
    lines: list[str] = []
    w = lines.append

    # 1. X integer check — sample a few thousand rows only.
    w("======== X integer check ========")
    n = 5000
    idx = np.random.default_rng(0).choice(a.n_obs, size=n, replace=False)
    idx.sort()
    sub = a[idx].X  # this triggers I/O but small
    if hasattr(sub, "toarray"):
        sub_dense = sub.toarray()
    else:
        sub_dense = np.asarray(sub)
    nz = sub_dense[sub_dense != 0]
    w(f"  sampled rows: {n}")
    w(f"  dtype: {sub_dense.dtype}")
    w(f"  shape: {sub_dense.shape}")
    w(f"  n nonzero: {nz.size:,}")
    w(f"  min/max: {nz.min()} / {nz.max()}")
    w(f"  mean/std: {nz.mean():.4f} / {nz.std():.4f}")
    frac_int = float(np.mean(np.isclose(nz, np.round(nz))))
    w(f"  fraction integer-valued: {frac_int:.6f}")
    w(f"  first 20 nonzero: {nz[:20]}")

    # 2. modality mix inside donor x region
    w("\n======== modality per donor x region ========")
    obs = a.obs[["donor", "region", "cell_source", "source"]].copy()
    obs = obs[obs["region"].isin(["LA", "RA"])]
    w(f"  LA+RA cells total: {len(obs):,}")
    w(f"  cell_source levels: {sorted(obs['cell_source'].dropna().unique().tolist())}")
    w(f"  source levels:      {sorted(obs['source'].dropna().unique().tolist())}")

    by = obs.groupby(["donor", "region", "cell_source"], observed=True).size().unstack("cell_source", fill_value=0)
    w("\n  donor x region -> cell_source counts (LA/RA only):")
    w(by.to_string())

    # how many donors have >1 cell_source within a single region?
    mix = (by > 0).sum(axis=1)
    w(f"\n  distribution of #cell_source per donor-region (n={len(mix)}):")
    w(mix.value_counts().sort_index().to_string())

    # 3. Used flag
    w("\n======== Used flag ========")
    used = a.obs[["Used", "region", "cell_type"]].copy()
    w(f"  overall Used value_counts:\n{used['Used'].value_counts(dropna=False).to_string()}")
    w("\n  Used Yes within LA+RA, by cell_type:")
    sub_used = used[used["region"].isin(["LA", "RA"])]
    ct = pd.crosstab(sub_used["cell_type"], sub_used["Used"], dropna=False)
    w(ct.to_string())

    # 4. Map to 4 top cell types — how many survive `Used=Yes` + region in {LA,RA}?
    w("\n======== mapping to 4 top cell types (Used=Yes, LA|RA) ========")
    mapping = {
        "Atrial_Cardiomyocyte": "atrial_cardiomyocyte",
        "Fibroblast": "fibroblast",
        "Endothelial": "endothelial_or_vascular",
        "Myeloid": "immune",
        "Lymphoid": "immune",
    }
    obs_m = a.obs[["donor", "region", "cell_type", "Used", "cell_source"]].copy()
    obs_m = obs_m[obs_m["region"].isin(["LA", "RA"]) & (obs_m["Used"] == "Yes")]
    obs_m["top_ct"] = obs_m["cell_type"].map(mapping)
    obs_m = obs_m[obs_m["top_ct"].notna()]
    w(f"  after filter n={len(obs_m):,}")
    w("\n  cells per (donor x region x top_ct):")
    agg = obs_m.groupby(["top_ct", "donor", "region"], observed=True).size().unstack("region", fill_value=0)
    w(agg.to_string())

    w("\n  combos passing n_cells >= 30:")
    # flatten
    flat = obs_m.groupby(["top_ct", "donor", "region"], observed=True).size().reset_index(name="n")
    passed = flat[flat["n"] >= 30]
    w(f"  passed={len(passed)}  dropped={len(flat) - len(passed)}")
    # paired-donor check per top_ct
    w("\n  paired-donor count per top_ct (has both LA and RA with n>=30):")
    pivot = passed.pivot_table(index=["top_ct", "donor"], columns="region", values="n", fill_value=0)
    paired = (pivot > 0).all(axis=1)
    pair_count = paired.groupby(level="top_ct").sum()
    w(pair_count.to_string())

    out = "\n".join(lines)
    REPORT.write_text(out)
    print(out)
    print(f"\n[wrote] {REPORT}")


if __name__ == "__main__":
    main()
