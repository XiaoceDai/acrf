"""
S1 — explore global_raw.h5ad schema without loading X.

Goal:
  * confirm obs columns that map to donor / region / cell_type / modality
  * confirm value distributions (LA/RA? 4 top cell types? multimodal?)
  * do NOT load the expression matrix
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

H5AD = Path("/home/bio/Project_P02/global_raw.h5ad")
REPORT = Path("/home/bio/Project_P02/reports/s1_schema.txt")
REPORT.parent.mkdir(parents=True, exist_ok=True)


def hline(title: str) -> str:
    return f"\n{'=' * 8} {title} {'=' * 8}"


def main() -> None:
    lines: list[str] = []
    w = lines.append

    w(f"file: {H5AD}  ({H5AD.stat().st_size / 1e9:.2f} GB)")

    # 1. Peek raw HDF5 top-level structure first (cheap, no decoding).
    w(hline("HDF5 top-level layout"))
    with h5py.File(H5AD, "r") as f:
        for k in f.keys():
            obj = f[k]
            kind = "group" if isinstance(obj, h5py.Group) else f"dataset {obj.shape} {obj.dtype}"
            w(f"  /{k:<12} {kind}")
        # encoding hints
        enc = dict(f.attrs) if len(f.attrs) else {}
        w(f"  root attrs: {enc}")

        # X info
        if "X" in f:
            x = f["X"]
            if isinstance(x, h5py.Group):  # sparse
                w(f"  /X (sparse): keys={list(x.keys())}, attrs={dict(x.attrs)}")
                for sk in x.keys():
                    d = x[sk]
                    w(f"    /X/{sk}: shape={d.shape} dtype={d.dtype}")
            else:
                w(f"  /X (dense): shape={x.shape} dtype={x.dtype}")

        # layers / raw
        for g in ("layers", "raw", "obsm", "varm", "uns"):
            if g in f:
                w(f"  /{g} keys: {list(f[g].keys())}")

    # 2. Open anndata in backed mode — does NOT load X.
    w(hline("anndata backed='r'"))
    a = ad.read_h5ad(H5AD, backed="r")
    w(f"  shape: n_obs={a.n_obs:,}  n_vars={a.n_vars:,}")
    w(f"  X file-backed: type={type(a.X).__name__}")

    # 3. obs schema
    w(hline("obs.columns (name | dtype | n_unique | sample values)"))
    for c in a.obs.columns:
        s = a.obs[c]
        try:
            nu = s.nunique(dropna=True)
        except Exception as e:
            nu = f"?({e})"
        sample = s.dropna().unique()[:6]
        w(f"  - {c:<30} {str(s.dtype):<15} nunique={nu}  sample={list(sample)}")

    # 4. var schema
    w(hline("var.columns"))
    for c in a.var.columns:
        s = a.var[c]
        w(f"  - {c:<30} {str(s.dtype):<15}")
    w(f"  var.index (first 10): {list(a.var.index[:10])}")
    w(f"  n_vars (genes/features): {a.n_vars:,}")

    # 5. obsm / varm / uns / layers (names only)
    w(hline("obsm / varm / uns / layers keys"))
    w(f"  obsm:   {list(a.obsm.keys())}")
    w(f"  varm:   {list(a.varm.keys())}")
    w(f"  layers: {list(a.layers.keys())}")
    try:
        w(f"  uns:    {list(a.uns.keys())}")
    except Exception as e:
        w(f"  uns:    <err {e}>")

    # 6. candidate metadata columns — value_counts for anything that looks
    # like donor / region / chamber / cell_type / modality / assay / sample
    w(hline("value_counts on candidate metadata columns"))
    candidates_regex = (
        "donor|sample|subject|patient|region|chamber|location|"
        "cell_type|celltype|cluster|annotation|ann_|leiden|louvain|label|"
        "modality|assay|tech|batch|dataset|study"
    )
    cand = [c for c in a.obs.columns if pd.Series([c]).str.contains(candidates_regex, case=False, regex=True).iloc[0]]
    w(f"  candidate columns: {cand}")
    for c in cand:
        s = a.obs[c]
        try:
            nu = s.nunique(dropna=False)
        except Exception:
            nu = -1
        w(f"\n  -- {c}  (nunique={nu}) --")
        if isinstance(nu, int) and 0 < nu <= 60:
            vc = s.value_counts(dropna=False)
            for k, v in vc.items():
                w(f"     {k!r:<40} {v:>10}")
        elif isinstance(nu, int):
            # too many distinct → show top 20 + size
            vc = s.value_counts(dropna=False).head(20)
            w(f"     (top 20 of {nu})")
            for k, v in vc.items():
                w(f"     {k!r:<40} {v:>10}")

    # 7. if there is a region/chamber col and a cell_type col, cross-tab
    w(hline("region × cell_type cross-tab (if identifiable)"))
    def pick(patterns: list[str]) -> str | None:
        for c in a.obs.columns:
            lc = c.lower()
            if any(p in lc for p in patterns):
                return c
        return None

    region_col = pick(["region", "chamber", "location", "tissue"])
    ct_col = pick(["cell_type", "celltype", "cell.type", "annotation"])
    donor_col = pick(["donor", "subject", "patient", "sample"])
    modality_col = pick(["modality", "assay", "tech"])
    w(f"  picked:  region={region_col!r}  cell_type={ct_col!r}  donor={donor_col!r}  modality={modality_col!r}")

    if region_col and ct_col:
        ct = pd.crosstab(a.obs[region_col], a.obs[ct_col], dropna=False)
        w(f"\n  crosstab region × cell_type (shape={ct.shape}):")
        w(ct.to_string())

    if donor_col and region_col:
        ct2 = pd.crosstab(a.obs[donor_col], a.obs[region_col], dropna=False)
        w(f"\n  crosstab donor × region (shape={ct2.shape}):")
        w(ct2.to_string())

    a.file.close()

    out = "\n".join(lines)
    REPORT.write_text(out)
    print(out)
    print(f"\n[wrote] {REPORT}")


if __name__ == "__main__":
    main()
