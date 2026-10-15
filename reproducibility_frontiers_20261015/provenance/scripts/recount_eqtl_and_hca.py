#!/usr/bin/env python3
"""Recount cached cis-eQTL rows and record the local HCA matrix identity.

Does not download data, does not rewrite Project_P6/cache, and does not
recompute family scores.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np

ROOT = Path("/home/bio/Project_P02")
OUT = ROOT / "revision/acrf_provenance_check_20261015/outputs"
GIT = ROOT / "revision/acrf_provenance_check_20261015/git_v020/results"
CACHE = ROOT / "Project_P6/cache"

ANCHORS = [
    ("PRELP", "ENSG00000188783", 1, 203475828, 203491352),
    ("OMD", "ENSG00000127083", 9, 92414245, 92424461),
    ("IDH2", "ENSG00000182054", 15, 90083045, 90102504),
    ("NDUFA8", "ENSG00000119421", 9, 122144058, 122159819),
    ("NTRK2", "ENSG00000148053", 9, 84668522, 85029699),
    ("NTF3", "ENSG00000185652", 12, 5404888, 5495026),
    ("FRS2", "ENSG00000166225", 12, 69464533, 69513650),
    ("BMP10", "ENSG00000163217", 2, 68869054, 68886842),
    ("SMAD9", "ENSG00000120693", 13, 36879326, 36925087),
    ("SMAD6", "ENSG00000137834", 15, 66703041, 66782082),
    ("PITX2", "ENSG00000164093", 4, 110610001, 110633149),
    ("KCNN3", "ENSG00000143603", 1, 154837802, 154961963),
]
# column positions in the headerless eQTL Catalogue slice used by h4*.R
C_TRAIT, C_POS, C_REF, C_ALT, C_VARIANT = 0, 2, 3, 4, 5
C_MAF, C_P, C_BETA, C_SE, C_GENE, C_RSID = 7, 8, 9, 10, 16, 18


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fnum(text: str):
    try:
        return float(text)
    except ValueError:
        return None


def load_stored() -> dict:
    cons = list(csv.DictReader((GIT / "figure6_mr_anchors_v3_consolidated.tsv").open(), delimiter="\t"))
    pos = list(csv.DictReader((GIT / "figure6_mr_positive_controls.tsv").open(), delimiter="\t"))
    out = {}
    for row in cons:
        out[row["anchor"]] = row
    for row in pos:
        if row["tissue"] == "atrial_appendage" and row["anchor"] not in out:
            out[row["anchor"]] = row
        elif row["anchor"] == "KCNN3" and row["tissue"] == "atrial_appendage":
            out["KCNN3"] = row
    return out


def recount() -> None:
    stored = load_stored()
    rows = []
    for symbol, gene, chrom, start, end in ANCHORS:
        lo, hi = start - 1_000_000, end + 1_000_000
        path = CACHE / f"eqtl_{symbol}_chr{chrom}_{lo}_{hi}.tsv.gz"
        rec = {
            "anchor": symbol,
            "gene_id": gene,
            "cache_file": str(path.relative_to(ROOT)) if path.exists() else "",
            "cache_bytes": path.stat().st_size if path.exists() else "",
            "cache_sha256": sha256(path) if path.exists() else "",
            "status": "cache_missing" if not path.exists() else "recounted_from_cache",
        }
        if not path.exists():
            rows.append(rec)
            continue
        raw = gene_n = filt_n = 0
        variants, poskeys, rsids = set(), set(), set()
        pos_to_rs = {}
        lead = None
        with gzip.open(path, "rt") as handle:
            for line in handle:
                raw += 1
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= C_GENE:
                    continue
                gid = parts[C_GENE]
                trait = parts[C_TRAIT]
                if not (gid.startswith(gene) or trait.startswith(gene)):
                    continue
                gene_n += 1
                p, beta, se, maf = fnum(parts[C_P]), fnum(parts[C_BETA]), fnum(parts[C_SE]), fnum(parts[C_MAF])
                if None in (p, beta, se, maf) or not (se > 0 and 0 < maf < 0.5):
                    continue
                filt_n += 1
                variants.add(parts[C_VARIANT])
                poskeys.add((parts[C_POS], parts[C_REF], parts[C_ALT]))
                rs = parts[C_RSID] if len(parts) > C_RSID else ""
                rsids.add(rs)
                pos_to_rs.setdefault(parts[C_POS], set()).add(rs)
                lead = p if lead is None else min(lead, p)
        multi_rs = sum(1 for rsset in pos_to_rs.values() if len(rsset) > 1)
        hist = stored.get(symbol, {})
        n_hist = hist.get("n_eqtl", "")
        thr_hist = hist.get("bonf_threshold", "")
        thr_now = 0.05 / filt_n if filt_n else ""
        passed_now = "" if lead is None or not filt_n else str(lead < thr_now)
        rows.append({
            **rec,
            "raw_rows": raw,
            "rows_after_gene_filter": gene_n,
            "rows_after_pvalue_beta_se_maf_filter": filt_n,
            "unique_variant_ids": len(variants),
            "unique_position_ref_alt": len(poskeys),
            "unique_rsids": len(rsids),
            "positions_with_multiple_rsids": multi_rs,
            "n_eqtl_meaning": "filtered cis-eQTL rows in the cached slice, not unique SNPs",
            "lead_p_recount": "" if lead is None else f"{lead:.8g}",
            "lead_p_historical": hist.get("lead_eqtl_p", ""),
            "n_eqtl_historical": n_hist,
            "n_eqtl_match": str(str(filt_n) == str(n_hist)),
            "bonferroni_threshold_0.05_over_filtered_rows": thr_now,
            "bonferroni_threshold_historical": thr_hist,
            "historical_threshold_equals_0.05_over_historical_n": (
                "" if not thr_hist or not n_hist else str(abs(float(thr_hist) - 0.05 / int(n_hist)) < 1e-18)
            ),
            "bonf_passed_from_recount": passed_now,
            "bonf_passed_historical": hist.get("bonf_passed", "not_in_source_table"),
            "is_eGene_atrial_historical": hist.get("is_eGene_atrial", "not_in_source_table"),
            "source_table": (
                "figure6_mr_anchors_v3_consolidated.tsv"
                if "bonf_threshold" in hist
                else "figure6_mr_positive_controls.tsv"
            ),
        })
    fields = list(rows[0].keys())
    # later rows may have extra keys
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    dest = OUT / "egene_cache_recount.tsv"
    with dest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print("wrote", dest)


def hca() -> None:
    path = ROOT / "global_raw.h5ad"
    checked = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    sha_line = (ROOT / "revision/acrf_provenance_check_20261015/logs/global_raw.h5ad.sha256").read_text().split()[0]
    with h5py.File(path, "r") as handle:
        n_obs = int(handle["X"]["indptr"].shape[0] - 1)
        n_var = int(handle["var"]["_index"].shape[0])
        obs_names = handle["obs"]["_index"][:]
        var_names = handle["var"]["_index"][:]

        def cats(name):
            raw = handle["obs"]["__categories"][name][:]
            return [x.decode() if isinstance(x, bytes) else str(x) for x in raw]

        def codes(name):
            return np.asarray(handle["obs"][name][:])

        donor, region, source, cell_source, used, cell_type = (
            cats(n) for n in ("donor", "region", "source", "cell_source", "Used", "cell_type")
        )
        keep_types = {
            "Atrial_Cardiomyocyte", "Fibroblast", "Endothelial", "Pericytes",
            "Smooth_muscle_cells", "Adipocytes", "Mesothelial", "Neuronal",
            "Myeloid", "Lymphoid",
        }
        src, reg, use, ct = codes("source"), codes("region"), codes("Used"), codes("cell_type")
        don = codes("donor")
        nuclei = src == source.index("Nuclei")
        la_ra = np.isin(reg, [region.index("LA"), region.index("RA")])
        used_yes = use == used.index("Yes")
        typed = np.isin(ct, [cell_type.index(x) for x in keep_types])
        mask = nuclei & la_ra & used_yes & typed
        donors_kept = sorted({donor[i] for i in np.unique(don[mask])})
        obs_fp = hashlib.sha256(b"\n".join(
            x.encode() if isinstance(x, str) else bytes(x) for x in obs_names
        )).hexdigest()
        var_fp = hashlib.sha256(b"\n".join(
            x.encode() if isinstance(x, str) else bytes(x) for x in var_names
        )).hexdigest()
        summary = {
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256_this_check": sha_line,
            "checked_at_utc": checked,
            "root_attrs_empty": len(handle.attrs) == 0,
            "n_obs": n_obs,
            "n_vars": n_var,
            "donors": donor,
            "regions": region,
            "source": source,
            "cell_source": cell_source,
            "cell_type": cell_type,
            "cells_after_nuclei_lara_used_mapped_type": int(mask.sum()),
            "donors_remaining_after_that_filter": donors_kept,
            "obs_name_sha256": obs_fp,
            "var_name_sha256": var_fp,
        }
    (OUT / "HCA_input_identity.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("hca cells", summary["cells_after_nuclei_lara_used_mapped_type"], "donors", donors_kept)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    recount()
    hca()
