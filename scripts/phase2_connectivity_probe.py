#!/usr/bin/env python3
"""Step 0 — probe network endpoints needed for Fig 6 Phase 2."""
from __future__ import annotations

import urllib.error
import urllib.request
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "phase2_connectivity_probe.tsv"

PROBES = [
    ("gwas_tabix_tbi", "GET",
     "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
     "GCST90624001-GCST90625000/GCST90624411/harmonised/GCST90624411.h.tsv.gz.tbi"),
    ("eqtl_atrial_tabix_tbi", "GET",
     "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/"
     "QTD000251/QTD000251.all.tsv.gz.tbi"),
    ("eqtl_blood_tabix_tbi", "GET",
     "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/"
     "QTD000356/QTD000356.all.tsv.gz.tbi"),
    ("gwas_catalog_api", "GET",
     "https://www.ebi.ac.uk/gwas/docs/"),
]


def probe(name: str, method: str, url: str) -> dict:
    req = urllib.request.Request(url, method=method, headers={"User-Agent": "ProjectP02/phase2"})
    try:
        with urllib.request.urlopen(req, timeout=20) as resp:
            return {
                "probe": name,
                "url": url,
                "ok": True,
                "status": resp.status,
                "note": resp.headers.get("Content-Type", ""),
            }
    except urllib.error.HTTPError as e:
        return {"probe": name, "url": url, "ok": False, "status": e.code, "note": str(e.reason)}
    except Exception as e:
        return {"probe": name, "url": url, "ok": False, "status": "", "note": str(e)}


def main() -> None:
    rows = [probe(n, m, u) for n, m, u in PROBES]
    df = pd.DataFrame(rows)
    df["probe_date"] = date.today().isoformat()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, sep="\t", index=False)
    n_ok = int(df["ok"].sum())
    print(f"wrote {OUT.relative_to(ROOT)} ({n_ok}/{len(df)} probes OK)")
    if n_ok < len(df):
        print("WARNING: some endpoints failed — Phase 2 will use curated gene-set fallback.")


if __name__ == "__main__":
    main()
