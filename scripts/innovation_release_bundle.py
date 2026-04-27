#!/usr/bin/env python3
"""
Build a reproducible innovation release bundle.

This script runs the full innovation pipeline in a fixed order and then
writes a manifest (with SHA256 checksums) for core outputs.

Outputs
-------
results/innovation_release_manifest.tsv
results/innovation_release_summary.tsv
"""
from __future__ import annotations

import hashlib
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

OUT_MANIFEST = ROOT / "results" / "innovation_release_manifest.tsv"
OUT_SUMMARY = ROOT / "results" / "innovation_release_summary.tsv"


def run(cmd: list[str]) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=ROOT)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    py = sys.executable

    # 1) Core innovation panel + robustness
    run([py, "scripts/innovation_ncp_consensus.py"])
    run([py, "scripts/innovation_ncp_robustness.py"])

    # 2) v2/v3 transferability benchmark
    run([
        py, "scripts/innovation_recompute_acai_panel.py",
        "--panel", "results/fig1c_ncp_panel_v2.tsv",
        "--tag", "innovation_ncp_v2",
    ])
    run([
        py, "scripts/innovation_recompute_acai_panel.py",
        "--panel", "results/innovation_ncp_v3_panel.tsv",
        "--tag", "innovation_ncp_v3",
    ])
    run([py, "scripts/v3/innovation_v2_vs_v3_benchmark.py"])

    # 3) Fig2 v3 pack
    run([py, "scripts/v3/innovation_generate_fig2_v3_pack.py"])

    output_files = [
        ("result", "results/innovation_ncp_consensus_long.tsv"),
        ("result", "results/innovation_ncp_v3_panel.tsv"),
        ("result", "results/innovation_ncp_v3_summary.tsv"),
        ("result", "results/innovation_ncp_bootstrap_frequency.tsv"),
        ("result", "results/innovation_ncp_bootstrap_summary.tsv"),
        ("result", "results/innovation_ncp_v2_acai_long.tsv"),
        ("result", "results/innovation_ncp_v2_acai_summary.tsv"),
        ("result", "results/innovation_ncp_v3_acai_long.tsv"),
        ("result", "results/innovation_ncp_v3_acai_summary.tsv"),
        ("result", "results/innovation_v2_vs_v3_metrics.tsv"),
        ("result", "results/innovation_v2_vs_v3_axis_effects.tsv"),
        ("result", "results/innovation_v2_vs_v3_decision_table.tsv"),
        ("result", "results/innovation_ncp_v3_axis_families.tsv"),
        ("result", "results/fig2c_chamber_decomposition_v3.tsv"),
        ("result", "results/fig2d_family_forest_v3.tsv"),
        ("result", "results/fig2e_replication_lollipop_v3.tsv"),
        ("figure", "figures/innovation_ncp_bootstrap_frequency.pdf"),
        ("figure", "figures/innovation_v2_vs_v3_benchmark.pdf"),
        ("figure", "figures/fig2b_donor_acai_heatmap_v3.pdf"),
        ("figure", "figures/fig2c_chamber_decomposition_v3.pdf"),
        ("figure", "figures/fig2d_family_forest_v3.pdf"),
        ("figure", "figures/fig2e_replication_lollipop_v3.pdf"),
    ]

    rows: list[dict[str, object]] = []
    for kind, rel in output_files:
        p = ROOT / rel
        if not p.exists():
            raise FileNotFoundError(f"Expected output missing: {rel}")
        st = p.stat()
        rows.append(
            {
                "kind": kind,
                "path": rel,
                "bytes": int(st.st_size),
                "sha256": sha256_file(p),
                "mtime_utc": datetime.fromtimestamp(st.st_mtime, tz=timezone.utc).isoformat(),
            }
        )

    manifest = pd.DataFrame(rows)
    manifest.to_csv(OUT_MANIFEST, sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {
                "released_at_utc": datetime.now(tz=timezone.utc).isoformat(),
                "n_files": int(len(manifest)),
                "n_result_files": int((manifest["kind"] == "result").sum()),
                "n_figure_files": int((manifest["kind"] == "figure").sum()),
            }
        ]
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    print(f"wrote {OUT_MANIFEST.relative_to(ROOT)} ({len(manifest)} files)")
    print(f"wrote {OUT_SUMMARY.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
