#!/usr/bin/env python3
"""
One-click generator for the innovation Figure 2 v3 pack.

Steps
-----
1) Build family map for v3 panel.
2) Render fig2b donor heatmap (v3).
3) Render fig2c chamber decomposition (v3).
4) Render fig2d family forest (v3).
5) Render fig2e replication lollipop (v3).

Outputs
-------
results/innovation_ncp_v3_axis_families.tsv
figures/fig2b_donor_acai_heatmap_v3.pdf
results/fig2c_chamber_decomposition_v3.tsv
figures/fig2c_chamber_decomposition_v3.pdf
results/fig2d_family_forest_v3.tsv
figures/fig2d_family_forest_v3.pdf
results/fig2e_replication_lollipop_v3.tsv
figures/fig2e_replication_lollipop_v3.pdf
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def run(cmd: list[str]) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=ROOT)


def main() -> None:
    py = sys.executable
    run([py, "scripts/v3/innovation_build_v3_axis_families.py"])
    run([
        py, "scripts/fig2b_donor_acai_heatmap.py",
        "--acai-tsv", "results/innovation_ncp_v3_acai_long.tsv",
        "--fam-tsv", "results/innovation_ncp_v3_axis_families.tsv",
        "--out-pdf", "figures/fig2b_donor_acai_heatmap_v3.pdf",
    ])
    run([
        py, "scripts/fig2c_chamber_decomposition.py",
        "--acai-tsv", "results/innovation_ncp_v3_acai_long.tsv",
        "--fam-tsv", "results/innovation_ncp_v3_axis_families.tsv",
        "--out-tsv", "results/fig2c_chamber_decomposition_v3.tsv",
        "--out-pdf", "figures/fig2c_chamber_decomposition_v3.pdf",
    ])
    run([
        py, "scripts/fig2d_family_forest.py",
        "--acai-tsv", "results/innovation_ncp_v3_acai_long.tsv",
        "--fam-tsv", "results/innovation_ncp_v3_axis_families.tsv",
        "--out-tsv", "results/fig2d_family_forest_v3.tsv",
        "--out-pdf", "figures/fig2d_family_forest_v3.pdf",
    ])
    run([py, "scripts/v3/innovation_fig2e_replication_lollipop_v3.py"])
    print("done: innovation fig2 v3 pack")


if __name__ == "__main__":
    main()
