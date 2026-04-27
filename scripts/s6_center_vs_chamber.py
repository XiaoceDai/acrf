"""
S6 — center effect vs chamber effect, descriptive only (no DA test).

Per top-4 cell type: each donor contributes a paired (LA, RA) dot, connected
by a line. Donors coloured by center: D-prefix = Sanger, H-prefix = Harvard.

Output:
    figures/fig_center_vs_chamber.pdf
    figures/center_vs_chamber_stats.tsv
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from plot_style import COL_DOUBLE, MM, panel_label, set_nature_style, strip_spines

H5AD = Path("/home/bio/Project_P02/global_raw.h5ad")
FIG_DIR = Path("/home/bio/Project_P02/figures")
FIG_DIR.mkdir(parents=True, exist_ok=True)

TOP4 = {
    "Atrial_Cardiomyocyte":    "atrial_cardiomyocyte",
    "Fibroblast":              "fibroblast",
    "Endothelial":             "endothelial_or_vascular",
    "Myeloid":                 "immune",
    "Lymphoid":                "immune",
}
TOP4_ORDER = ["atrial_cardiomyocyte", "fibroblast", "endothelial_or_vascular", "immune"]
TOP4_PRETTY = {
    "atrial_cardiomyocyte":    "atrial CM",
    "fibroblast":              "fibroblast",
    "endothelial_or_vascular": "endo / vascular",
    "immune":                  "immune",
}

CENTER_COLORS = {
    "Sanger":  "#1f6fb4",
    "Harvard": "#e07b2a",
}


def center_of(d: str) -> str:
    return "Sanger" if str(d).startswith("D") else "Harvard"


def main() -> None:
    set_nature_style()

    a = ad.read_h5ad(H5AD, backed="r")
    obs = a.obs[["donor", "region", "cell_type", "source", "Used"]].copy()
    obs = obs[
        (obs["source"].astype(str) == "Nuclei")
        & obs["region"].isin(["LA", "RA"])
        & (obs["Used"].astype(str) == "Yes")
    ]
    obs["cell_type"] = obs["cell_type"].astype(str)
    obs["top_ct"] = obs["cell_type"].map(TOP4)
    obs = obs[obs["top_ct"].notna()]
    obs["region"] = obs["region"].astype(str)
    obs["center"] = obs["donor"].astype(str).map(center_of)

    g = (obs.groupby(["donor", "region", "top_ct"], observed=True)
            .size().rename("n").reset_index())
    g["pct"] = g["n"] / g.groupby(["donor", "region"], observed=True)["n"].transform("sum") * 100
    g["center"] = g["donor"].map(center_of)

    have_both = (g.groupby("donor", observed=True)["region"]
                    .apply(lambda r: {"LA", "RA"}.issubset(set(r))))
    paired_donors = sorted(have_both[have_both].index.tolist(),
                           key=lambda d: (center_of(d), d))
    g = g[g["donor"].isin(paired_donors)].copy()

    fig = plt.figure(figsize=(COL_DOUBLE, 62 * MM))
    gs = fig.add_gridspec(1, 4, wspace=0.35,
                          left=0.06, right=0.98, top=0.78, bottom=0.28)

    stats_rows = []

    for k, ct in enumerate(TOP4_ORDER):
        ax = fig.add_subplot(gs[0, k])
        sub = g[g["top_ct"] == ct]
        w = sub.pivot(index="donor", columns="region", values="pct").reindex(index=paired_donors)
        w_center = pd.Series({d: center_of(d) for d in w.index})

        jitter = {"Sanger": -0.06, "Harvard": +0.06}
        x_pos = {"LA": 0.0, "RA": 1.0}

        for donor in w.index:
            c = w_center[donor]
            ax.plot(
                [x_pos["LA"] + jitter[c], x_pos["RA"] + jitter[c]],
                [w.loc[donor, "LA"],      w.loc[donor, "RA"]],
                color=CENTER_COLORS[c], alpha=0.55, lw=0.7, zorder=1,
            )

        for c, colr in CENTER_COLORS.items():
            mask = w_center == c
            if mask.sum() == 0:
                continue
            for region in ("LA", "RA"):
                ax.scatter(
                    np.repeat(x_pos[region] + jitter[c], mask.sum()),
                    w.loc[mask, region].values,
                    s=14, color=colr, edgecolor="white", linewidth=0.4,
                    zorder=3, label=c if (region == "LA" and k == 0) else None,
                )

        for c, colr in CENTER_COLORS.items():
            mask = w_center == c
            if mask.sum() == 0:
                continue
            for region in ("LA", "RA"):
                m = w.loc[mask, region].mean()
                x0 = x_pos[region] + jitter[c]
                ax.plot([x0 - 0.07, x0 + 0.07], [m, m],
                        color=colr, lw=1.4, solid_capstyle="round", zorder=4)

        ax.set_xticks([0, 1]); ax.set_xticklabels(["LA", "RA"])
        ax.set_xlim(-0.35, 1.35)
        ax.set_ylim(bottom=0)
        ax.set_title(TOP4_PRETTY[ct])
        if k == 0:
            ax.set_ylabel("% of top-4 cells\n(within donor × chamber)")
        strip_spines(ax)
        panel_label(ax, "abcd"[k], x=-0.12)

        chamber_diff = (w["RA"] - w["LA"]).abs().median()
        center_diff = abs(w.loc[w_center == "Sanger"].mean().mean()
                          - w.loc[w_center == "Harvard"].mean().mean())
        stats_rows.append({
            "top_ct": ct,
            "median_|RA-LA|_within_donor_pct": round(float(chamber_diff), 2),
            "|Sanger_mean - Harvard_mean|_pct": round(float(center_diff), 2),
        })

        ratio = (stats_rows[-1]["|Sanger_mean - Harvard_mean|_pct"]
                 / max(stats_rows[-1]["median_|RA-LA|_within_donor_pct"], 1e-9))
        ax.text(0.5, -0.22,
                f"center Δ / chamber Δ ≈ {ratio:.1f}×",
                transform=ax.transAxes, ha="center", va="top",
                fontsize=5.5, color="#444")

    fig.text(0.06, 0.935, "Center effect vs chamber effect (top-4 classes)",
             fontsize=8, fontweight="bold")
    fig.legend(
        handles=[plt.Line2D([0], [0], marker="o", color="w",
                            markerfacecolor=CENTER_COLORS[c],
                            markeredgecolor="white", markersize=5, label=c)
                 for c in ("Sanger", "Harvard")],
        loc="upper right", bbox_to_anchor=(0.98, 0.96),
        frameon=False, ncol=2, fontsize=6, handlelength=0.8, columnspacing=0.8,
    )
    fig.text(0.06, 0.035,
             "Each dot = one donor × chamber. Lines connect paired LA/RA within a "
             "donor. Short dashes = center means. Slope = chamber effect; vertical "
             "offset between Sanger and Harvard = center effect. Descriptive only; "
             "no statistical testing (descriptive-only panel).",
             fontsize=5.5, color="#444")

    out = FIG_DIR / "fig_center_vs_chamber.pdf"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")

    pd.DataFrame(stats_rows).to_csv(
        FIG_DIR / "center_vs_chamber_stats.tsv", sep="\t", index=False
    )
    print(f"wrote {FIG_DIR / 'center_vs_chamber_stats.tsv'}")


if __name__ == "__main__":
    main()
