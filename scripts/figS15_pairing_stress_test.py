#!/usr/bin/env python3
"""
FigS15: pairing stress test on the same healthy baseline cohort.

Question:
    Under pairing break / chamber imbalance, does naive chamber contrast drift
    more than paired aCAI?

Design:
    - Use the same GSE112339 healthy cohort and the same NCP_v2 panel as Fig1B.
    - Compute a per-program full-pair reference estimate across all 8 donors.
    - Simulate sampling scenarios:
        1) 8v8 matched (A=B)                  -> sanity check
        2) 6v6 broken pairing (A, B random)   -> moderate break
        3) 6v4 broken pairing                 -> imbalance + break
        4) 4v4 broken pairing                 -> severe break
    - Compare median absolute error vs reference for:
        naive  = mean(LA_A) - mean(RA_B)
        paired = mean(LA_d - RA_d) over overlap donors d in (A ∩ B)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from plot_style import COL_DOUBLE, FG_COLOR, LINE_COLOR, MUTED_COLOR, SUPPORT_COLOR, set_nature_style

set_nature_style()
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)

NCP_TSV = ROOT / "results" / "fig1c_ncp_panel_v2.tsv"
COUNTS = ROOT / "bulk" / "Project_P01" / "GSE112339" / "GSE112339_HTSeq_hg19_Ensembl.txt.gz"
GMT_HALL = ROOT / "bulk" / "Project_P01" / "gene_sets" / "h.all.v2026.1.Hs.symbols.gmt"
GMT_REAC = ROOT / "bulk" / "Project_P01" / "gene_sets" / "c2.cp.reactome.v2026.1.Hs.symbols.gmt"

OUT_LONG = ROOT / "results" / "figS15_pairing_stress_test.tsv"
OUT_SUM = ROOT / "results" / "figS15_pairing_stress_test_summary.tsv"
OUT_PDF = ROOT / "figures" / "figS15_pairing_stress_test.pdf"

N_RUNS = 2500
RNG_SEED = 20260422

SCENARIOS = [
    ("8v8 matched", 8, 8, True),
    ("6v6 broken", 6, 6, False),
    ("6v4 broken", 6, 4, False),
    ("4v4 broken", 4, 4, False),
]


def load_gmt(path: Path, prefix: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            out[f"{prefix}__{parts[0]}"] = [g for g in parts[2:] if g]
    return out


def load_log_cpm() -> pd.DataFrame:
    raw = pd.read_csv(COUNTS, sep="\t", index_col=0)
    raw = raw.groupby(raw.index).max()
    keep = [c for c in raw.columns if c.startswith(("LA_", "RA_"))]
    raw = raw[keep]
    cpm = raw.div(raw.sum(axis=0), axis=1) * 1e6
    return np.log2(cpm + 1.0)


def build_program_matrix(log_cpm: pd.DataFrame, pathways: dict[str, list[str]], ncp: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    la_cols = sorted([c for c in log_cpm.columns if c.startswith("LA_")])
    donors = [c[3:] for c in la_cols]

    rows = []
    for _, r in ncp.iterrows():
        gset = [g for g in pathways.get(str(r["pathway"]), []) if g in log_cpm.index]
        if len(gset) < 5:
            continue
        sub = log_cpm.loc[gset]
        row: dict[str, float | str] = {"pathway": str(r["pathway"])}
        for d in donors:
            row[f"LA_{d}"] = float(sub[f"LA_{d}"].mean())
            row[f"RA_{d}"] = float(sub[f"RA_{d}"].mean())
            row[f"D_{d}"] = float(row[f"LA_{d}"] - row[f"RA_{d}"])
        rows.append(row)
    return pd.DataFrame(rows), donors


def run_stress(df: pd.DataFrame, donors: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(RNG_SEED)
    ref = df[[f"D_{d}" for d in donors]].mean(axis=1).to_numpy(float)

    long_rows: list[dict[str, object]] = []
    for scenario, n_la, n_ra, matched in SCENARIOS:
        for run_idx in range(1, N_RUNS + 1):
            if matched:
                a = rng.choice(donors, size=n_la, replace=False)
                b = a.copy()
            else:
                a = rng.choice(donors, size=n_la, replace=False)
                b = rng.choice(donors, size=n_ra, replace=False)

            a = [str(x) for x in a]
            b = [str(x) for x in b]
            overlap = sorted(set(a) & set(b))
            if not overlap:
                continue

            naive = (
                df[[f"LA_{d}" for d in a]].mean(axis=1).to_numpy(float)
                - df[[f"RA_{d}" for d in b]].mean(axis=1).to_numpy(float)
            )
            paired = df[[f"D_{d}" for d in overlap]].mean(axis=1).to_numpy(float)

            err_naive = float(np.median(np.abs(naive - ref)))
            err_paired = float(np.median(np.abs(paired - ref)))

            long_rows.append(
                {
                    "scenario": scenario,
                    "run_id": run_idx,
                    "method": "naive",
                    "median_abs_err": err_naive,
                    "overlap_n": len(overlap),
                    "n_la": n_la,
                    "n_ra": n_ra,
                }
            )
            long_rows.append(
                {
                    "scenario": scenario,
                    "run_id": run_idx,
                    "method": "paired",
                    "median_abs_err": err_paired,
                    "overlap_n": len(overlap),
                    "n_la": n_la,
                    "n_ra": n_ra,
                }
            )

    long_df = pd.DataFrame(long_rows)
    summary = (
        long_df.groupby(["scenario", "method"], as_index=False)
        .agg(
            n_runs=("median_abs_err", "size"),
            median_abs_err=("median_abs_err", "median"),
            p90_abs_err=("median_abs_err", lambda x: float(np.quantile(x, 0.9))),
            mean_overlap=("overlap_n", "mean"),
        )
        .sort_values(["scenario", "method"])
        .reset_index(drop=True)
    )
    return long_df, summary


def plot(long_df: pd.DataFrame, summary: pd.DataFrame) -> None:
    order = [s[0] for s in SCENARIOS]
    fig, (ax1, ax2) = plt.subplots(
        1,
        2,
        figsize=(COL_DOUBLE * 1.06, COL_DOUBLE * 0.48),
        gridspec_kw={"width_ratios": [1.35, 1.0], "wspace": 0.22},
    )

    # Panel a: scenario-wise error distribution
    for i, scen in enumerate(order):
        dn = long_df[(long_df["scenario"] == scen) & (long_df["method"] == "naive")]["median_abs_err"].to_numpy(float)
        dp = long_df[(long_df["scenario"] == scen) & (long_df["method"] == "paired")]["median_abs_err"].to_numpy(float)
        bp = ax1.boxplot(
            [dn, dp],
            positions=[i - 0.18, i + 0.18],
            widths=0.28,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": FG_COLOR, "linewidth": 1.1},
            boxprops={"linewidth": 0.8},
            whiskerprops={"linewidth": 0.8},
            capprops={"linewidth": 0.8},
        )
        bp["boxes"][0].set_facecolor("#D8D8D8")
        bp["boxes"][1].set_facecolor(SUPPORT_COLOR)

    ax1.set_xticks(np.arange(len(order)))
    ax1.set_xticklabels(order, rotation=20, ha="right", fontsize=6.8)
    ax1.set_ylabel("median absolute error vs full-pair reference", fontsize=7.2)
    ax1.set_title("a  Pairing-break stress test", loc="left", fontsize=9.4, fontweight="bold")
    ax1.grid(axis="y", alpha=0.2)
    for s in ("top", "right"):
        ax1.spines[s].set_visible(False)

    # Annotate medians
    for i, scen in enumerate(order):
        s_na = summary[(summary["scenario"] == scen) & (summary["method"] == "naive")].iloc[0]
        s_pa = summary[(summary["scenario"] == scen) & (summary["method"] == "paired")].iloc[0]
        ax1.text(i - 0.18, float(s_na["median_abs_err"]) + 0.002, f"{float(s_na['median_abs_err']):.3f}", fontsize=6.3, ha="center", va="bottom", color=MUTED_COLOR)
        ax1.text(i + 0.18, float(s_pa["median_abs_err"]) + 0.002, f"{float(s_pa['median_abs_err']):.3f}", fontsize=6.3, ha="center", va="bottom", color=SUPPORT_COLOR)

    # Panel b: error vs overlap
    sub = long_df[long_df["scenario"] != "8v8 matched"].copy()
    g = (
        sub.groupby(["method", "overlap_n"], as_index=False)["median_abs_err"]
        .median()
        .sort_values(["method", "overlap_n"])
    )
    for method, color, marker in [("naive", "#666666", "o"), ("paired", SUPPORT_COLOR, "s")]:
        gg = g[g["method"] == method]
        ax2.plot(gg["overlap_n"], gg["median_abs_err"], marker=marker, color=color, lw=1.4, ms=4.5, label=method)

    ax2.set_xlabel("overlap donors between LA/RA sets", fontsize=7.2)
    ax2.set_ylabel("median absolute error", fontsize=7.2)
    ax2.set_title("b  Error drops with overlap; paired remains lower", loc="left", fontsize=9.2, fontweight="bold")
    ax2.grid(alpha=0.22)
    ax2.legend(frameon=False, fontsize=6.8, loc="upper right")
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)

    fig.text(
        0.012,
        0.01,
        "Naive = mean(LA_A)-mean(RA_B); Paired = mean(LA_d-RA_d) over overlap donors d in A∩B.",
        fontsize=6.4,
        color=MUTED_COLOR,
        ha="left",
        va="bottom",
    )
    fig.subplots_adjust(left=0.07, right=0.99, top=0.92, bottom=0.16, wspace=0.22)
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, bbox_inches=None, pad_inches=0.04)
    plt.close(fig)


def main() -> None:
    for p in [NCP_TSV, COUNTS, GMT_HALL, GMT_REAC]:
        if not p.exists():
            raise FileNotFoundError(f"Missing input: {p.relative_to(ROOT)}")

    ncp = pd.read_csv(NCP_TSV, sep="\t")
    paths = {}
    paths.update(load_gmt(GMT_HALL, "HALLMARK"))
    paths.update(load_gmt(GMT_REAC, "REACTOME"))
    log_cpm = load_log_cpm()
    prog_df, donors = build_program_matrix(log_cpm, paths, ncp)

    long_df, summary = run_stress(prog_df, donors)
    long_df.to_csv(OUT_LONG, sep="\t", index=False, float_format="%.6f")
    summary.to_csv(OUT_SUM, sep="\t", index=False, float_format="%.6f")
    plot(long_df, summary)

    print(f"wrote {OUT_LONG.relative_to(ROOT)}")
    print(f"wrote {OUT_SUM.relative_to(ROOT)}")
    print(f"wrote {OUT_PDF.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
