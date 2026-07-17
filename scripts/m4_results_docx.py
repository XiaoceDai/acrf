#!/usr/bin/env python3
"""
Build the Nature Communications-style Results Word document.

Output: manuscript_drafts/Results.docx

Style is shared with Methods.docx via m0_doc_style.

Phase 4A/4A-P: Results preamble + Figure 1–2 are built at runtime and locked.
Phase 4B/4B-P: Figure 3–4 narratives are built at runtime and locked.
Phase 4C: Figure 5–6 narratives are built at runtime from final source tables.
Post-Figure-6 Results text (Specificity) remains frozen pending later migration.
"""
from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Pt, RGBColor

from m0_doc_style import apply_base_style, add_body

OUT = ROOT / "manuscript_drafts" / "Results.docx"
OUT.parent.mkdir(parents=True, exist_ok=True)


def _fid(n: int) -> str:
    return f"F{int(n)}"


def _uminus(x: float, nd: int = 3) -> str:
    return f"{x:.{nd}f}".replace("-", "−")


def build_results_preamble_fig1_fig2() -> list[tuple[str, list[str]]]:
    """Load final Fig1/Fig2 sources, assert locks, return Overview+Fig1+Fig2 sections."""
    import pandas as pd

    mag_path = ROOT / "results/ncp_v3_singscore_reanalysis/fig1b_healthy_reference_direction_magnitude_summary.tsv"
    conc_path = ROOT / "results/ncp_v3_singscore_reanalysis/fig1b_healthy_reference_concordance_summary.tsv"
    qc_path = ROOT / "results/fig1c_selection_provenance_qc.md"
    dep_path = ROOT / "results/ncp_v3_singscore_reanalysis/fig2a_bulk_deployment_summary.tsv"
    st_sum_path = ROOT / "results/ncp_v3_singscore_reanalysis/fig2c_family_cohort_states_summary.tsv"
    st_src_path = ROOT / "results/ncp_v3_singscore_reanalysis/fig2c_family_cohort_states_source.tsv"
    cons_sum_path = ROOT / "results/ncp_v3_singscore_reanalysis/fig2d_family_effects_consistency_summary.tsv"
    cons_src_path = ROOT / "results/ncp_v3_singscore_reanalysis/fig2d_family_effects_consistency_source.tsv"
    fam_path = ROOT / "results/ncp_v3_family_labels.tsv"

    mag_lines = mag_path.read_text(encoding="utf-8").splitlines()
    mag_row = dict(zip(mag_lines[0].split("\t"), mag_lines[1].split("\t")))
    n_programs = int(mag_row["n_programs"])
    n_healthy = int(mag_row["n_healthy_donors"])
    n_la = int(mag_row["n_la_biased"])
    n_ra = int(mag_row["n_ra_biased"])
    n_conc = int(mag_row["n_programs_8_of_8_concordant"])
    median_m_full = mag_row["median_M"]
    if not median_m_full.startswith("0.028354720506272"):
        raise AssertionError(f"unexpected median_M: {median_m_full}")
    median_m = float(median_m_full)
    min_m = float(mag_row["minimum_M"])
    max_m = float(mag_row["maximum_M"])
    if (n_programs, n_la, n_ra, n_conc, n_healthy) != (30, 22, 8, 30, 8):
        raise AssertionError("Figure 1B composition lock failed")

    conc = pd.read_csv(conc_path, sep="\t").iloc[0]
    if int(conc["n_programs"]) != 30 or int(conc["n_families"]) != 16:
        raise AssertionError("Figure 1 concordance/topology lock failed")
    if int(conc["n_full_concordance"]) != 30:
        raise AssertionError("Figure 1 concordance count lock failed")

    qc = qc_path.read_text(encoding="utf-8")
    n_strict = int(re.search(r"strict-selected programs:\s*(\d+)", qc).group(1))
    n_fallback = int(re.search(r"fallback-selected programs:\s*(\d+)", qc).group(1))
    n_cand = int(re.search(r"unique healthy candidates:\s*(\d+)", qc).group(1))
    n_elig = int(re.search(r"strict-eligible candidates:\s*(\d+)", qc).group(1))
    if (n_strict, n_fallback) != (25, 5):
        raise AssertionError(f"strict/fallback lock failed: {n_strict}/{n_fallback}")
    if (n_cand, n_elig) != (80, 38):
        raise AssertionError("candidate/strict-pool lock failed")
    if n_strict + n_fallback != 30:
        raise AssertionError("25+5 != 30")

    dep = pd.read_csv(dep_path, sep="\t")
    tot = dep.loc[dep["scope"] == "total"].iloc[0]
    n_paired = int(tot["paired_donors"])
    n_af = int(tot["AF_donors"])
    n_sr = int(tot["SR_donors"])
    if (n_paired, n_af, n_sr) != (51, 25, 26):
        raise AssertionError("Figure 2A donor lock failed")
    if bool(tot["cohort_specific_reselection"]):
        raise AssertionError("unexpected pathway reselection flag")

    st_src = pd.read_csv(st_src_path, sep="\t")
    if st_src.duplicated(["family_id", "cohort"]).sum() != 0:
        raise AssertionError("duplicate family×cohort state rows")
    sc = Counter(st_src["state"].tolist())
    expected_states = {"Amplified": 5, "Preserved": 16, "Attenuated": 15, "Reversed": 12}
    if dict(sc) != expected_states or sum(sc.values()) != 48:
        raise AssertionError(f"state count lock failed: {dict(sc)}")
    st_sum = pd.read_csv(st_sum_path, sep="\t")
    for _, r in st_sum.iterrows():
        if int(r["n_calls"]) != sc[r["state"]]:
            raise AssertionError("state summary/source mismatch")

    cons_src = pd.read_csv(cons_src_path, sep="\t")
    pooled_rows = cons_src[cons_src["estimate_type"] == "pooled"].copy()
    if pooled_rows.duplicated("family_id").sum() != 0:
        raise AssertionError("duplicate pooled family rows")
    if set(pooled_rows["pooling_model"].dropna().unique()) != {"fixed-effect_inverse_variance"}:
        raise AssertionError("primary pooling model is not fixed-effect inverse-variance")
    by = {
        k: sorted(_fid(x) for x in g["family_id"].tolist())
        for k, g in pooled_rows.groupby("consistency_class")
    }
    high_order = ["F1", "F9", "F10", "F12", "F14", "F15"]
    mod_order = ["F2", "F3", "F5", "F6", "F8", "F11", "F13", "F16"]
    low_order = ["F4", "F7"]
    if by.get("High") != sorted(high_order):
        raise AssertionError(f"High membership mismatch: {by.get('High')}")
    if by.get("Moderate") != sorted(mod_order):
        raise AssertionError(f"Moderate membership mismatch: {by.get('Moderate')}")
    if by.get("Low") != sorted(low_order):
        raise AssertionError(f"Low membership mismatch: {by.get('Low')}")

    cons_sum = pd.read_csv(cons_sum_path, sep="\t").copy()
    cons_sum["fid"] = cons_sum["family_id"].map(_fid)
    high_g = cons_sum.set_index("fid")
    expected_g = {
        "F1": -0.608779,
        "F9": -0.630704,
        "F10": -0.799070,
        "F12": -0.702093,
        "F14": -0.925194,
        "F15": -1.116983,
    }
    for fid, eg in expected_g.items():
        if abs(float(high_g.loc[fid, "pooled_g"]) - eg) >= 1e-5:
            raise AssertionError(f"pooled g mismatch for {fid}")
    if not (high_g.loc[high_order, "pooled_g"] < 0).all():
        raise AssertionError("not all High pooled effects are negative")

    f15 = high_g.loc["F15"]
    g15 = float(f15["pooled_g"])
    lo15 = float(f15["pooled_ci_lower"])
    hi15 = float(f15["pooled_ci_upper"])
    if (
        abs(g15 - (-1.116983)) >= 1e-5
        or abs(lo15 - (-1.708325)) >= 1e-5
        or abs(hi15 - (-0.525642)) >= 1e-5
    ):
        raise AssertionError("F15 pooled estimate/CI mismatch")

    fam = pd.read_csv(fam_path, sep="\t")
    lab = {
        _fid(r.family_id): str(r.family_label).replace("signalling", "signaling")
        for _, r in fam.iterrows()
    }

    median_disp = f"{median_m:.4f}"
    if median_disp != "0.0284":
        raise AssertionError(f"median display lock failed: {median_disp}")
    range_text = f"approximately {min_m:.5f} to {max_m:.3f}"

    # Donor-unit lock: 30/30 programs concordant in all 8 healthy donors (not "donor pairs").
    if not (n_programs == 30 and n_conc == 30 and n_healthy == 8 and n_la == 22 and n_ra == 8):
        raise AssertionError("Figure 1B donor-unit / concordance lock failed")
    if median_disp != "0.0284":
        raise AssertionError("median M display lock failed during prose emit")

    overview = (
        "Overview",
        [
            (
                "We organized the analyses around a single frozen reference, the Normal Chamber "
                "Prior (NCP): 30 programs organized into 16 families, calibrated in healthy donors "
                "with paired atrial samples and carried forward with fixed definitions across "
                "disease datasets. Subsequent sections examine the healthy reference, paired bulk "
                "deployment across AF and sinus-rhythm cohorts, LA single-nucleus localization, "
                "cross-sectional duration context, cross-layer triangulation, and genetic "
                "calibration as distinct evidence layers with layer-specific interpretation "
                "boundaries."
            )
        ],
    )

    fig1 = (
        "A frozen healthy chamber reference defines the NCP",
        [
            (
                "Figure 1A defines the donor-paired estimand used throughout. Rhythm is "
                "assigned at the donor level, whereas chamber is measured within donor. "
                "The primary paired-bulk analysis unit is therefore a donor-level LA–RA program "
                "contrast. Each contrast is aligned to the frozen healthy chamber direction "
                "and scaled by the corresponding healthy reference magnitude, yielding a "
                "relative aligned chamber asymmetry that permits comparison across programs. "
                "Figure 1A is a definition of the measurement framework, not an AF association result."
            ),
            (
                "In the healthy reference (Figure 1B), the 30 frozen programs were scored in eight "
                "healthy donors with paired LA and RA samples. Twenty-two programs were LA-biased "
                "and eight were RA-biased. All 30 programs matched the canonical chamber direction "
                "in every donor. Healthy reference magnitudes varied substantially across programs "
                f"(median M = {median_disp}; range {range_text}), motivating "
                "program-specific scaling. These results calibrate the healthy chamber reference; "
                "they are not AF discovery findings."
            ),
            (
                f"Figure 1C summarizes selection provenance. {n_cand} candidate pathways were "
                f"evaluated, and {n_elig} met the strict gate and entered the strict pool. "
                f"After diversity-aware redundancy control, {n_strict} programs were retained "
                f"from the strict pool, and {n_fallback} additional programs were admitted "
                f"under the prespecified fallback rule to complete the frozen {n_programs}-program "
                f"panel. Fallback admission was prespecified. Disease cohorts did not contribute "
                f"to candidate scoring or reselection, and AF cohorts did not alter program "
                f"membership or chamber direction."
            ),
            (
                f"Figure 1D maps the frozen topology used downstream: the {n_programs} programs "
                f"are organized into {int(conc['n_families'])} families in canonical display order, "
                f"retaining the {n_la} LA-biased and {n_ra} RA-biased programs defined above. "
                "These families span translation and ribosome quality control, "
                "extracellular-matrix and glycosaminoglycan biology, growth-factor and receptor "
                "signaling, metabolic programs, and related atrial-remodeling themes. Figure 1 "
                "therefore fixes the estimand, healthy calibration, selection provenance, and "
                "program-to-family topology used in all subsequent analyses."
            ),
        ],
    )

    # Assert High membership and F15 estimate before emitting polished Fig2D prose.
    if high_order != ["F1", "F9", "F10", "F12", "F14", "F15"]:
        raise AssertionError("High membership changed unexpectedly")
    if mod_order != ["F2", "F3", "F5", "F6", "F8", "F11", "F13", "F16"]:
        raise AssertionError("Moderate membership changed unexpectedly")
    if low_order != ["F4", "F7"]:
        raise AssertionError("Low membership changed unexpectedly")
    if _uminus(g15) != "−1.117" or _uminus(lo15) != "−1.708" or _uminus(hi15) != "−0.526":
        raise AssertionError("F15 display estimate/CI changed unexpectedly")

    fig2 = (
        "Paired-chamber bulk cohorts reveal selective cross-cohort remodeling",
        [
            (
                f"We next deployed the frozen NCP in three paired-chamber bulk cohorts "
                f"totaling {n_paired} paired donors ({n_af} AF; {n_sr} sinus rhythm). "
                f"Program definitions and leading-edge gene sets remained fixed. For each "
                f"donor, paired chamber contrasts were direction-aligned and scaled by "
                f"healthy reference magnitudes, then summarized as family-level relative "
                f"aligned chamber asymmetry index (aCAIrel) without pathway reselection in disease "
                f"cohorts (Figure 2A). Figure 2B shows donor-level heterogeneity together with "
                f"cohort and family structure, with AF-associated shifts superimposed on these "
                f"sources of variation. The three cohorts were analyzed as parallel "
                f"paired-chamber deployments, without designating any cohort as a "
                f"discovery or replication arm."
            ),
            (
                f"At the family-by-cohort level (Figure 2C), descriptive state assignments "
                f"across {sum(sc.values())} family×cohort calls comprised "
                f"{sc['Amplified']} Amplified, {sc['Preserved']} Preserved, "
                f"{sc['Attenuated']} Attenuated, and {sc['Reversed']} Reversed "
                f"({sc['Amplified']} + {sc['Preserved']} + {sc['Attenuated']} + "
                f"{sc['Reversed']} = {sum(sc.values())}). These labels summarize observed "
                f"remodeling patterns and are not independent significance tests; a shared "
                f"state label does not imply a unique causal trajectory."
            ),
            (
                "Pooled AF-minus-SR family effects were estimated by fixed-effect "
                "inverse-variance weighting and graded for cross-cohort consistency "
                "(Figure 2D). Consistency classes comprised six High, eight Moderate, and two "
                "Low families. The High group included three LA-biased translation/RQC "
                "families—F1 Translation / RQC, F12 rRNA processing, and F14 Translation "
                "initiation—and three RA-biased matrix/developmental families—F9 HS-GAG / bone "
                "/ GPCR, F10 GAG disease / TGF-β, and F15 NTRK2 / TrkB. Moderate-consistency "
                "families were F2, F3, F5, F6, F8, F11, F13, and F16; Low-consistency families "
                "were F4 and F7. All six High-consistency pooled AF-minus-SR effects were "
                "negative on the aligned family scale. The largest-magnitude High estimate was "
                f"observed for F15 NTRK2 / TrkB (Hedges’ g = {_uminus(g15)}; 95% CI "
                f"{_uminus(lo15)} to {_uminus(hi15)}). Negative pooled shifts indicate "
                "reproducible AF-minus-SR displacement on the aligned scale and should not "
                "be translated automatically into a single mechanism such as reversal."
            ),
            (
                "Together, these results show that cross-cohort reproducibility was selective "
                "rather than uniform. High-consistency remodeling was not confined to "
                "RA-biased families; it included both LA-biased translation/RQC families and "
                "RA-biased matrix/developmental families. Metabolic and mitochondrial families "
                "were more often Moderate or Low, but broad cohort-dependence claims for "
                "translational families as a class are not supported. Consistency classes "
                "are evidence summaries, not causal claims."
            ),
        ],
    )

    return [overview, fig1, fig2]



def build_results_fig3_fig4() -> list[tuple[str, list[str]]]:
    """Load final Fig3/Fig4 sources, assert locks, return Figure 3–4 Results sections."""
    import pandas as pd

    fid_classes = ROOT / "results/fig3c_family_fidelity_classes.tsv"
    carrier = ROOT / "results/fig3_v3_carrier_fidelity_class.tsv"
    bridge = ROOT / "results/fig3a_localization_bridge_source.tsv"
    fam_lab = ROOT / "results/ncp_v3_family_labels.tsv"
    fig3d_sum = ROOT / "results/ncp_v3_singscore_reanalysis/fig3d_relay_extension_summary.tsv"
    fig3c_sum = ROOT / "results/ncp_v3_singscore_reanalysis/fig3c_fidelity_relay_summary.tsv"
    design4 = ROOT / "results/ncp_v3_singscore_reanalysis/fig4a_duration_design_summary.tsv"
    chron = ROOT / "results/ncp_v3_singscore_reanalysis/fig4_v3_chronicity_class.tsv"
    laa = ROOT / "results/ncp_v3_singscore_reanalysis/gse41177_family_duration_laa_primary.tsv"
    lapv = ROOT / "results/ncp_v3_singscore_reanalysis/gse41177_family_duration_lapv_sensitivity.tsv"
    fig4d = ROOT / "results/ncp_v3_singscore_reanalysis/fig4d_representative_duration_examples_summary.tsv"

    bridge_df = pd.read_csv(bridge, sep="\t")
    bmap = {
        str(r["item"]): str(r["value"])
        for _, r in bridge_df.iterrows()
    }
    if bmap.get("dataset") != "GSE255612":
        raise AssertionError("Figure 3 dataset lock failed")
    if int(bmap["n_AF"]) != 18 or int(bmap["n_controls"]) != 16:
        raise AssertionError("Figure 3 donor locks failed")
    if bmap.get("chamber") != "LA-only snRNA-seq":
        raise AssertionError("Figure 3 LA-only lock failed")
    if int(bmap["n_compartments"]) != 9:
        raise AssertionError("Figure 3 compartment-count lock failed")

    classes = pd.read_csv(fid_classes, sep="\t")
    carrier_df = pd.read_csv(carrier, sep="\t")
    by_class = {
        k: sorted(int(x) for x in g["family_id"].tolist())
        for k, g in classes.groupby("fidelity_class")
    }
    if by_class.get("Fidelity") != [10]:
        raise AssertionError(f"Fidelity lock failed: {by_class.get('Fidelity')}")
    if by_class.get("Relay") != [11, 15]:
        raise AssertionError(f"Relay lock failed: {by_class.get('Relay')}")
    if by_class.get("Extension") != [16]:
        raise AssertionError(f"Extension lock failed: {by_class.get('Extension')}")
    if len(by_class.get("Below_detection", [])) != 12:
        raise AssertionError("Below-detection count lock failed")
    non_below = sorted(
        by_class["Fidelity"] + by_class["Relay"] + by_class["Extension"]
    )
    if non_below != [10, 11, 15, 16]:
        raise AssertionError(f"non-Below set lock failed: {non_below}")

    # Cross-source class agreement
    by_carrier = {
        k: sorted(int(x) for x in g["family_id"].tolist())
        for k, g in carrier_df.groupby("fidelity_class")
    }
    for key in ("Fidelity", "Relay", "Extension", "Below_detection"):
        if sorted(by_class[key]) != sorted(by_carrier[key]):
            raise AssertionError(f"Figure 3 class disagreement for {key}")

    def supported_comp(fid: int) -> str:
        row = classes.loc[classes["family_id"] == fid].iloc[0]
        return str(row["snRNA_supported_compartment"])

    f10_comp = supported_comp(10)
    f11_comp = supported_comp(11)
    f15_comp = supported_comp(15)
    f16_comp = supported_comp(16)
    if f10_comp != "fibroblast":
        raise AssertionError(f"F10 compartment lock failed: {f10_comp}")
    if f11_comp != "adipocyte":
        raise AssertionError(f"F11 compartment lock failed: {f11_comp}")
    if f15_comp != "endothelial":
        raise AssertionError(f"F15 compartment lock failed: {f15_comp}")
    if set(f16_comp.split(",")) != {"atrialCM", "fibroblast"}:
        raise AssertionError(f"F16 compartment lock failed: {f16_comp}")

    # Carrier table agreement for supported compartments
    c15 = str(carrier_df.loc[carrier_df.family_id == 15, "sig_compartments"].iloc[0])
    c16 = str(carrier_df.loc[carrier_df.family_id == 16, "sig_compartments"].iloc[0])
    c10 = str(carrier_df.loc[carrier_df.family_id == 10, "sig_compartments"].iloc[0])
    c11 = str(carrier_df.loc[carrier_df.family_id == 11, "sig_compartments"].iloc[0])
    if c15 != "endothelial" or c10 != "fibroblast" or c11 != "adipocyte":
        raise AssertionError("carrier fidelity compartment disagreement")
    if set(c16.split(",")) != {"atrialCM", "fibroblast"}:
        raise AssertionError("carrier F16 compartment disagreement")

    # Panel C/D displayed entities
    csum = pd.read_csv(fig3c_sum, sep="\t")
    dsum = pd.read_csv(fig3d_sum, sep="\t")
    panel_c_fams = sorted(set(int(x) for x in csum["family_id"].tolist()))
    panel_d_disp = dsum.loc[dsum["displayed"] == True]
    panel_d_fams = sorted(set(int(x) for x in panel_d_disp["family_id"].tolist()))
    if panel_c_fams != [10, 11]:
        raise AssertionError(f"Figure 3C displayed families unexpected: {panel_c_fams}")
    if panel_d_fams != [15, 16]:
        raise AssertionError(f"Figure 3D displayed families unexpected: {panel_d_fams}")
    if (dsum["family_id"] == 4).any():
        raise AssertionError("F4 BCAA unexpectedly present in Figure 3D source")

    labels = pd.read_csv(fam_lab, sep="\t")
    lab = {
        int(r.family_id): str(r.family_label).replace("signalling", "signaling")
        for _, r in labels.iterrows()
    }
    # Publication-facing short labels required by phase brief
    pub = {
        10: "GAG disease / TGF-β",
        11: "BMP signaling",
        15: "NTRK2 / TrkB",
        16: "Keratan sulfate",
        2: "Mito Ca²⁺ transport",
        9: "HS-GAG / bone / GPCR",
    }

    # Figure 4 locks
    d4 = pd.read_csv(design4, sep="\t").iloc[0]
    if int(d4["primary_eligible_af_laa"]) != 12:
        raise AssertionError("Figure 4 primary N lock failed")
    if bool(d4["uses_acairel"]):
        raise AssertionError("Figure 4 unexpectedly uses aCAIrel")
    if str(d4["lapv_role"]) != "regional_sensitivity_only":
        raise AssertionError("Figure 4 LAPV role lock failed")
    if str(d4["input_layer"]) != "family_activity_score":
        raise AssertionError("Figure 4 input layer lock failed")

    chron_df = pd.read_csv(chron, sep="\t")
    chron_by = {
        k: sorted(f"F{int(x)}" for x in g["family_id"].tolist())
        for k, g in chron_df.groupby("chronicity_class")
    }
    strong = ["F2", "F9"]
    weak = ["F3", "F7", "F8", "F10", "F13", "F15"]
    stable = ["F1", "F4", "F5", "F6", "F11", "F12", "F14", "F16"]
    if chron_by.get("Strong_trend") != strong:
        raise AssertionError(f"Strong-trend lock failed: {chron_by.get('Strong_trend')}")
    if chron_by.get("Weak_trend") != sorted(weak):
        raise AssertionError(f"Weak-trend lock failed: {chron_by.get('Weak_trend')}")
    if chron_by.get("Stable") != sorted(stable):
        raise AssertionError(f"Stable lock failed: {chron_by.get('Stable')}")
    all_fams = set(chron_by["Strong_trend"]) | set(chron_by["Weak_trend"]) | set(chron_by["Stable"])
    if len(all_fams) != 16:
        raise AssertionError("Figure 4 classes do not cover 16 families")
    if len(chron_by["Strong_trend"]) + len(chron_by["Weak_trend"]) + len(chron_by["Stable"]) != 16:
        raise AssertionError("Figure 4 class overlap or count error")

    laa_df = pd.read_csv(laa, sep="\t").set_index("family_id")
    lapv_df = pd.read_csv(lapv, sep="\t").set_index("family_id")
    fig4d_df = pd.read_csv(fig4d, sep="\t")
    disp4 = [int(x) for x in fig4d_df.sort_values("facet_rank")["family_id"].tolist()]
    if disp4 != [9, 2, 15]:
        raise AssertionError(f"Figure 4D displayed entities unexpected: {disp4}")

    def rho3(fid: int, table, col: str) -> float:
        return float(table.loc[fid, col])

    r9_laa = rho3(9, laa_df, "rho_spearman")
    r2_laa = rho3(2, laa_df, "rho_spearman")
    r15_laa = rho3(15, laa_df, "rho_spearman")
    r9_lapv = rho3(9, lapv_df, "rho_spearman")
    r2_lapv = rho3(2, lapv_df, "rho_spearman")
    r15_lapv = rho3(15, lapv_df, "rho_spearman")
    # expected approximate locks
    if abs(r9_laa - (-0.420)) > 0.001 or abs(r9_lapv - (-0.161)) > 0.001:
        raise AssertionError("F9 rho lock failed")
    if abs(r2_laa - (-0.564)) > 0.001 or abs(r2_lapv - 0.196) > 0.001:
        raise AssertionError("F2 rho lock failed")
    if abs(r15_laa - (-0.385)) > 0.001 or abs(r15_lapv - (-0.186)) > 0.001:
        raise AssertionError("F15 rho lock failed")

    # regional signs from chronicity table
    sign = chron_df.set_index("family_id")["regional_sign_agreement"].to_dict()
    if bool(sign[9]) is not True or bool(sign[2]) is not False or bool(sign[15]) is not True:
        raise AssertionError("regional sign locks failed")

    def uminus(x: float, nd: int = 3) -> str:
        return f"{x:.{nd}f}".replace("-", "−")

    # Lock display rhos to approved three-decimal values before emitting polished prose.
    if (
        uminus(r9_laa) != "−0.420"
        or uminus(r9_lapv) != "−0.161"
        or uminus(r2_laa) != "−0.564"
        or uminus(r2_lapv) != "0.196"
        or uminus(r15_laa) != "−0.385"
        or uminus(r15_lapv) != "−0.186"
    ):
        raise AssertionError("Figure 4 display rho lock failed during Phase 4B-P emit")

    fig3 = (
        "LA single-nucleus data resolve selective family localization",
        [
            (
                "To place bulk remodeling in a cellular context, we analyzed LA single-nucleus "
                "RNA-seq from GSE255612 (18 AF and 16 non-fibrillating control donors) using "
                "donor-aware compartment pseudobulks across nine LA compartments: atrial "
                "cardiomyocytes, fibroblasts, endothelial cells, macrophages, lymphocytes, "
                "pericytes, vascular smooth muscle cells, adipocytes, and neuronal cells "
                "(Figure 3A). Frozen NCP program and family definitions were carried forward "
                "without reselection, and localization was evaluated by nominal permutation "
                "testing. Because this resource contains LA samples only, this layer was used to "
                "localize AF-associated family shifts within LA compartments rather than to retest "
                "the paired-chamber remodeling identified in bulk data."
            ),
            (
                "At the family level (Figure 3B–C), localization was selective. Final "
                "carrier-relative localization classes comprised one Fidelity family "
                "(F10 GAG disease / TGF-β), two Relay families (F11 BMP signaling and "
                "F15 NTRK2 / TrkB), one Extension family (F16 Keratan sulfate), and "
                "12 Below-detection families. F10 showed nominal Fidelity localization in "
                "fibroblasts; F11 showed nominal Relay localization in adipocytes; F15 showed "
                "nominal Relay localization in endothelial cells only; and F16 showed Extension "
                "involving fibroblasts and atrial cardiomyocytes. These localization calls are "
                "nominal and are not FDR-controlled discovery. Below-detection indicates that no "
                "included LA compartment met the preset support criteria and does not imply "
                "absence of biological remodeling."
            ),
            (
                "Figure 3C summarizes the resulting class architecture, whereas Figure 3D provides "
                "donor-level examples for F15 and F16, illustrating Relay and Extension without "
                "expanding the non-Below set. Together, panels B–D show that LA compartment "
                "localization was selective: four of 16 families were non-Below-detection, spanning "
                "Fidelity, Relay, and Extension patterns. This layer therefore supplies cellular "
                "context for selected families rather than bi-chamber replication or causal "
                "validation of the bulk remodeling results."
            ),
        ],
    )

    fig4 = (
        "Cross-sectional duration analysis identifies selective chronicity trends",
        [
            (
                "We next asked whether frozen family activity scores tracked AF duration in "
                "GSE41177. The primary analysis included left-atrial-appendage (LAA) samples from "
                "12 AF donors with finite numeric month-coded duration (Figure 4A). Because this "
                "was not a paired LA–RA duration design, the analysis used frozen family activity "
                "scores rather than aCAIrel. Left atrial–pulmonary vein junction (LAPV) samples "
                "provided regional sensitivity only, and sinus-rhythm samples were not included in "
                "the final main-figure duration panels. Duration associations were interpreted as "
                "cross-sectional chronicity context."
            ),
            (
                "Chronicity classes were assigned from the full-precision absolute LAA Spearman "
                "correlation, |ρ_LAA|; displayed rounding, permutation p values, bootstrap "
                "intervals, and LAPV estimates did not determine class (Figure 4B–C). The "
                "Strong-trend group comprised F2 Mito Ca²⁺ transport and "
                "F9 HS-GAG / bone / GPCR. Six families—F3, F7, F8, F10, F13, and F15—were "
                "Weak-trend, and the remaining eight were Stable. Strong-, Weak-, and Stable-trend "
                "labels summarize trend magnitude rather than statistical significance."
            ),
            (
                "Figure 4D displays three representative families. F9 showed a Strong negative LAA "
                "duration trend (Spearman ρ = −0.420) and a same-sign negative LAPV sensitivity "
                "estimate (ρ = −0.161). F2 showed a Strong negative LAA trend (ρ = −0.564) but an "
                "opposite-sign positive LAPV estimate (ρ = 0.196), illustrating regional sign "
                "discordance. F15 showed a Weak-trend negative LAA association (ρ = −0.385) and a "
                "same-sign negative LAPV estimate (ρ = −0.186). Regional sign is a sensitivity "
                "descriptor only and is not interpreted as replication or failure."
            ),
            (
                "Overall, duration associations were selective rather than global. Trend magnitude "
                "and regional sign agreement are distinct dimensions, and regional direction was "
                "not uniform across families. The analysis supplies cross-sectional chronicity "
                "context only and does not establish temporal progression, disease evolution, "
                "therapeutic response, or recurrence risk."
            ),
        ],
    )

    return [fig3, fig4]



def build_results_fig5_fig6() -> list[tuple[str, list[str]]]:
    """Load final Fig5/Fig6 sources, assert locks, return Figure 5–6 Results sections."""
    import pandas as pd

    def fids(xs):
        return sorted((f"F{int(x)}" for x in xs), key=lambda s: int(s[1:]))

    tri_path = ROOT / "results/ncp_v3_singscore_reanalysis/fig5b_cross_layer_triangulation_source.tsv"
    contract_path = ROOT / "results/ncp_v3_singscore_reanalysis/fig5a_triangulation_contract_source.tsv"
    geomx_sum = ROOT / "results/ncp_v3_singscore_reanalysis/fig5d_geomx_roi_spatial_evidence_summary.tsv"
    geomx_full = ROOT / "Project_P5/GSE261363/roi_scores_v3.tsv"
    comp5 = ROOT / "results/ncp_v3_singscore_reanalysis/figure5_composite_manifest.tsv"
    prov_path = ROOT / "results/provenance/gwas_gcst90624411/gwas_gcst90624411_provenance.tsv"
    a_path = ROOT / "results/fig6a_anchor_coloc_calibration_source.tsv"
    b_path = ROOT / "results/fig6b_posterior_lead_distance_source.tsv"
    c_path = ROOT / "results/fig6c_family_af_overlap_source.tsv"
    csum_path = ROOT / "results/fig6c_family_af_overlap_summary.tsv"
    d_path = ROOT / "results/fig6d_non_nested_traceability_source.tsv"
    dsum_path = ROOT / "results/fig6d_non_nested_traceability_summary.tsv"
    comp6 = ROOT / "results/figure6_composite_manifest.tsv"

    tri = pd.read_csv(tri_path, sep="\t")
    if len(tri) != 16:
        raise AssertionError("Figure 5 triangulation must cover 16 families")
    by = {int(k): fids(g.family_id) for k, g in tri.groupby("support_count")}
    if by.get(3) != ["F10", "F15"]:
        raise AssertionError(f"3/3 lock failed: {by.get(3)}")
    if by.get(2) != ["F9"]:
        raise AssertionError(f"2/3 lock failed: {by.get(2)}")
    if by.get(1) != ["F1", "F2", "F3", "F7", "F8", "F11", "F12", "F13", "F14", "F16"]:
        raise AssertionError(f"1/3 lock failed: {by.get(1)}")
    if by.get(0) != ["F4", "F5", "F6"]:
        raise AssertionError(f"0/3 lock failed: {by.get(0)}")
    if sum(len(v) for v in by.values()) != 16:
        raise AssertionError("support-count coverage failed")
    path_a = fids(tri.loc[tri["Path_A"] == True, "family_id"])
    if path_a != ["F10", "F15"]:
        raise AssertionError(f"Path A lock failed: {path_a}")

    one = tri[tri.support_count == 1].copy()

    def gates(r):
        g = []
        if r.bulk_consistency == "High":
            g.append("bulk")
        if r.snRNA_class != "Below_detection":
            g.append("snRNA")
        if r.duration_class != "Stable":
            g.append("duration")
        return g

    one["gates"] = one.apply(gates, axis=1)
    bulk_only = fids(one[one.gates.apply(lambda g: g == ["bulk"])].family_id)
    sn_only = fids(one[one.gates.apply(lambda g: g == ["snRNA"])].family_id)
    dur_only = fids(one[one.gates.apply(lambda g: g == ["duration"])].family_id)
    if bulk_only != ["F1", "F12", "F14"]:
        raise AssertionError(f"bulk-only 1/3 failed: {bulk_only}")
    if sn_only != ["F11", "F16"]:
        raise AssertionError(f"snRNA-only 1/3 failed: {sn_only}")
    if dur_only != ["F2", "F3", "F7", "F8", "F13"]:
        raise AssertionError(f"duration-only 1/3 failed: {dur_only}")

    contract = {
        str(r["metric"]): str(r["value"])
        for _, r in pd.read_csv(contract_path, sep="\t").iterrows()
    }
    if contract.get("n_Path_A_gates") != "3":
        raise AssertionError("Path A gate count lock failed")
    for key in ("hca_is_gate", "dominant_state_is_gate", "geomx_is_gate", "genetics_is_gate", "priorityscore_opened"):
        if contract.get(key) not in {"False", "false", "0"}:
            raise AssertionError(f"{key} must be False")
    if "Figure_1D" not in contract.get("canonical_order_source", ""):
        raise AssertionError("canonical Figure 1D order lock failed")

    panels5 = set(pd.read_csv(comp5, sep="\t")["panel"].astype(str))
    if panels5 != {"A", "B", "C", "D"}:
        raise AssertionError(f"Figure 5 panel structure not A–D: {panels5}")

    gx = pd.read_csv(geomx_sum, sep="\t")
    f10 = gx.loc[gx["name"] == "F10"].iloc[0]
    f15 = gx.loc[gx["name"] == "F15"].iloc[0]
    f9 = gx.loc[gx["name"] == "F9"].iloc[0]
    if abs(float(f10["rho"]) - 0.439130) > 1e-5 or abs(float(f10["p"]) - 0.031804) > 1e-5:
        raise AssertionError("GeoMx F10 stats mismatch")
    if abs(float(f15["rho"]) - 0.431304) > 1e-5 or abs(float(f15["p"]) - 0.035350) > 1e-5:
        raise AssertionError("GeoMx F15 stats mismatch")
    if abs(float(f9["delta"]) - (-0.033333)) > 1e-5 or abs(float(f9["p"]) - 0.048130) > 1e-5:
        raise AssertionError("GeoMx F9 stats mismatch")
    if int(f10["n_rois"]) != 24 or int(f10["n_patients"]) != 4:
        raise AssertionError("GeoMx eligible ROI/patient lock failed")

    full = pd.read_csv(geomx_full, sep="\t")
    n_full = int(full["gsm"].nunique())
    n_elig = int(full.loc[full["eligible_for_5d"] == True, "gsm"].nunique())
    if n_full != 48 or n_elig != 24:
        raise AssertionError(f"GeoMx full/eligible ROI lock failed: {n_full}/{n_elig}")

    # Figure 6
    prov = {
        str(r["key"]): str(r["value"])
        for _, r in pd.read_csv(prov_path, sep="\t").iterrows()
    }
    if prov.get("accession") != "GCST90624411":
        raise AssertionError("GWAS accession lock failed")
    if int(prov["locked_af_cases"]) != 252438 or int(prov["locked_controls"]) != 1959739:
        raise AssertionError("GWAS sample-size lock failed")
    if str(prov.get("coloc_rerun_required", "")).upper() != "FALSE":
        raise AssertionError("coloc_rerun_required must be FALSE")
    if "Yuan" not in prov.get("publication", ""):
        raise AssertionError("Yuan publication attribution lock failed")

    a = pd.read_csv(a_path, sep="\t")
    if len(a) != 12:
        raise AssertionError("Figure 6A must have 12 anchors")
    if int((a["anchor_type"] == "framework").sum()) != 10:
        raise AssertionError("framework anchor count lock failed")
    calib = sorted(a.loc[a["anchor_type"] == "calibration", "gene"].astype(str))
    if calib != ["KCNN3", "PITX2"]:
        raise AssertionError(f"calibration controls lock failed: {calib}")
    if int(a["local_AF_GWAS_support"].sum()) != 7:
        raise AssertionError("GWAS 7/12 lock failed")
    if int(a["atrial_eGene_bonferroni"].sum()) != 1:
        raise AssertionError("eGene 1/12 lock failed")
    if a.loc[a["atrial_eGene_bonferroni"] == True, "gene"].tolist() != ["SMAD9"]:
        raise AssertionError("eGene gene lock failed")
    if int(a["PP4_ge_0.70"].sum()) != 0:
        raise AssertionError("PP4 0/12 lock failed")
    if int(a["PP3_dominant"].sum()) != 12:
        raise AssertionError("PP3 12/12 lock failed")

    b = pd.read_csv(b_path, sep="\t")
    if len(b) != 12 or b["delta_mb"].isna().any():
        raise AssertionError("Figure 6B lead-distance completeness failed")
    if int(b["pp3_dominant"].sum()) != 12 or int(b["pp4_support"].sum()) != 0:
        raise AssertionError("Figure 6B PP3/PP4 lock failed")
    dmin, dmax = float(b["delta_mb"].min()), float(b["delta_mb"].max())
    if abs(dmin - 0.459018) > 1e-6 or abs(dmax - 1.923937) > 1e-5:
        raise AssertionError(f"lead-distance range lock failed: {dmin}-{dmax}")
    kcnn3 = float(b.loc[b["gene"] == "KCNN3", "delta_mb"].iloc[0])
    if abs(kcnn3 - 0.459018) > 1e-6:
        raise AssertionError("KCNN3 distance lock failed")

    c = pd.read_csv(c_path, sep="\t")
    csum = pd.read_csv(csum_path, sep="\t").iloc[0]
    if int(csum["n_tests"]) != 64 or int(csum["n_overlap_events"]) != 5:
        raise AssertionError("Figure 6C test/nonzero lock failed")
    if int(csum["n_fdr_significant"]) != 0 or int(csum["universe_n"]) != 14635:
        raise AssertionError("Figure 6C FDR/universe lock failed")
    nz = c[c["n_overlap"] > 0][
        ["family_id", "external_set_id", "overlap_genes"]
    ].sort_values(["family_id", "external_set_id"])
    expected_nz = {
        (3, "yuan2025_cross_population", "CAV1"),
        (3, "roselli2018_loci", "CAV1"),
        (13, "yuan2025_cross_population", "CAV1"),
        (13, "roselli2018_loci", "CAV1"),
        (9, "selewa2023_regulatory", "HEY2"),
    }
    got_nz = {
        (int(r.family_id), str(r.external_set_id), str(r.overlap_genes))
        for _, r in nz.iterrows()
    }
    if got_nz != expected_nz:
        raise AssertionError(f"Figure 6C nonzero cells mismatch: {got_nz}")

    dsum = pd.read_csv(dsum_path, sep="\t").set_index("dimension")
    expected_d = {
        "local_AF_GWAS": 7,
        "atrial_eGene": 1,
        "same_variant_coloc": 0,
        "common_risk": 2,
        "regulatory": 2,
        "rare_functional": 0,
    }
    for dim, n in expected_d.items():
        if int(dsum.loc[dim, "present"]) != n or int(dsum.loc[dim, "n"]) != 12:
            raise AssertionError(f"Figure 6D {dim} lock failed")

    panels6 = {
        str(x).replace("panel_", "").replace("_svg", "")
        for x in pd.read_csv(comp6, sep="\t")["artifact"]
        if str(x).startswith("panel_")
    }
    if panels6 != {"A", "B", "C", "D"}:
        # tolerate panel_A_svg style already handled
        if not {"A", "B", "C", "D"}.issubset(
            {
                str(x).split("_")[1]
                for x in pd.read_csv(comp6, sep="\t")["artifact"]
                if str(x).startswith("panel_")
            }
        ):
            raise AssertionError(f"Figure 6 panel structure not A–D: {panels6}")

    fig5 = (
        "Rule-based cross-layer triangulation prioritizes two NCP families",
        [
            (
                "After evaluating the paired-bulk, LA single-nucleus, and duration layers "
                "separately, we placed their corrected family-level summaries on the frozen "
                "16-family NCP topology in the canonical Figure 1D order, without fitting a merged "
                "latent model (Figure 5A–B). The healthy LA carrier prior from the Human Cell Atlas "
                "and the dominant bulk rewiring state were retained as contextual annotations. "
                "Path A used exactly three binary gates: High bulk cross-cohort consistency, a "
                "non-Below-detection LA single-nucleus class, and a non-Stable LAA duration class. "
                "The contextual annotations, GeoMx, and genetics did not contribute to the support "
                "count."
            ),
            (
                "Cross-layer support was selective (Figure 5B–C). F10 GAG disease / TGF-β and "
                "F15 NTRK2 / TrkB satisfied all three gates and were assigned Path A. F9 HS-GAG / "
                "bone / GPCR satisfied the bulk and duration gates but remained Below-detection in "
                "the LA single-nucleus layer, yielding 2/3 support. Ten families satisfied one "
                "gate: F1, F12, and F14 through bulk consistency; F11 and F16 through LA "
                "single-nucleus localization; and F2, F3, F7, F8, and F13 through duration. "
                "F4, F5, and F6 satisfied none of the three gates. The resulting architecture "
                "therefore comprised two 3/3, one 2/3, ten 1/3, and three 0/3 families. Path A is "
                "a rule-based triage label across non-exchangeable evidence layers, not a joint "
                "statistical test or causal-validation claim."
            ),
            (
                "GeoMx provided independent spatial context (Figure 5D). The full GSE261363 "
                "resource contained 48 regions of interest (ROIs) from four patients; the final "
                "panel used 24 eligible cardiac atrial-appendage DNA-positive ROIs. Across the "
                "eligible ROIs, F10 family-signature scores showed a nominal positive correlation "
                "with fibroblast-panel scores (Spearman ρ = 0.439, nominal p = 0.0318), and F15 "
                "scores showed a nominal positive correlation with macrophage-panel scores "
                "(ρ = 0.431, nominal p = 0.0354). F9 showed a nominal ROI-level AF-minus-SR "
                "difference (Δ = −0.033, Welch p = 0.0481). Multiple ROIs arose from each patient, "
                "and no patient-clustered or mixed-effects inference was applied. These nominal "
                "ROI-level summaries provide spatial context only and did not modify the Path A "
                "classification."
            ),
        ],
    )

    fig6 = (
        "Genetic follow-up reveals calibration-limited, non-nested traceability",
        [
            (
                "We next evaluated 10 framework anchors together with PITX2 and KCNN3 AF "
                "calibration controls using the Yuan et al. 2025 cross-population AF GWAS "
                "(GWAS Catalog accession GCST90624411; corrected publication sample size: "
                "252,438 AF cases and 1,959,739 controls) and GTEx left-atrial-appendage cis-eQTL "
                "data (Figure 6A). Seven of the 12 anchor windows contained local genome-wide AF "
                "associations. Only SMAD9 met the corrected atrial eGene criterion, and no anchor "
                "met the prespecified same-variant colocalization threshold of PP4 ≥ 0.70. SMAD9 "
                "therefore illustrates atrial eGene evidence without GWAS–eQTL colocalization. "
                "The same lack of colocalization support at PITX2 and KCNN3 helped calibrate the "
                "interpretation of the overall non-support."
            ),
            (
                "Continuous posterior diagnostics reinforced this interpretation (Figure 6B). "
                "All 12 anchors were PP3-dominant, whereas PP4 remained below the support "
                "threshold. The frozen lead AF GWAS and atrial eQTL variants were separated by "
                "0.46–1.92 Mb across the tested anchors. Lead separation was a descriptive "
                "companion diagnostic and was not used to calculate the colocalization posterior. "
                "Together, the high-PP3, low-PP4 pattern is consistent with distinct local "
                "association signals rather than supported same-variant bulk-LAA total-expression "
                "mediation."
            ),
            (
                "At the family level, we compared the frozen leading-edge genes of 16 NCP "
                "families with four prespecified external AF gene sets, yielding 64 family–set "
                "tests against a background of 14,635 NCP-eligible genes (Figure 6C). Only five "
                "family–set cells contained any member overlap, involving two unique genes across "
                "three families. Four cells reflected CAV1 membership in both F3 Adipo-lipid axis "
                "and F13 Amino-acid metabolism across two common-risk or prioritized sets; the "
                "fifth reflected HEY2 membership in F9 HS-GAG / bone / GPCR and a regulatory set. "
                "No family–set enrichment survived BH-FDR correction. The repeated CAV1 cells "
                "reflect shared set membership rather than independent genetic observations, so "
                "these overlaps are traceability cues rather than family-level enrichment."
            ),
            (
                "Finally, the non-nested traceability matrix summarized independent evidence "
                "dimensions without converting them into an ordinal tier or composite score "
                "(Figure 6D). Local AF GWAS, atrial eGene, and same-variant colocalization flags "
                "were present for 7/12, 1/12, and 0/12 anchors, respectively. Membership in risk "
                "or prioritized gene sets and regulatory annotations was present for 2/12 anchors "
                "in each category, whereas rare or functional membership was present for 0/12. "
                "These categories are not exchangeable and cannot be ordered or summed. Genetics "
                "did not enter Path A and is interpreted here as calibration-limited traceability "
                "for follow-up rather than causal validation or evidence excluding untested genetic "
                "or cis-regulatory mechanisms."
            ),
        ],
    )

    return [fig5, fig6]


SECTIONS_AFTER_FIG6: list[tuple[str, list[str]]] = [
    (
        'Specificity, robustness, and methodological validity',
        [
                (
                    'Because the dual-axis interpretation rests on a small set of '
                    'pre-frozen programs, we examined its specificity and methodological '
                    'validity through a four-component falsification supplement and two '
                    'methodology supplements (Supplementary Figs S8, S10, S11). Random '
                    'size-matched gene-set permutations (1,000 per program) produced a '
                    'non-trivial empirical p for the strongest single program (ECM '
                    'proteoglycans empirical p = 0.10), and an axis-membership '
                    'permutation (1,000 random 3:3 partitions of the six frozen support '
                    'programs) yielded |gSTRUCT| empirical p = 0.33; both axes pooled '
                    'negatively under any 3:3 partition. Six pre-registered '
                    'negative-control pathways (HALLMARK_INFLAMMATORY_RESPONSE, '
                    'KRAS_SIGNALING_UP, E2F_TARGETS for axis-specificity; '
                    'REACTOME_OLFACTORY_SIGNALING_PATHWAY, KERATINIZATION, '
                    'SENSORY_PERCEPTION_OF_TASTE for biological unrelatedness) sometimes '
                    'matched a single statistic in magnitude (for example, '
                    'REACTOME_KERATINIZATION pooled g = −1.05) but were accompanied by '
                    'high cross-cohort heterogeneity (I² = 78 %) and a much weaker '
                    'LAA-duration ρ (|ρ| ≤ 0.29). Specificity therefore arises from the '
                    'joint combination of pooled effect direction, low cross-cohort '
                    'heterogeneity, and a coherent LAA-duration association rather than '
                    'from any single statistic. A leave-one-program-out (LOPO) test '
                    'preserved the direction of both axes across all six removals (6 / '
                    '6).'
                ),
                (
                    "The framework's pipeline was further validated through "
                    'internal-consistency checks (Supplementary Fig. S10): NCP '
                    'self-recovery on the healthy bulk reference produced a mean aCAIrel '
                    '≈ 1 for all thirty NCP pathway programs; donor-level direction '
                    'stability was 8 / 8 healthy donors with aCAI > 0 for every panel '
                    'program; cross-method per-donor correlation between singscore and a '
                    'mean-z scoring ranged from 0.52 to 0.93; and axis-pooled g was '
                    'unchanged when singscore was replaced with mean-z (structural g = '
                    '−0.57 / −0.79; both with I² = 0 %). External-style validation '
                    '(Supplementary Fig. S11) showed that the pipeline recovered known '
                    'axis effects with at least 92 % detection at a true effect of 0.3 in '
                    'aCAIrel units, that the cohort-donor bootstrap of the dual-axis '
                    'pooled g produced 95 % intervals consistent with and tighter than '
                    'the Hartung–Knapp intervals from the axis meta-analysis '
                    '(results/figure6_meta_axis.tsv), and that a classical limma '
                    'per-chamber AF-vs-SR analysis agreed in direction with the framework '
                    'in 25 / 30 program × cohort cells (83 %; 100 % in GSE115574), with '
                    'discordant cells corresponding to near-zero classical effects. A '
                    'reduced 5-gene-per-axis panel reproduced the full-axis donor-level '
                    'composite (Spearman ρ = 0.55–0.90) and the cohort-level SR → AF '
                    'centroid direction; the reduced metabolic panel did not, however, '
                    'replace OXPHOS for chronicity and is therefore reported as a '
                    'framework-feasibility check rather than a clinical classifier '
                    '(Supplementary Fig. S9).'
                ),
                (
                    'Sex-stratified sensitivity was constrained by data availability: '
                    'GSE115574 lacked per-donor sex in the GEO record, GSE128188 '
                    'contained only male donors, and only GSE79768 contained both sexes '
                    'with sufficient AF and SR representation. Within GSE79768 the '
                    "structural axis was directionally consistent in both sexes (Hedges' "
                    'g = −0.76 in male donors, −0.27 in female donors), whereas the '
                    'metabolic axis was too small to call. A cross-cohort sex-stratified '
                    'meta-analysis was therefore underpowered (k ≤ 2 per stratum) and is '
                    'not presented as a main-text result.'
                ),
        ],
    ),
]

def main() -> None:
    doc = Document()
    apply_base_style(doc)

    doc.add_heading("Results", level=1)

    sections = (
        build_results_preamble_fig1_fig2()
        + build_results_fig3_fig4()
        + build_results_fig5_fig6()
        + SECTIONS_AFTER_FIG6
    )
    for header, paragraphs in sections:
        doc.add_heading(header, level=2)
        for p in paragraphs:
            para = doc.add_paragraph(p)
            para.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY

    note = doc.add_paragraph()
    note.paragraph_format.space_before = Pt(18)
    run = note.add_run(
        "End-of-Results note (not for typesetter): all numerical "
        "claims in this section are reproducible from the "
        "tab-separated deliverable files released with the manuscript "
        "(see Code availability) and are summarised per figure in the "
        "Methods text and in the multiple-testing-correction summary "
        "(Supplementary Table S2)."
    )
    run.italic = True
    run.font.size = Pt(9)
    run.font.color.rgb = RGBColor(0x66, 0x66, 0x66)

    doc.save(OUT)
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
