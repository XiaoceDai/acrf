"""Central pathway-name abbreviation table for Fig 1 panels.

Used by fig1c, fig1d, fig1e, fig1f to make every pathway label
consistent and short enough to be readable at 5-6 pt fonts in
publication-size panels.
"""
from __future__ import annotations

# Map full MSigDB pathway name (with the HALLMARK__ / REACTOME__ prefix
# stripped) to a short, biologist-readable label.  Add new entries here
# as more pathways need to be displayed.
PATHWAY_SHORT = {
    # Mitochondrial / OXPHOS family
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION":          "OXPHOS",
    "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT":     "ETC",
    "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT":
        "Aerobic resp / ETC",
    "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE":        "TCA cycle",
    "REACTOME_COMPLEX_I_BIOGENESIS":               "Complex I biogenesis",
    "REACTOME_COMPLEX_III_ASSEMBLY":               "Complex III assembly",
    "REACTOME_COMPLEX_IV_ASSEMBLY":                "Complex IV assembly",
    "REACTOME_MITOCHONDRIAL_TRANSLATION":          "Mito translation",
    "REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT":       "Mito protein import",
    "REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION":  "Mito proteostasis",
    "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION":
        "Mito β-oxidation",
    "REACTOME_MITOCHONDRIAL_TRNA_AMINOACYLATION":  "Mito tRNA aminoacyl.",
    "REACTOME_MATURATION_OF_TCA_ENZYMES_AND_REGULATION_OF_TCA_CYCLE":
        "TCA enzyme maturation",
    "REACTOME_MITOCHONDRIAL_CALCIUM_ION_TRANSPORT":
        "Mito Ca²⁺ transport",

    # Translation family
    "REACTOME_TRANSLATION":                                      "Translation (cyto)",
    "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION":                "Eukaryotic translation elong.",
    "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION":                "Eukaryotic translation init.",
    "REACTOME_RIBOSOME_ASSOCIATED_QUALITY_CONTROL":              "Ribosome QC (RQC)",
    "REACTOME_RIBOSOME_QUALITY_CONTROL_RQC_COMPLEX_EXTRACTION":  "RQC complex extraction",
    "REACTOME_TRNA_AMINOACYLATION":                              "tRNA aminoacylation",
    "REACTOME_NONSENSE_MEDIATED_DECAY_NMD":                      "NMD",
    "REACTOME_RRNA_PROCESSING":                                  "rRNA processing",
    "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY":
        "GCN2 / amino-acid response",
    "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE":
        "SRP cotranslational targeting",
    "REACTOME_PROTEIN_LOCALIZATION":                             "Protein localization",
    "REACTOME_INFLUENZA_VIRAL_RNA_TRANSCRIPTION_AND_REPLICATION":
        "Influenza viral RNA",
    "REACTOME_INFLUENZA_INFECTION":                              "Influenza infection",
    "REACTOME_HOST_INTERACTIONS_OF_HIV_FACTORS":                 "HIV host interactions",
    "REACTOME_SARS_COV_1_MODULATES_HOST_TRANSLATION_MACHINERY":
        "SARS-CoV-1 host translation",
    "REACTOME_SARS_COV_2_MODULATES_HOST_TRANSLATION_MACHINERY":
        "SARS-CoV-2 host translation",
    "REACTOME_ACTIVATION_OF_THE_MRNA_UPON_BINDING_OF_THE_CAP_BINDING_COMPLEX_AND_EIFS_AND_SUBSEQUENT_BINDING_TO_43S":
        "mRNA cap activation",
    "REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS":
        "SLIT/ROBO expression reg.",
    "REACTOME_FORMATION_OF_TUBULIN_FOLDING_INTERMEDIATES_BY_CCT_TRIC":
        "CCT/TRiC tubulin folding",

    # Lipid / fatty acid / adipose family
    "REACTOME_TRIGLYCERIDE_METABOLISM":             "TG metabolism",
    "REACTOME_TRIGLYCERIDE_CATABOLISM":             "TG catabolism",
    "HALLMARK_FATTY_ACID_METABOLISM":               "FA metabolism",
    "REACTOME_BRANCHED_CHAIN_AMINO_ACID_CATABOLISM":"BCAA catabolism",
    "HALLMARK_ADIPOGENESIS":                        "Adipogenesis",
    "HALLMARK_PEROXISOME":                          "Peroxisome",
    "REACTOME_PEROXISOMAL_PROTEIN_IMPORT":          "Peroxisomal protein import",
    "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES":
        "Amino-acid metabolism",
    "REACTOME_FATTY_ACID_METABOLISM":               "FA metabolism (Reactome)",
    "REACTOME_SELENOAMINO_ACID_METABOLISM":         "Selenoamino-acid metab.",

    # Matrix / GAG / ECM family
    "REACTOME_HEPARAN_SULFATE_HEPARIN_HS_GAG_METABOLISM":
        "HS-GAG metabolism",
    "REACTOME_HS_GAG_BIOSYNTHESIS":                 "HS-GAG biosynthesis",
    "REACTOME_HS_GAG_DEGRADATION":                  "HS-GAG degradation",
    "REACTOME_DISEASES_ASSOCIATED_WITH_GLYCOSAMINOGLYCAN_METABOLISM":
        "GAG metabolism diseases",
    "REACTOME_KERATAN_SULFATE_BIOSYNTHESIS":        "Keratan sulfate biosynth.",
    "REACTOME_CHONDROITIN_SULFATE_DERMATAN_SULFATE_METABOLISM":
        "CS / DS metabolism",
    "REACTOME_LAMININ_INTERACTIONS":                "Laminin interactions",
    "REACTOME_FORMATION_OF_THE_DYSTROPHIN_GLYCOPROTEIN_COMPLEX":
        "Dystrophin glycoprotein cx",
    "REACTOME_NON_INTEGRIN_MEMBRANE_ECM_INTERACTIONS":
        "Non-integrin membrane-ECM",
    "REACTOME_ASSEMBLY_AND_CELL_SURFACE_PRESENTATION_OF_NMDA_RECEPTORS":
        "NMDA receptor assembly",

    # Developmental / signalling family
    "REACTOME_RUNX2_REGULATES_BONE_DEVELOPMENT":    "RUNX2 / bone dev.",
    "REACTOME_RUNX2_REGULATES_OSTEOBLAST_DIFFERENTIATION":
        "RUNX2 / osteoblast diff.",
    "REACTOME_SIGNALING_BY_BMP":                    "BMP signalling",
    "REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS":    "TGF-β family signalling",
    "HALLMARK_TGF_BETA_SIGNALING":                  "TGF-β signalling (Hall)",
    "REACTOME_DOWNREGULATION_OF_SMAD2_3_SMAD4_TRANSCRIPTIONAL_ACTIVITY":
        "SMAD2/3/4 down-reg.",
    "REACTOME_SIGNALING_BY_NTRK2_TRKB":             "NTRK2 / TrkB signalling",
    "REACTOME_CLASS_B_2_SECRETIN_FAMILY_RECEPTORS": "Class B/2 Secretin receptors",
    "REACTOME_EPIGENETIC_REGULATION_OF_GENE_EXPRESSION":
        "Epigenetic gene regulation",
    "REACTOME_EPIGENETIC_REGULATION_BY_POLYCOMB":   "Polycomb epigenetic reg.",
    "REACTOME_DEGRADATION_OF_DVL":                  "Dishevelled (DVL) deg.",
    "REACTOME_AUF1_HNRNP_D0_BINDS_AND_DESTABILIZES_MRNA":
        "AUF1 mRNA destabilization",
    "REACTOME_COOPERATION_OF_PREFOLDIN_AND_TRIC_CCT_IN_ACTIN_AND_TUBULIN_FOLDING":
        "Prefoldin–TRiC cooperation",
    "REACTOME_SNRNP_ASSEMBLY":                      "snRNP assembly",
    "REACTOME_CELLULAR_RESPONSE_TO_STARVATION":     "Starvation response",
}


def short_pathway_label(full_name: str, max_chars: int = 30) -> str:
    """Return a short biologist-readable label for an MSigDB pathway name.

    Looks up the central abbreviation table first; falls back to a clean
    sentence-case version of the original with redundant prefixes
    stripped.  Always truncates to max_chars so long entries from the
    abbreviation table don't break tight figure layouts.
    """
    s = full_name
    s = s.replace("HALLMARK__HALLMARK_", "HALLMARK_")
    s = s.replace("REACTOME__REACTOME_", "REACTOME_")
    if s in PATHWAY_SHORT:
        out = PATHWAY_SHORT[s]
    else:
        out = (s.replace("HALLMARK_", "").replace("REACTOME_", "")
                .replace("_", " ").lower().capitalize())
    if len(out) > max_chars:
        out = out[: max_chars - 1] + "…"
    return out


def family_short_label(family_id: int, representative_pathway: str) -> str:
    """Compact 'F1 OXPHOS' style label for an axis-family group."""
    return f"F{family_id}  {short_pathway_label(representative_pathway, 22)}"
