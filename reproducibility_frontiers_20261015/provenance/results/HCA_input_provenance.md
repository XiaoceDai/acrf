# HCA input provenance

Existing input. The analysis path is `/home/bio/Project_P02/global_raw.h5ad`. It is not a symlink. `scripts/p_hca_pseudobulk_chamber.py` reads that path in backed mode. The same file is on disk.

Evidence linking it to Heart Cell Atlas v1. Shape 486,134 × 33,538, 14 donors (D1–D7, D11, H2–H7), six regions (AX, LA, LV, RA, RV, SP), and cell sources Harvard-Nuclei, Sanger-Nuclei, Sanger-Cells, and Sanger-CD45 match `reports/s1_schema.txt` and the published v1 description (Litviňuková et al., 2020; HCA project ad98d3cd-26fb-4ee3-99c9-8a2ab085e737). The HDF5 root attributes are empty, so the file itself does not store a version string.

Public-file byte identity. Not verified. No official checksum and no saved download URL were found. The SHA256 from this check is a hash of the local file, not a historical download checksum.

Download date. Not recovered. The filesystem modification time is 2026-04-17 09:26:52 +0800 and is not a download date.

Filter. Nuclei, LA or RA, Used = Yes, and the mapped cell types leave 77,969 cells from all 14 donors. Paired-donor counts in the existing pseudobulk table are compartment-specific (1 to 13) and were not recomputed.
