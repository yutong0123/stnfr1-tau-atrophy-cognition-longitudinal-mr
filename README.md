# stnfr1-tau-atrophy-cognition-longitudinal-mr
Pipeline for sTNFR1 to tau pathology, brain atrophy, and cognitive decline (longitudinal + MR) in aging cohorts (HRS/ADNI).
This repository contains the primary analysis code used in the manuscript (HRS/ADNI analyses and the main Mendelian randomization pipeline).

**HRS/ADNI (controlled-access data):** 
We provide the primary analysis scripts (modeling/estimation) only. Dataset-specific preprocessing/cleaning steps and variable mappings are not included due to controlled-access data restrictions and release-to-release differences.

**MR (summary statistics):**
The main MR pipeline uses `TwoSampleMR`/`ieugwasr` and queries OpenGWAS using dataset IDs specified in the script. Users must set `OPENGWAS_JWT` as an environment variable (no credentials stored).

