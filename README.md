# Drought-induced dispersal limitation fragments rhizosphere metacommunity assembly while preserving plant growth-promoting functional capacity in maize

**Olanrewaju OS, Ayangbenro AS, Akanmu AO, Bezuidenhout CC, and Babalola OO**

**

---

## Overview

This repository contains all analysis scripts supporting the manuscript.
---

## Repository Structure

    maize-rhizosphere-drought-assembly/
    |
    |-- Analysis-pipeline.R          # Main analysis pipeline 
    |-- Theme_constants.R              # Global color grammar and typography
    |
    |-- scripts/
    |   └-- CHPC/
    |       |-- run_icamp.pbs          # PBS job script for iCAMP (CHPC cluster)
    |       |-- run_icamp.R            # iCAMP community assembly analysis
    |       |-- fix_rcbray.pbs         # PBS job script for RC_bray correction
    |       └-- fix_rcbray.R           # RC_bray correction script
    |
    |-- data/
    |   |-- Morphological_characterization.csv          # Corrected morphological data
    |   └-- Morphological_characterization_original.csv # Original uncorrected data
    |
    |-- .gitignore
    └-- README.md

---

## Requirements

### R version
R v4.5.2 (2025-10-31)

### Key packages

Bioconductor:
- phyloseq v1.54.0
- ANCOMBC v2.10.1
- ALDEx2 v1.40.0

CRAN:
- vegan v2.7.2
- lme4 v1.1.37
- lmerTest v3.1.3
- ranger v0.17.0
- pROC v1.19.0.1
- betapart v1.6.1
- picante v1.8.2
- NST v3.1.10
- igraph v2.2.2
- ggVennDiagram v1.5.4

GitHub:
- SpiecEasi v1.99.0
- Install with: remotes::install_github("zdk123/SpiecEasi")

### External tools
- PICRUSt2 v2.5.2 — conda environment (bioconda channel)
- FAPROTAX v1.2.11 — Python script (collapse_table.py)
- iCAMP v1.3.4 — Run on CHPC Lengau (see scripts/CHPC/)
- iCAMP v1.5.12 — used for local Mac runs

---

## iCAMP on HPC

Community assembly mechanism analysis (iCAMP) was run on the Centre for High Performance Computing (CHPC) Lengau cluster due to computational requirements (~30 min with 20 parallel workers, 90 GB RAM).

Before submitting PBS scripts, replace all placeholders in the .pbs files:

| Placeholder | Replace with |
|-------------|-------------|
| YOUR_PROJECT_CODE | Your HPC project allocation code |
| YOUR_USERNAME | Your HPC username |
| YOUR_EMAIL | Your email |
| YOUR_QUEUE | Your HPC queue |
| YOUR_NCPUS | Number of CPUs to request |
| YOUR_MEM | Memory to request (e.g. 100gb) |
| YOUR_WALLTIME | Walltime limit (e.g. 24:00:00) |

---

## Data Availability

Raw 16S rRNA amplicon sequencing data are deposited in the NCBI Sequence Read Archive (SRA):

BioProject: PRJNA1422815
URL: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1422815

---

## Morphological Data Corrections

Four morphological measurements with biologically implausible values were identified and replaced with NA prior to statistical analysis:

| Sample | Trait | Original | Corrected |
|--------|-------|----------|-----------|
| Watered/Day0/Susceptible/Rep3 | Shoot_length | 0.58 cm | NA |
| Watered/Day2/Resistant/Rep2 | Shoot_fresh_weight | 11.64 g | NA |
| Watered/Day2/Resistant/Rep3 | Shoot_fresh_weight | 6.08 g | NA |
| Watered/Day5/Resistant/Rep3 | Shoot_length | 0.00 cm | NA |

Both the corrected (Morphological_characterization.csv) and original (Morphological_characterization_original.csv) files are provided for full transparency.

---

## Reproducibility

All stochastic analyses use random seed 42. R session information including all package versions is saved automatically to tables/session_info.txt during each pipeline run.

---

## Citation

Olanrewaju OS, Ayangbenro AS, Akanmu AO, Bezuidenhout CC, and Babalola OO (2025). Drought-induced dispersal limitation fragments rhizosphere metacommunity assembly while preserving plant growth-promoting functional capacity in maize.

---

## License

Scripts are released under the MIT License.
Data files are released under CC BY 4.0.

---

## Contact

Author: Dr. Oluwaseyi Samuel Olanrewaju
Olusam777@gmail.com

