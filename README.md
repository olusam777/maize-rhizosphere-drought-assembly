# Drought-induced dispersal limitation fragments rhizosphere metacommunity assembly while preserving plant growth-promoting functional capacity in maize

**Olanrewaju OS, Ayangbenro AS, Akanmu AO, Bezuidenhout CC, and Babalola OO**

*Manuscript under review — Microbiome (2025)*

---

## Overview

This repository contains all analysis scripts and data files supporting the manuscript. The study examines how progressive drought stress restructures rhizosphere bacterial community assembly mechanisms and functional capacity in two maize genotypes (drought-resistant and drought-susceptible) across a six-day greenhouse experiment (n = 108 samples).

**Central finding:** Drought doubles dispersal limitation (16.0% to 29.5%) and collapses homogenizing dispersal (22.4% to 7.0%), fragmenting rhizosphere metacommunity assembly, while plant growth-promoting functional capacity is maintained through functional redundancy.

---

## Repository Structure

    maize-rhizosphere-drought-assembly/
    |
    |-- drought_pipeline_v2.r          # Main analysis pipeline (~6,500 lines, 8 analytical pillars)
    |-- figure_assembly.r              # Figure compilation script
    |-- theme_constants.r              # Global color grammar and typography
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

## Pipeline Structure

drought_pipeline_v2.r is organized into 8 analytical pillars executed sequentially:

| Pillar | Analysis |
|--------|----------|
| 1 | Data import, quality filtering, rarefaction, alpha and beta diversity |
| 2 | Community assembly mechanisms (iCAMP, βNTI, RC_bray) |
| 3 | Differential abundance (ANCOM-BC2), functional profiling (PICRUSt2, FAPROTAX) |
| 4 | Plant phenotype-microbiome coupling (Procrustes, Mantel) |
| 5 | Predictive modelling (Random Forest) |
| 6 | Extended community composition and core microbiome |
| 7 | Co-occurrence network analysis (SPIEC-EASI, Zi-Pi keystone classification) |
| 8 | Phylogenetic signal analysis and supplementary figure assembly |

Estimated full pipeline runtime: 4-6 hours on a MacBook with 16 GB RAM. Computationally intensive steps (SPIEC-EASI x4, ANCOM-BC2 x5, Random Forest) are cached as RDS files — delete the relevant RDS file to force recomputation.

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
| YOUR_EMAIL@institution.ac.za | Your institutional email |
| YOUR_QUEUE | Your HPC queue name |
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

Olanrewaju OS, Ayangbenro AS, Akanmu AO, Bezuidenhout CC, and Babalola OO (2025). Drought-induced dispersal limitation fragments rhizosphere metacommunity assembly while preserving plant growth-promoting functional capacity in maize. Microbiome (under review).

---

## License

Scripts are released under the MIT License.
Data files are released under CC BY 4.0.

---

## Contact

Corresponding author: Dr. Oluwaseyi Samuel Olanrewaju
Unit for Environmental Sciences and Management
North-West University, Potchefstroom, South Africa
ORCID: [your ORCID here]
