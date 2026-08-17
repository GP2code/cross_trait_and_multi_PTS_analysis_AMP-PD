# Cross-trait and multi-polytranscriptomic score analysis of Parkinson's disease identifies novel associations and improves prediction

`GP2 ❤️ Open Science 😍`

[![DOI](https://zenodo.org/badge/DOI/nnnnn/zenodo.nnnnn.svg)](https://doi.org/nnnnn/zenodo.nnnnn) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Last updated:** July 2026

## Summary
This is the online repository for the manuscript titled **"Cross-trait and multi-polytranscriptomic score analysis of Parkinson's disease identifies novel associations and improves prediction"**.

## Data statement

The individual level data used in this project are hosted in collaboration with the Accelerating Medicines Partnership in Parkinson’s disease, and are available via application on the website (https://amp-pd.org/register-for-amp-pd). Tier 1 data can be accessed by completing a form on the Accelerating Medicines Partnership in Parkinson’s Disease (AMP®-PD) website (https://amp-pd.org/register-for-amp-pd). Tier 2 data access requires approval and a Data Use Agreement signed by your institution.

## Repository structure
```bash
THIS_REPO/
├── LICENSE
├── README.md
└── analyses
    ├── AMP_PD_analyses
    │   ├── 01_AMP_PD_clinical_and_RNA_data_prep.ipynb
    │   ├── 02_AMP_PD_univariate_TWAS_SMR_PTS_scores.ipynb
    │   ├── 03_AMP_PD_univariate_PTS_sensivity_analyses.ipynb
    │   ├── 04_AMP_PD_penalised_regression_multi_PTS_models.ipynb
    │   ├── 05_AMP_PD_XGBoost_multi_PTS_models.ipynb
    │   └── 06_AMP_PD_PRS_analyses.ipynb
    └── sumstats_analyses
        ├── 01_munge
        │   ├── MungeSumstats_for_TRS_project.R
        │   └── submit_MungeSumstats_for_TRS_project.sh
        ├── 02a_TWAS
        │   ├── FUSION_conditional_TWAS_post_process_loop_eQTLGen.sh
        │   └── FUSION_loop_eQTLGen.sh
        ├── 02b_SMR
        │   ├── SMR_loop_eQTLGen.sh
        │   └── clump_SMR_and_SMR_multi_results.R
        ├── 03_GSORA
        │   └── GSORA_TRS.R
        ├── 04_PGS
        │   ├── SBayesRC_calc_hapmap_sumstats_array.sh
        │   └── SBayesRC_impute_hapmap_sumstats_array.sh
        └── 05_MR
            ├── Extract_PD_instruments.R
            ├── Extract_exposure_instruments.R
            ├── MR_clust_PD_LBD.R
            ├── MR_clust_enrichment.R
            ├── MR_get_proxies_run_IVW.R
            ├── MR_get_proxies_run_IVW_reverse.R
            └── MR_senstivity_analyses.R
```

## Analysis Notebooks
### Languages: Python, bash, and R
The repository is divided into two main analysis components:

| **Directory**       | **Notebooks**                                            | **Description**                                                                                  | 
|:--------------------|:---------------------------------------------------------|:-------------------------------------------------------------------------------------------------|
|`AMP_PD_analyses/`   |                                                          | analyses performed using individual-level data from the AMP-PD                                   |
|`AMP_PD_analyses/`   | `01_AMP_PD_clinical_and_RNA_data_prep.ipynb`             | data preparation for running PTS analyses in AMP-PD                                              |
|`AMP_PD_analyses/`   | `02_AMP_PD_univariate_TWAS_SMR_PTS_scores.ipynb`         | calculate and test associations between PTS scores from FUSION, SMR and SMR-multi                |
|`AMP_PD_analyses/`   | `03_AMP_PD_univariate_PTS_sensivity_analyses.ipynb`      | senstivity analyses using conditional/clumped PTS scores                                         |
|`AMP_PD_analyses/`   | `04_AMP_PD_penalised_regression_multi_PTS_models.ipynb`  | LASSO and elastic net models                                                                     |
|`AMP_PD_analyses/`   | `05_AMP_PD_XGBoost_multi_PTS_models.ipynb`               | XGBoost models                                                                                   |
|`AMP_PD_analyses/`   | `06_AMP_PD_PRS_analyses.ipynb`                           | polygenic scoring analyses and machine learning models                                           |
|`sumstats_analyses`  |                                                          | analyses performed using GWAS summary statistics including:                                      |
|`01_munge`           | `MungeSumstats_for_TRS_project.R`                        | summary statistic cleaning                                                                       |
|`01_munge`           | `submit_MungeSumstats_for_TRS_project.sh`                | summary statistic cleaning                                                                       |
|`02a_TWAS`           | `FUSION_conditional_TWAS_post_process_loop_eQTLGen.sh`   | transcriptome-wide association study                                                             |
|`02a_TWAS`           | `FUSION_loop_eQTLGen.sh`                                 | transcriptome-wide association study                                                             |
|`02b_SMR`            | `clump_SMR_and_SMR_multi_results.R`                      | summary-based Mendelian randomization                                                            |
|`02b_SMR`            | `SMR_loop_eQTLGen.sh`                                    | summary-based Mendelian randomization                                                            |
|`03_GSORA`           | `GSORA_TRS.R`                                            | gene set overrepresentation analysis                                                             |
|`04_PGS`             | `SBayesRC_calc_hapmap_sumstats_array.sh`                 | polygenic scoring                                                                                |
|`04_PGS`             | `SBayesRC_impute_hapmap_sumstats_array.sh`               | polygenic scoring                                                                                |
|`05_MR`              | `Extract_exposure_instruments.R`                         | Mendelian randomisation                                                                          |
|`05_MR`              | `Extract_PD_instruments.R`                               | Mendelian randomisation                                                                          |
|`05_MR`              | `MR_clust_enrichment.R`                                  | Mendelian randomisation                                                                          |
|`05_MR`              | `MR_clust_PD_LBD.R`                                      | Mendelian randomisation                                                                          |
|`05_MR`              | `MR_get_proxies_run_IVW_reverse.R`                       | Mendelian randomisation                                                                          |
|`05_MR`              | `MR_get_proxies_run_IVW.R`                               | Mendelian randomisation                                                                          |
|`05_MR`              | `MR_senstivity_analyses.R`                               | Mendelian randomisation                                                                          |

## Computing environments

Summary-statistics analyses were performed on the **Apocrita High Performance Computing (HPC)** cluster at **Queen Mary University of London**.
Individual-level analyses using **AMP-PD** data were performed on the **Terra** cloud platform.

- Apocrita HPC: https://docs.hpc.qmul.ac.uk/
- Terra: https://terra.bio/

## Software
| **Software**                        | **Version** | **Resource URL**                              | **RRID**        | **Notes**                                                    | 
|:------------------------------------|:------------|:----------------------------------------------|:----------------|:-------------------------------------------------------------|
| R Project for Statistical Computing | 4.4.1         | http://www.r-project.org/                     | RRID:SCR_001905 | Used for majority of statistical analyses                    |
| Python Programming Language         | 3.10        | http://www.python.org/                        | RRID:SCR_008394 | Used primarily for running XGBoost                           |
| SMR                                 | 1.3.1       | https://yanglab.westlake.edu.cn/software/smr/ | RRID:SCR_026042 | Used to run summary-based Mendelian randomisation            |
| GCTB                                | 2.5.5       | https://gctbhub.cloud.edu.au/software/gctb/   | NA              | Used for SBayesRC for PGS SNP-weighting                      |
| PLINK                               | 2.0         | https://www.cog-genomics.org/plink/2.0/       | RRID:SCR_001757 | Used for genetic QC and PGS scoring in individual level data |
