# MoBa Trio Analyses

This repository contains the analysis code used in the manuscript preprint:

> *The Norwegian Mother, Father, and Child cohort study (MoBa) genotyping data resource: MoBaPsychGen pipeline v.1*  
> Elizabeth C. Corfield, Oleksandr Frei, Alexey A. Shadrin, Zillur Rahman, Aihua Lin, Lavinia Athanasiu, Bayram Cevdet Akdeniz, Laurie Hannigan, Robyn E. Wootton, Chloe Austerberry, Amanda Hughes, Martin Tesli, Lars T. Westlye, Hreinn Stefánsson, Kári Stefánsson, Pål R. Njølstad, Per Magnus, Neil M. Davies, Vivek Appadurai, Gibran Hemani, Eivind Hovig, Tetyana Zayats, Helga Ask, Ted Reichborn-Kjennerud, Ole A. Andreassen, Alexandra Havdahl  
> bioRxiv 2022.06.23.496289; doi: [10.1101/2022.06.23.496289](https://doi.org/10.1101/2022.06.23.496289)

The preprint and published manuscript provide the scientific context and example results.

Scripts in this repository are organised by analysis type:

- `data-prep-universal` — prepare phenotypes and covariates used across analyses
- `trio-GWAS` — trio-GWAS code and scripts to estimate SNP heritabilities (h2) and genetic correlations (rG)
- `pgs-analyses` — trio polygenic score (PGS) and PTDT analyses
- `trio-GCTA` — trio GCTA analyses
- `plotting` — plotting scripts used for figures in the manuscript

The QC pipeline used upstream is available at
[MoBaPsychGen-QC-pipeline](https://github.com/psychgen/MoBaPsychGen-QC-pipeline).

## Contents

- `data-prep-universal/`
  - `scripts/`
    - `prep-phenos-covars_moba-qc.R`
    - `prep-raw-moba-trios-summary-data.R`
- `pgs-analyses/`
  - `scripts/`
    - `01_data_prep.R`
    - `02_ptdt_heightex.R`
    - `03_trioPGS.R`
    - `04_ptdt_eduex.R`
    - `README.md`
    - `readme.Rmd`
- `plotting/`
  - `01_qc_figure.R`
  - `02_pgs_figure.R`
  - `03_gwas_figure.R`
  - `04_gcta_figure.R`
  - `set_aes_pars.R`
- `trio-GCTA/`
  - `README.md`
  - `scripts/`
    - `compute_grm.jl`
    - `data_prep_pheno.R`
    - `fit_trio_cluster_exam_filter.jl`
    - `fit_trio_cluster_height_filter.jl`
    - `post_edit.jl`
    - `trio_output_process.R`
- `trio-GWAS/`
  - `README.md`
  - `scripts/`
    - `filter-summary-data-info-score.R`
    - `ldsc-rg-exam-qc.sh`
    - `ldsc-rg-height-qc.sh`
    - `munge-ldsc-submit.sh`
    - `prep-phenos-covars_moba-qc.R`
    - `prep-raw-moba-trios-summary-data.R`

## License

See `LICENSE` in the repository root for license details.

## Acknowledgements

See the manuscript for contributors, funding, and other acknowledgements.
