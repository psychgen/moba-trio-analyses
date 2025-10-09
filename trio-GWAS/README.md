# Trio GWAS results files (MoBa QC paper)

## LDSC Regression - Genetic correlations

- Note that the 'population' (vanilla) GWAS results are from the model performed on the trios sample only (NOT the full MoBa sample). 
- As such, the sample size is the same for the vanilla and trios GWAS 

## Description of data files used

File path: `qc_paper/trio-GWAS/sumstats`

'Raw' MoBa trio GWAS summary files are saved here: `qc_paper/trio-GWAS/sumstats/raw`

1) File name: model06_trios_phenotype*.basic

    Contents:
    Main association results from non-within/population GWAS on the trio offsprings only (i.e. basic linear analysis not conditional on parent genotype).

    Columns included:
    Chromosome; Predictor; Basepair; A1; A2; Wald_Stat; Wald_P; Effect; SD; A1_Mean; MAF

2) File name: model06_trios_phneotype.trios 

    Contents: 
    Main association results from trios GWAS (i.e., effect of offspring genotype conditional on parents� genotype).

    Columns included:
    Chromosome; Predictor; Basepair; A1; A2; Offspring_Effect; Offspring_SD; Offspring_Z; Offspring_P; Father_Effect; Father_SD; Father_Z; Father_P; Mother_Effect; Mother_SD; Mother_Z; Mother_P

    *phenotype = _height/_exam.

## Prepared summary data files for LDSC are saved here: qc_paper/trio-GWAS/sumstats

1) Prepared summary data files for MoBa trio GWAS results:
   - File names: `exam_pop_ldsc` / `exam_con_ldsc` / `height_pop_ldsc` / `height_con_ldsc`
   - Contents: Reduced summary files including only columns necessary for LDSC (rsID; CHR; BPl A1; A2; Wald_Stat (Z score); Wald_P (Pvalue))

2) Munged MoBa summary data (see also `qc_paper/trio-GWAS/log/` for log files):
   - File names: `exam_pop_munged.sumstats.gz` / `exam_con_munged.sumstats.gz` / `height_pop_munged.sumstats.gz` / `height_con_munged.sumstats.gz`

3) Munged prior GWAS summary data (see also `qc_paper/trio-GWAS/log/` for log files):
   - File names: `EA3.sumstats.gz` / `eversmoke.sumstats.gz` / `lockebmi.sumstats.gz` / `woodbmi.sumstats.gz`

## Description of scripts

File path: `qc_paper/trio-GWAS/scripts`

1) File name: `prep-raw-moba-trios-summary-data.R`
   - Purpose: prepares moba trio gwas summary data files

2) File name: `munge-ldsc-submit.sh`
   - Purpose: munge previous GWAS summary data

3) File name: `ldsc-rg-height-qc.sh`
   - Purpose: munge moba trio GWAS summary data on height and compute genetic correlations

4) File name: `ldsc-rg-exam-qc.sh`
   - Purpose: munge moba trio GWAS summary data on exam performance and compute genetic correlations

## Results files

File path: `qc_paper/trio-GWAS/results`

1) Genetic correlations between MoBa population/vanilla GWAS on height/exam and prior gwas:
   - `exam/heightpop_EA3`
   - `exam/heightpop_eversmoker`
   - `exam/heightpop_lockebmi`
   - `exam/heightpop_woodbmi`

2) Genetic correlations between MoBa trio/conditional GWAS on height/exam and prior gwas:
   - `exam/heightcon_EA3`
   - `exam/heightcon_eversmoker`
   - `exam/heightcon_lockebmi`
   - `exam/heightcon_woodbmi`
