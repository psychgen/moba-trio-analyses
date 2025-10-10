#!/bin/bash

# Activate the conda environment for LDSC
#module unload Python/3.10.4-GCCcore-11.3.0
#module unload SciPy-bundle/2022.05-foss-2022a
source /cluster/projects/p471/ldsc-activation

# Set the working directory
cd /tsd/p471/data/durable/projects/qc_paper/trio-GWAS

# Munge data
munge_sumstats.py \
--out scratch/lockebmi \
--merge-alleles eur_w_ld_chr/w_hm3.snplist \
--chunksize 500000 \
--signed-sumstats BETA,0 \
--N-col N \
--a1 A1 \
--a2 A2 \
--frq Freq_Tested_Allele_in_HRS \
--p P \
--snp SNP \
--sumstats scratch/Meta-analysis_Locke_et_al_UKBiobank_2018_UPDATED.txt

munge_sumstats.py \
--out scratch/woodbmi \
--merge-alleles eur_w_ld_chr/w_hm3.snplist \
--chunksize 500000 \
--signed-sumstats BETA,0 \
--N-col N \
--a1 Tested_Allele \
--a2 Other_Allele \
--frq Freq_Tested_Allele_in_HRS \
--p P \
--snp SNP \
--sumstats scratch/Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz

munge_sumstats.py \
--out scratch/EA3 \
--merge-alleles eur_w_ld_chr/w_hm3.snplist \
--chunksize 500000 \
--signed-sumstats BETA,0 \
--N-col N \
--a1 EA \
--a2 OA \
--frq EAF \
--p P \
--snp rsID \
--sumstats scratch/EA3_excl_23andMe_MOBA.meta

munge_sumstats.py \
--out scratch/eversmoke \
--merge-alleles eur_w_ld_chr/w_hm3.snplist \
--chunksize 500000 \
--signed-sumstats Beta,0 \
--N 1251809 \
--a1 A1 \
--a2 A2 \
--frq EAF_A1 \
--p Pval \
--snp MarkerName \
--sumstats scratch/EVER_SMOKER_GWAS_MA_UKB_TAG.txt.nodup.gz

munge_sumstats.py \
--out scratch/danermdd \
--merge-alleles eur_w_ld_chr/w_hm3.snplist \
--chunksize 500000 \
--signed-sumstats OR,1 \
--N-cas-col Nca \
--N-con-col Nco \
--a1 A1 \
--a2 A2 \
--frq FRQ_A_410782 \
--p P \
--sumstats scratch/daner_pgc_mdd_ex.23aME.loo.no.MoBa.gz

# Deactivate the conda environment
source deactivate