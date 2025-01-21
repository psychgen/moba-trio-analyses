#!/bin/bash
#
# Job name (this is what the cluster will show as running - keep short):
#SBATCH --job-name=trio_gwas
#Project:
#SBATCH --account=p471
#Wall clock limit (change based on resources required):
#SBATCH --time=2:00:00 
# Output filename customization
# This specification is Jobname_User_JobID
#SBATCH --output=./output/%x_%u_%j.out
#SBATCH --error=./error/trio_gwas_submission.txt
#SBATCH --ntasks=1
# Max memory usage (change based on resources required):
#SBATCH --mem-per-cpu=8G
#Number of CPUs per task. 
#SBATCH --cpus-per-task=1

# Activate the conda environment for LDSC
source /cluster/projects/p471/ldsc-activation

# Set the working directory
cd /tsd/p471/cluster/projects/trio_gwas

# Munge data

#munge_sumstats.py \
#--out height_pop_munged \
#--merge-alleles eur_w_ld_chr/w_hm3.snplist \
#--N 21261 \
#--chunksize 500000 \
#--sumstats moba-qc-rg/height_pop_ldsc.txt

#munge_sumstats.py \
#--out height_con_munged \
#--merge-alleles eur_w_ld_chr/w_hm3.snplist \
#--N 21261 \
#--chunksize 500000 \
#--sumstats moba-qc-rg/height_con_ldsc.txt

munge_sumstats.py \
--out lockebmi \
--merge-alleles eur_w_ld_chr/w_hm3.snplist \
--chunksize 500000 \
--N-col N \
--a1 Tested_Allele \
--a2 Other_Allele \
--sumstats moba-qc-rg/Meta-analysis_Locke_et_al_UKBiobank_2018_UPDATED.txt.nodup.gz

munge_sumstats.py \
--out woodbmi \
--merge-alleles eur_w_ld_chr/w_hm3.snplist \
--chunksize 500000 \
--N-col N \
--a1 Tested_Allele \
--a2 Other_Allele \
--sumstats moba-qc-rg/Meta-analysis_Wood_et_al+UKBiobank_2018.txt

munge_sumstats.py \
--out EA3 \
--merge-alleles eur_w_ld_chr/w_hm3.snplist \
--chunksize 500000 \
--a1 EA \
--a2 OA \
--ignore Z_unadj,BETA \
--N-col N \
--sumstats moba-qc-rg/EA3_excl_23andMe_MOBA.meta

munge_sumstats.py \
--out eversmoke \
--merge-alleles eur_w_ld_chr/w_hm3.snplist \
--chunksize 500000 \
--N 1251809 \
--sumstats moba-qc-rg/EVER_SMOKER_GWAS_MA_UKB_TAG.txt.nodup.gz

# Deactivate the conda environment
source deactivate
