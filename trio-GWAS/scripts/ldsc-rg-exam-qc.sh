#!/bin/bash
#SBATCH --job-name=trio_gwas
#SBATCH --account=pXXX
#SBATCH --time=2:00:00 
#SBATCH --output=./output/%x_%u_%j.out # this specification is jobname_user_jobID
#SBATCH --error=./error/trio_gwas_submission.txt
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1

# Activate the conda environment for LDSC
source /cluster/projects/ldsc-activation

# Set the working directory
cd /cluster/projects/trio_gwas/moba-qc-rg

# Munge data

munge_sumstats.py \
--out exam_pop_munged \
--merge-alleles eur_w_ld_chr/w_hm3.snplist \
--N 41825 \
--chunksize 500000 \
--sumstats sumstats/exam_pop_filtered_ldsc.txt

munge_sumstats.py \
--out exam_con_munged \
--merge-alleles eur_w_ld_chr/w_hm3.snplist \
--N 41825 \
--chunksize 500000 \
--sumstats sumstats/exam_con_filtered_ldsc.txt

# LD score regression - compute rgs 

# Population trios gwas 
ldsc.py \
--rg exam_pop_munged.sumstats.gz,lockebmi.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out exampop_lockebmi

ldsc.py \
--rg exam_pop_munged.sumstats.gz,EA3.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out exampop_EA3

ldsc.py \
--rg exam_pop_munged.sumstats.gz,woodbmi.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out exampop_woodbmi

ldsc.py \
--rg exam_pop_munged.sumstats.gz,eversmoke.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out exampop_eversmoker

# Conditional trios gwas 
ldsc.py \
--rg exam_con_munged.sumstats.gz,lockebmi.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out examcon_lockebmi

ldsc.py \
--rg exam_con_munged.sumstats.gz,EA3.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out examcon_EA3

ldsc.py \
--rg exam_con_munged.sumstats.gz,woodbmi.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out examcon_woodbmi

ldsc.py \
--rg exam_con_munged.sumstats.gz,eversmoke.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out examcon_eversmoker

# Deactivate the conda environment
source deactivate
