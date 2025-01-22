#!/bin/bash
#SBATCH --job-name=trio_gwas
#SBATCH --account=p471
#SBATCH --time=2:00:00 
#SBATCH --output=./output/%x_%u_%j.out
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
--out height_pop_munged \
--merge-alleles eur_w_ld_chr/w_hm3.snplist \
--N 23810 \
--chunksize 500000 \
--sumstats sumstats/height_pop_filtered_ldsc.txt

munge_sumstats.py \
--out height_con_munged \
--merge-alleles eur_w_ld_chr/w_hm3.snplist \
--N 23810 \
--chunksize 500000 \
--sumstats sumstats/height_con_filtered_ldsc.txt

# LD score regression - compute rgs 

# Population trios gwas 
ldsc.py \
--rg height_pop_munged.sumstats.gz,lockebmi.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out heightpop_lockebmi

ldsc.py \
--rg height_pop_munged.sumstats.gz,EA3.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out heightpop_EA3

ldsc.py \
--rg height_pop_munged.sumstats.gz,woodbmi.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out heightpop_woodbmi

ldsc.py \
--rg height_pop_munged.sumstats.gz,eversmoke.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out heightpop_eversmoker

# Conditional trios gwas 
ldsc.py \
--rg height_con_munged.sumstats.gz,lockebmi.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out heightcon_lockebmi

ldsc.py \
--rg height_con_munged.sumstats.gz,EA3.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out heightcon_EA3

ldsc.py \
--rg height_con_munged.sumstats.gz,woodbmi.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out heightcon_woodbmi

ldsc.py \
--rg height_con_munged.sumstats.gz,eversmoke.sumstats.gz \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out heightcon_eversmoker

# Deactivate the conda environment
source deactivate
