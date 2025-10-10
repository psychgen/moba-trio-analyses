#!/bin/bash

# Usage: ./ldsc_rg_all_mqc.sh <pop_sumstats.tbl> <trios_sumstats.tbl> <out_prefix>
#set -e

if [ $# -ne 3 ]; then
  echo "Usage: $0 <pop_sumstats.tbl> <trios_sumstats.tbl> <output_prefix>"
  exit 1
fi

# Activate the conda environment for LDSC
#module unload Python/3.10.4-GCCcore-11.3.0
#module unload SciPy-bundle/2022.05-foss-2022a
source /cluster/projects/p471/ldsc-activation

# Set the working directory
cd /tsd/p471/data/durable/projects/qc_paper/trio-GWAS

POP_SUMSTATS=$1
TRIOS_SUMSTATS=$2
OUT_PREFIX=$3

MERGE_ALLELES="eur_w_ld_chr/w_hm3.snplist"
REF_LD_CHR="eur_w_ld_chr/"
W_LD_CHR="eur_w_ld_chr/"


OUT1="${POP_SUMSTATS%.*}"
OUT2="${TRIOS_SUMSTATS%.*}"

# Munge both summary stats
munge_sumstats.py \
--out "${OUT1}" \
--merge-alleles "${MERGE_ALLELES}" \
--chunksize 500000 \
--signed-sumstats Effect,0 \
--N-col Neff \
--a1 A1 \
--a2 A2 \
--p P \
--ignore SNP \
--sumstats "${POP_SUMSTATS}"

munge_sumstats.py \
--out "${OUT2}" \
--merge-alleles "${MERGE_ALLELES}" \
--chunksize 500000 \
--signed-sumstats Effect,0 \
--N-col Neff \
--a1 A1 \
--a2 A2 \
--p P \
--ignore SNP \
--sumstats "${TRIOS_SUMSTATS}"

# Run univariate LDSC
ldsc.py \
  --h2 "${OUT1}.sumstats.gz" \
  --ref-ld-chr "${REF_LD_CHR}" \
  --w-ld-chr   "${W_LD_CHR}" \
  --out "results/ldsc_aug25/test/${OUT_PREFIX}_pop_ldsc.ldsc"

ldsc.py \
  --h2 "${OUT2}.sumstats.gz" \
  --ref-ld-chr "${REF_LD_CHR}" \
  --w-ld-chr   "${W_LD_CHR}" \
  --out "results/ldsc_aug25/test/${OUT_PREFIX}_trio_ldsc.ldsc"

# Run rg
ldsc.py \
  --rg "${OUT1}.sumstats.gz,${OUT2}.sumstats.gz" \
  --ref-ld-chr "${REF_LD_CHR}" \
  --w-ld-chr   "${W_LD_CHR}" \
  --out "results/ldsc_aug25/test/${OUT_PREFIX}_rg_pop_v_trio"

# Population gwas 
ldsc.py \
--rg "${OUT1}.sumstats.gz,scratch/lockebmi.sumstats.gz" \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out "results/ldsc_aug25/test/${OUT_PREFIX}_pop_lockebmi"

ldsc.py \
--rg "${OUT1}.sumstats.gz,scratch/woodheight.sumstats.gz" \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out "results/ldsc_aug25/test/${OUT_PREFIX}_pop_woodheight"

ldsc.py \
--rg "${OUT1}.sumstats.gz,scratch/EA3.sumstats.gz" \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out "results/ldsc_aug25/test/${OUT_PREFIX}_pop_EA3"

ldsc.py \
--rg "${OUT1}.sumstats.gz,scratch/eversmoke.sumstats.gz" \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out "results/ldsc_aug25/test/${OUT_PREFIX}_pop_eversmoker"

ldsc.py \
--rg "${OUT1}.sumstats.gz,scratch/danermdd.sumstats.gz" \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out "results/ldsc_aug25/test/${OUT_PREFIX}_pop_mdd"

# Trios gwas 
ldsc.py \
--rg "${OUT2}.sumstats.gz,scratch/lockebmi.sumstats.gz" \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out "results/ldsc_aug25/test/${OUT_PREFIX}_trio_lockebmi"

ldsc.py \
--rg "${OUT2}.sumstats.gz,scratch/woodheight.sumstats.gz" \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out "results/ldsc_aug25/test/${OUT_PREFIX}_trio_woodheight"

ldsc.py \
--rg "${OUT2}.sumstats.gz,scratch/EA3.sumstats.gz" \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out "results/ldsc_aug25/test/${OUT_PREFIX}_trio_EA3"

ldsc.py \
--rg "${OUT2}.sumstats.gz,scratch/eversmoke.sumstats.gz" \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out "results/ldsc_aug25/test/${OUT_PREFIX}_trio_eversmoker"

ldsc.py \
--rg "${OUT2}.sumstats.gz,scratch/danermdd.sumstats.gz" \
--ref-ld-chr eur_w_ld_chr/ \
--w-ld-chr eur_w_ld_chr/ \
--out "results/ldsc_aug25/test/${OUT_PREFIX}_trio_mdd"

cat "Finished"

# Deactivate the conda environment
source deactivate