#!/bin/bash

CHR=$1
RES=$2
POP=$3

INP_DIR="/oak/stanford/groups/candes/popstruct/analysis/knockoffs"
OUT_DIR="/oak/stanford/groups/candes/matteo/transfer_knockoffs/data"
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR"/knockoffs"
OUT_FILE=$OUT_DIR"/knockoffs/ukb_gen_chr"$CHR"_res"$RES"_"$POP

RUN_SPLIT_KNOCKOFFS=1
RUN_COMPUTE_MAF=1

if [[ $RUN_SPLIT_KNOCKOFFS -eq 1 ]]; then
  GENO_FILE=$INP_DIR"/ukb_gen_chr"$CHR"_ibd1_res"$RES
  plink \
    --bed $GENO_FILE".bed" \
    --bim $GENO_FILE".bim" \
    --fam $GENO_FILE".fam" \
    --keep $OUT_DIR"/qc/samples_"$POP".txt" \
    --keep-allele-order \
    --make-bed \
    --memory 39000 \
    --out $OUT_FILE
  rm -f $OUT_FILE".log"
  rm -f $OUT_FILE".nosex"
fi

if [[ $RUN_COMPUTE_MAF -eq 1 ]]; then
  plink \
    --bed $OUT_FILE".bed" \
    --bim $OUT_FILE".bim" \
    --fam $OUT_FILE".fam" \
    --freq \
    --memory 1000 \
    --out $OUT_FILE
  rm -f $OUT_FILE".log"
  rm -f $OUT_FILE".nosex"
fi

