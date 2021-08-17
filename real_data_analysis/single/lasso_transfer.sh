#!/bin/bash

POP=$1
RESOLUTION=$2
PHENO=$3

OUT_DIR="/oak/stanford/groups/candes/transfer_gwas/gwas/stats"
mkdir -p $OUT_DIR
ml boost/1.69.0  
ml R/3.5.1
export OPENBLAS_NUM_THREADS=1

OUT_FILE=$OUT_DIR"/stats/lasso_transfer_res"$RESOLUTION"_"$POPULATION"_"$PHENO".txt"
if [[ ! -f $OUT_FILE ]]; then
  Rscript --no-save lasso_transfer.R $POP $RESOLUTION $PHENO
fi

