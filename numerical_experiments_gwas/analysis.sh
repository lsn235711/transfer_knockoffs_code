#!/bin/bash

POP=$1
NCAUSAL=$2
SPEC=$3
SNR=$4
FOLD=$5

OUT_DIR="/oak/stanford/groups/candes/matteo/transfer_knockoffs"
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR"/stats"
mkdir -p $OUT_DIR"/discoveries"
ml boost/1.69.0  
ml R/3.5.1
export OPENBLAS_NUM_THREADS=1

OUT_FILE=$OUT_DIR"/stats/lasso_res5_n"$NCAUSAL"_snr"$SNR"_s"$SPEC"_"$POP"_fold"$FOLD".txt"
if [[ ! -f $OUT_FILE ]]; then
  Rscript --no-save lasso.R $POP $NCAUSAL $SPEC $SNR $FOLD
fi

OUT_FILE=$OUT_DIR"/discoveries/discoveries_res5_n"$NCAUSAL"_snr"$SNR"_s"$SPEC"_"$POP"_fold"$FOLD".txt"
if [[ ! -f $OUT_FILE ]]; then
  Rscript --no-save filter.R $POP $NCAUSAL $SPEC $SNR $FOLD 0
fi
