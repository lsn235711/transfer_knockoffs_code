#!/bin/bash

POP=$1
NCAUSAL=$2
SPEC=$3
SNR=$4
FOLD=$5

OUT_DIR="/oak/stanford/groups/candes/matteo/transfer_knockoffs"
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR"/discoveries_others"
ml boost/1.69.0  
ml R/3.5.1
export OPENBLAS_NUM_THREADS=1

OUT_FILE=$OUT_DIR"/discoveries_others/discoveries_transfer_2_res5_n"$NCAUSAL"_snr"$SNR"_s"$SPEC"_"$POP"_fold"$FOLD".txt"
if [[ ! -f $OUT_FILE ]]; then
  Rscript --no-save filter_linear.R $POP $NCAUSAL $SPEC $SNR $FOLD
fi

OUT_FILE=$OUT_DIR"/discoveries_others/discoveries_transfer_3_res5_n"$NCAUSAL"_snr"$SNR"_s"$SPEC"_"$POP"_fold"$FOLD".txt"
if [[ ! -f $OUT_FILE ]]; then
  Rscript --no-save filter_adaptive.R $POP $NCAUSAL $SPEC $SNR $FOLD
fi
