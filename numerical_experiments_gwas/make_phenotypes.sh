#!/bin/bash

POP=$1
NCAUSAL=$2
SPECIFICITY=$3
SNR=$4

mkdir -p "/oak/stanford/groups/candes/matteo/transfer_knockoffs/data"
mkdir -p "/oak/stanford/groups/candes/matteo/transfer_knockoffs/data/phenotypes"
ml boost/1.69.0  
ml R/3.5.1
export OPENBLAS_NUM_THREADS=1
Rscript --no-save make_phenotypes.R $POP $NCAUSAL $SPECIFICITY $SNR
