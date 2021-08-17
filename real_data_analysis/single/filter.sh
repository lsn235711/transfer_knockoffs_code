#!/bin/bash
RESOLUTION=$1
PHENO=$2

ml R/3.5

Rscript --vanilla filter.R $RESOLUTION $PHENO
