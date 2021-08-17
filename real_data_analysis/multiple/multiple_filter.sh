#!/bin/bash
RESOLUTION=$1

ml R/3.5

Rscript --vanilla multiple_filter.R $RESOLUTION
