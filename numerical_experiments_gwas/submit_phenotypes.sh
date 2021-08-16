#!/bin/bash
#
# Class: SLURM dispatcher
#
# Generate knockoffs
# 
# Author: Matteo Sesia
# Date:   06/17/2019

# Parameters
POP_LIST=("African" "Asian" "Indian" "European" "British")
NCAUSAL_LIST=(100)
SPEC_LIST=(0 50 100)
SNR_LIST=(5 10)

# Slurm parameters
PART=candes,hns,stat,normal  # Partition names
MEMO=10G                     # Memory required (10G)
TIME=00-00:10:00             # Time required (10m)
CORE=1                      # Cores required (1)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

OUT_DIR="/oak/stanford/groups/candes/matteo/transfer_knockoffs/data/phenotypes"
for POP in "${POP_LIST[@]}"; do    
  for NCAUSAL in "${NCAUSAL_LIST[@]}"; do    
    for SNR in "${SNR_LIST[@]}"; do    
      for SPEC in "${SPEC_LIST[@]}"; do    
        OUT_FILE=$OUT_DIR"/phenotypes_n"$NCAUSAL"_snr"$SNR"_s"$SPEC"_"$POP".txt"
        if [[ ! -f $OUT_FILE ]]; then
          # Script to be run
          SCRIPT="make_phenotypes.sh $POP $NCAUSAL $SPEC $SNR"
          # Define job name
          JOBN="pheno_"$POP"_"$SNR"_"$NCAUSAL"_"$SPEC
          OUTF=$LOGS"/"$JOBN".out"
          ERRF=$LOGS"/"$JOBN".err"
          # Assemble slurm order for this job
          ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
          # Print order
          echo $ORD
          # Submit order
          $ORD
          #./$SCRIPT
        fi
      done
    done
  done
done
