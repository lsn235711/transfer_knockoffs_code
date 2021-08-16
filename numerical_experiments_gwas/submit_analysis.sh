#!/bin/bash
#
# Class: SLURM dispatcher
#
# Generate knockoffs
# 
# Author: Matteo Sesia
# Date:   06/17/2019

# Parameters
POP_LIST=("African" "African-small" "Asian" "Asian-small" "Indian" "Indian-small" "European" "European-small" "British" "Everyone")
#POP_LIST=("British")
NCAUSAL_LIST=(100)
SPEC_LIST=(0 50 100)
SNR_LIST=(5 10)
FOLD_LIST=$(seq 1 10)

# Slurm parameters
PART=candes,hns,stat,normal  # Partition names
MEMO=20G                     # Memory required (20G)
TIME=00-2:00:00             # Time required (1h)
CORE=5                      # Cores required (5)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

OUT_DIR="/oak/stanford/groups/candes/matteo/transfer_knockoffs"
for FOLD in $FOLD_LIST; do
  for POP in "${POP_LIST[@]}"; do    
    for NCAUSAL in "${NCAUSAL_LIST[@]}"; do    
      for SNR in "${SNR_LIST[@]}"; do    
        for SPEC in "${SPEC_LIST[@]}"; do    
          OUT_FILE=$OUT_DIR"/discoveries/discoveries_res5_n"$NCAUSAL"_snr"$SNR"_s"$SPEC"_"$POP"_fold"$FOLD".txt"
          if [[ ! -f $OUT_FILE ]]; then
            # Script to be run
            SCRIPT="analysis.sh $POP $NCAUSAL $SPEC $SNR $FOLD"
            # Define job name
            JOBN="lasso_"$POP"_"$NCAUSAL"_"$SPEC"_"$SNR"_"$FOLD
            OUTF=$LOGS"/"$JOBN".out"
            ERRF=$LOGS"/"$JOBN".err"
            # Assemble slurm order for this job
            ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
            # Print order
            echo $ORD
            # Submit order
            #$ORD
            #./$SCRIPT
          fi
        done
      done
    done
  done
done
