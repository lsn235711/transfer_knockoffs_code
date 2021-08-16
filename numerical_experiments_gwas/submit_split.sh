#!/bin/bash
#
# Class: SLURM dispatcher
#
# Generate knockoffs
# 
# Author: Matteo Sesia
# Date:   06/17/2019

# Parameters
RES_LIST=(5)
POP_LIST=("African" "Asian" "Indian" "European" "British" "Everyone")
CHR_LIST=(1)

# Slurm parameters
PART=candes,hns,stat,normal  # Partition names
MEMO=40G                     # Memory required (40G)
TIME=00-00:10:00             # Time required (10m)
CORE=1                      # Cores required (1)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

OUT_DIR="/oak/stanford/groups/candes/matteo/transfer_knockoffs/data/knockoffs"
for CHR in "${CHR_LIST[@]}"; do    
  for RES in "${RES_LIST[@]}"; do    
    for POP in "${POP_LIST[@]}"; do    
      OUT_FILE=$OUT_DIR"/ukb_gen_chr"$CHR"_res"$RES"_"$POP".bed"
      if [[ ! -f $OUT_FILE ]]; then
        # Script to be run
        SCRIPT="split_data.sh $CHR $RES $POP"
        # Define job name
        JOBN="data_"$CHR"_"$RES"_"$POP
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
