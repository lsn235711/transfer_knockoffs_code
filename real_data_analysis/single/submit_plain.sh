#!/bin/bash

# Parameters
POP_LIST=("African" "Asian" "Indian")
#POP_LIST=("Asian")
#PHENO_LIST=("height" "bmi" "platelet" "diabetes" "sbp")
PHENO_LIST=("cvd" "sbp")


# Slurm parameters
PART=candes,hns,stat,normal  # Partition names
MEMO=20G                     # Memory required (50G)
TIME=00-2:00:00             # Time required (1d)
CORE=10                      # Cores required (10)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

# Loop through the configurations
OUT_DIR="/oak/stanford/groups/candes/transfer_gwas/gwas/stats"
for RESOLUTION in {0..6}; do
  for POP in "${POP_LIST[@]}"; do
    for PHENO in "${PHENO_LIST[@]}"; do
      OUT_FILE=$OUT_DIR"/lasso_transfer_res"$RESOLUTION"_"$POPULATION"_"$PHENO  
      if [[ ! -f $OUT_FILE ]]; then
          # Script to be run
          SCRIPT="lasso_plain.sh $POP $RESOLUTION $PHENO"
            
          # Define job name
          JOBN="lasso_plain_"$POP"_"$RESOLUTION"_"$PHENO
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
