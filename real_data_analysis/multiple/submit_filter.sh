#!/bin/bash
# Parameters
PHENO_LIST=("height" "bmi" "platelet")

# Slurm parameters
PART=candes,hns,stat,pilanci      # Partition names
MEMO=10G                         # Memory required (50GB)
TIME=00-04:00:00                  # Time required (1d)
CORE=10                           # Cores required (10)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS

# Loop through the configurations
for RESOLUTION in {0..6}; do
  for PHENO in "${PHENO_LIST[@]}"; do

    # Script to be run
    SCRIPT="filter.sh $RESOLUTION $PHENO"

    # Define job name for this chromosome
    JOBN="filter_"$RESOLUTION"_"$PHENO
    OUTF=$LOGS"/"$JOBN".out"
    ERRF=$LOGS"/"$JOBN".err"

    # Assemble slurm order for this job
    ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT

    #Print order
    echo $ORD

    # Submit order
    $ORD

  done
done


 
