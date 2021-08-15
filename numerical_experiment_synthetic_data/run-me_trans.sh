#!/bin/bash

for i in {1..300}
do
   sbatch --export=seed=$i script_trans.sh
done


  
