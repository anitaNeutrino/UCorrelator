#!/bin/bash
#SBATCH --job-name=specAvg 
#SBATCH --time=10:00:00
#SBATCH --account=pi-avieregg
#SBATCH --partition=broadwl
#SBATCH --cpus-per-task=2
#SBATCH --output=./log/%j.log 


echo $@
source /home/cozzyd/anita/env.sh

root.exe -b -q  macro/makeAvg.C\($1\) 


