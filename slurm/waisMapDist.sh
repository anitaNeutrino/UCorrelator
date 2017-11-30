#!/bin/bash
#SBATCH --job-name=waisMapDist 
#SBATCH --time=03:00:00
#SBATCH --account=pi-avieregg
#SBATCH --partition=broadwl
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4096
#SBATCH --output=./log/%j.log 


source /home/cozzyd/anita/env.sh


echo $@

root.exe -b -q macro/waisMapDistributions.C"($1,$2)"




