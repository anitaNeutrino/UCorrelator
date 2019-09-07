#!/bin/bash
#SBATCH --job-name=makeSourceMaps 
#SBATCH --time=04:00:00
#SBATCH --account=pi-avieregg
#SBATCH --partition=broadwl
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=4096
#SBATCH --output=./log/makeSourceMaps_%j.log 


source /home/cozzyd/anita/env.sh


echo $@

root.exe -b -q macro/makeSourceMap.C"($1,$2)" 




