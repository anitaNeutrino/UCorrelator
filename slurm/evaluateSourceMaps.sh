#!/bin/bash
#SBATCH --job-name=evalSourceMaps 
#SBATCH --time=06:00:00
#SBATCH --account=pi-avieregg
#SBATCH --partition=broadwl
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4096
#SBATCH --output=./log/%j.log 


source /home/cozzyd/anita/env.sh


echo $@
START=${1-300}
END=${2-310}
MC=${3-0}
DEC=${4-1}


root.exe -b -q macro/evaluateSourceMap.C"($START,$END,$MC,$DEC)" 




