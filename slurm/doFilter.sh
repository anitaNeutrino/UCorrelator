#!/bin/bash
#SBATCH --job-name=evaluteFilters 
#SBATCH --output=./log/%j.log 
#SBATCH --time=24:00:00
#SBATCH --account=pi-avieregg
#SBATCH --partition=broadwl
#SBATCH --cpus-per-task=6


echo $@
source /home/cozzyd/anita/env.sh

RUN=$1
N=${2-0}

bin/evaluateFiltersA3 $RUN $N 

