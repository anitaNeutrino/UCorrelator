#!/bin/bash
#SBATCH --job-name=doLDB 
#SBATCH --output=./log/%j.log 
#SBATCH --time=03:00:00
#SBATCH --account=pi-avieregg
#SBATCH --partition=broadwl
#SBATCH --cpus-per-task=8


echo $@
source /home/cozzyd/anita/env.sh

RUN=$1
N=${2-0}
FILTER=${3-0}


export OMP_NUM_THREADS=8 

bin/doLDB $RUN $N $FILTER

