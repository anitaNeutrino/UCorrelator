#!/bin/bash
#SBATCH --job-name=doWais 
#SBATCH --output=./log/%j.log 
#SBATCH --time=03:00:00
#SBATCH --account=kicp
#SBATCH --partition=kicp
#SBATCH --cpus-per-task=4


echo $@
source /home/cozzyd/anita/env.sh

RUN=$1
N=${2-0}
START=${3-0}
FILTER=${4-sinsub_10_3_ad_2}

export OMP_NUM_THREADS=4 
srun bin/doWais $RUN $N $START $FILTER

