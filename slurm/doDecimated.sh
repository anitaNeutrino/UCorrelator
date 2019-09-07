#!/bin/bash
#SBATCH --job-name=doDecimated 
#SBATCH --output=./log/decimated_%j.log 
#SBATCH --time=08:00:00
#SBATCH --account=pi-avieregg
#SBATCH --partition=broadwl
#SBATCH --cpus-per-task=4


echo $@
source /home/cozzyd/anita/env.sh

RUN=$1
N=${2-0}
START={${3-0}
FILTER=${4-adsinsub_2_10_3} 

srun bin/doDecimated $RUN $N $START $FILTER

