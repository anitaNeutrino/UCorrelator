#!/bin/bash
#SBATCH --job-name=doAll 
#SBATCH --output=./log/decimated_%j.log 
#SBATCH --time=36:00:00
#SBATCH --account=pi-avieregg
#SBATCH --partition=broadwl
#SBATCH --cpus-per-task=4


echo $@
source /home/cozzyd/anita/env.sh

RUN=$1
N=${2-0}
START={${3-0}
FILTER=${4-sinsub_10_3_ad_2} 

srun bin/doAll $RUN $N $START $FILTER

