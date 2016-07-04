#!/bin/bash
#SBATCH --job-name=doDecimated 
#SBATCH --output=./log/%j.log 
#SBATCH --time=04:00:00
#SBATCH --account=kicp
#SBATCH --partition=kicp
#SBATCH --cpus-per-task=4


echo $@
source /home/cozzyd/anita/env.sh

RUN=$1
N=${2-0}

srun bin/doDecimatedNoFilter $RUN $N

