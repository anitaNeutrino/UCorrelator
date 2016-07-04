#!/bin/bash
#SBATCH --job-name=doDecimated 
#SBATCH --output=./log/%j.log 
#SBATCH --time=04:00:00
#SBATCH --account=kicp
#SBATCH --partition=kicp
#SBATCH --cpus-per-task=2


echo $@
source /home/cozzyd/anita/env.sh

RUN=$1
N=${2-0}
SINSUB=${3-true}

srun bin/doDecimated $RUN $N SINSUB

