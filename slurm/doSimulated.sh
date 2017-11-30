#!/bin/bash
#SBATCH --job-name=doSimulated 
#SBATCH --output=./log/doSimulated%j.log 
#SBATCH --time=03:00:00
#SBATCH --account=pi-avieregg
#SBATCH --partition=broadwell
#SBATCH --cpus-per-task=4


echo $@
source /home/cozzyd/anita/env.sh

RUN=$1
N=${2-0}
START=${3-0}
OUTPUT=${4-simulated}
FILTER=${5-sinsub_10_3_ad_2}

export OMP_NUM_THREADS=4 
srun bin/doSimulated $RUN $N $START $OUTPUT $FILTER

