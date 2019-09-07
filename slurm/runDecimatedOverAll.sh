#!/bin/sh

#SBATCH --job-name=doAllDecimated 
#SBATCH --time=24:00:00
#SBATCH --output=./log/doAllDecimated_%j.log 
#SBATCH --nodes=24
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --exclusive
#SBATCH --account=pi-avieregg
#SBATCH --partition=broadwl

module load parallel 

START=${1-130}
END=${2-439}
source /home/cozzyd/anita/env.sh
srun="srun --exclusive -N1 -n1 -c$SLURM_CPUS_PER_TASK"
parallel="parallel --delay .2 -j $SLURM_NNODES --joblog runOverAll_$START_$END.log --resume"
$parallel "$srun bin/doDecimated {1} &> log/decimatedRunOver_{1}.log" ::: {$START..$END}

