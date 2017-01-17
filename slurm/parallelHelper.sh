#!/bin/sh

#SBATCH --job-name=anitaJob 
#SBATCH --time=36:00:00
#SBATCH --output=./log/anitaJob%j.log #
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive

module load parallel 

PROC=$1
NCPU=$2 
START=$3 
END=$4 

echo $PROC 
echo $NCPU 
echo $(seq $START $END) 

PCT=${2}00\% 
export OMP_NUM_THREADS=$2 
parallel --delay .2 --jobs ${PCT} --joblog log/parallel_helper${1}_${START}_${END}.log "bin/$1 {1} &> log/$1_{1}.log" ::: $(seq $START $END)














