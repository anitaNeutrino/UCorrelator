#!/bin/bash
#SBATCH --job-name=fakeAnitaEvent 
#SBATCH --time=03:00:00
#SBATCH --account=pi-avieregg
#SBATCH --partition=broadwl
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2048
#SBATCH --output=./log/fake_%j.log 


source /home/cozzyd/anita/env.sh


echo $@
RUN=${1-342}
AMP=${2-0}
root.exe -b -q macro/makeFakeEventPointingTree.C"($RUN,$AMP)"




