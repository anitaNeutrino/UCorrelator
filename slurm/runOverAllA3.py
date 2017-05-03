#! /usr/bin/env python

import os 
import sys 


binary = sys.argv[1]

threads_per_task=4
runs_per_task=4


partition = "broadwl" 
if len(sys.argv) > 2: 
    partition = sys.argv[2]



account = "pi-avieregg" 


if len(sys.argv) > 3: 
    account = sys.argv[3]


if len(sys.argv) > 4: 
    treads_per_task = int(sys.argv[4])

if len(sys.argv) > 5: 
    runs_per_task = int(sys.argv[5])


cpus = threads_per_task * runs_per_task 

for start in range(130,439,runs_per_task): 

    cmd="sbatch -p %s -A %s -c %d slurm/parallelHelper.sh %s %d %d %d" % (partition, account, cpus, binary, threads_per_task, start, start+runs_per_task -1) 

    print(cmd) 
    os.system(cmd)
            


