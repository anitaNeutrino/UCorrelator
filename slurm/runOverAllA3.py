#! /usr/bin/env python

import os 
import sys 


binary = sys.argv[1]

partition = "sandyb" 
if len(sys.argv) > 2: 
    partition = sys.argv[2]

account = "pi-avieregg" 

if len(sys.argv) > 3: 
    account = sys.argv[3]


for start in range(130,439,8): 
    print "sbatch -p %s -A %s slurm/parallelHelper.sh %s 2 %d %d" % (partition, account, binary, start, start+8) 
    os.system("sbatch -p %s -A %s slurm/parallelHelper.sh %s 2 %d %d" % (partition, account, binary, start, start+8) )
            


