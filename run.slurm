#!/bin/bash
#SBATCH -J Job0903
#SBATCH -N 1      # apply node number
#SBATCH --ntasks-per-node=8  # apply core number
#SBATCH -p kshcnormal    # queue

/public/home/jiaoyj/test20220903h200divide/karstprofile  # the path and the exe
