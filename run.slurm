#!/bin/bash
#SBATCH -J Job0903
#SBATCH -N 1      # 申请的节点数
#SBATCH --ntasks-per-node=8  # 申请的核数
#SBATCH -p kshcnormal    # 可访问队列

/public/home/jiaoyj/test20220903h200divide/karstprofile  # 可执行程序的路径
