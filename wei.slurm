#!/bin/bash
#SBATCH --ntasks-per-node=40
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 01:25:00
#SBATCH -A gsienkf
##SBATCH --partition=orion
#SBATCH --partition=bigmem
#SBATCH --job-name=getkf
#SBATCH --output=log.getkf
##SBATCH --mem=0

 source ~/intelenv

 ulimit -S unlimited
 ulimit -c unlimited

 python tile-dist.py --output=1

