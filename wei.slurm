#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 03:25:00
#SBATCH -A gsienkf
##SBATCH --partition=orion
#SBATCH --partition=bigmem
#SBATCH --job-name=python
#SBATCH --output=log.python
##SBATCH --mem=0

 source ~/intelenv

#ulimit -S unlimited
#ulimit -c unlimited

#python all-obs-jedi.py --output=1 | tee log.all-obs.out

 python tile-dist.py --output=1 | tee log.tile-dist.out

