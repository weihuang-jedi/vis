#!/bin/bash
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 03:25:00
#SBATCH -A gsienkf
#SBATCH --partition=bigmem
#SBATCH --job-name=plot
#SBATCH --output=log.plot-incr+obsdiff
##SBATCH --mem=0

 source ~/pythonenv

#ulimit -S unlimited
#ulimit -c unlimited

#python all-obs-jedi.py --output=1
#python tile-dist.py --output=1

 python plot-incr+obsdiff.py --output=1

 trimit

