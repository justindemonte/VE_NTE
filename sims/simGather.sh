#!/bin/bash

#SBATCH -N 1
#SBATCH --mem=2g
#SBATCH -n 1
#SBATCH -t 0:23:00
module load r/4.2.2

sbatch --mem=2g --time=0:23:00 --job-name=simGather R CMD BATCH --no-save --no-restore simGather.R simGather.Rout