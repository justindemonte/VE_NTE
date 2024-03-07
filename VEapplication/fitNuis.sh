#!/bin/bash
 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 0:05:00
#SBATCH --mem=1g
module add anaconda/2023.03
module add r/4.2.2
source /nas/longleaf/home/demjus01/virtualEnv/bin/activate

sbatch --mem=128g --time=1:00:00 --job-name=fitNuis R CMD BATCH --no-save --no-restore  fitNuisModels.R fitNuisModels.Rout

sleep 10