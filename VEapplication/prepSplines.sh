#!/bin/bash
 
#SBATCH -N 1
#SBATCH --mem=1g
#SBATCH -n 1
#SBATCH -t 0:05:00
module add anaconda/2023.03
module add r/4.2.2
source /nas/longleaf/home/demjus01/virtualEnv/bin/activate

sbatch --mem=32g --time=0:30:00 --job-name=prepSplines R CMD BATCH --no-save --no-restore prepSplineBasis4deli.R prepSplines.Rout

sleep 10