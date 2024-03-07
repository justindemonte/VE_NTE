#!/bin/bash
 
#SBATCH -N 1
#SBATCH --mem=1g
#SBATCH -n 1
#SBATCH -t 0:05:00
module add anaconda/2023.03
module add r/4.2.2
source /nas/longleaf/home/demjus01/virtualEnv/bin/activate

sbatch --mem=16g --time=0:23:00 --job-name=deriveDeli R CMD BATCH --no-save --no-restore  deriveDeliData.R deriveDeliData.Rout

sleep 10
