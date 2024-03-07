#!/bin/bash
 
#SBATCH -N 1
#SBATCH --mem=1g
#SBATCH -n 1
#SBATCH -t 0:01:00
module add anaconda/2023.03
module add r/4.2.2
source /nas/longleaf/home/demjus01/virtualEnv/bin/activate
 
#define variables.  
numTrials=12
numTrialsMinus1=11

sbatch --mem=1g --time=0:05:00 --array=0-$numTrialsMinus1 --job-name=RderiveTrials_CD R CMD BATCH --no-save --no-restore  deriveTrialData.R deriveTrialData.Rout

sleep 10