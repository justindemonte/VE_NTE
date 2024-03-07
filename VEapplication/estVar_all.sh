#!/bin/bash
 
#SBATCH -N 1
#SBATCH --mem=1g
#SBATCH -n 1
#SBATCH -t 0:10:00
module add anaconda/2023.03
module add r/4.2.2
source /nas/longleaf/home/demjus01/virtualEnv/bin/activate
 
#define variables.  
numTrials=12
numTrialsMinus1=11

testHypothesis=FALSE
sbatch -p bigmem --qos=bigmem_access --mem=2000g --time=3:00:00 --array=1-$numTrialsMinus1 --job-name=estVar_all R CMD BATCH --no-save --no-restore "--args $testHypothesis" estimateVariance.R estimateVar_all.Rout

sleep 10
