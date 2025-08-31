#!/bin/bash
 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 23g 
#SBATCH --time=0:45:00 
#SBATCH --job-name=truM
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --array=1-3
module add anaconda/2023.03
module add r/4.4.0
source ../virtualEnv/bin/activate
 
#define variables
numTrials=13
numTimePts=21

simulation=main
EM=0

N=10000000

R CMD BATCH --no-save --no-restore "--args $N $numTrials $numTimePts $EM $simulation" 202getTruth.R 202getTruth_m.Rout
