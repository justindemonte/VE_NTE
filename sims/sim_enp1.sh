#!/bin/bash
 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 19g 
#SBATCH --time=0:45:00 
#SBATCH --job-name=est_enp1
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --array=1-3000
module add anaconda/2023.03
module add r/4.4.0
source ../virtualEnv/bin/activate
 
#define variables
numTrials=13
numTimePts=21

simulation=extension
EM=.2
analysis=naive
spline=FALSE
scenario=1

N=50000
NSIM=3000

R CMD BATCH --no-save --no-restore "--args $scenario $analysis $N $numTrials $numTimePts $NSIM $EM $simulation $spline" 203estVE.R 203estVE_enp1.Rout
