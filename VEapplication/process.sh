#!/bin/bash
 
#SBATCH -N 1
#SBATCH --mem=1g
#SBATCH -n 1
#SBATCH -t 0:05:00
module add anaconda/2023.03
module add r/4.2.2
source /nas/longleaf/home/demjus01/virtualEnv/bin/activate
 
#define variables.  
numTrials=12

sbatch --mem=4g --time=0:23:00 --job-name=process_CD R CMD BATCH --no-save --no-restore "--args $numTrials" processItalyData.R processItalyData.Rout

sleep 10
