#!/bin/bash
 
#SBATCH -N 1
#SBATCH --mem=1g
#SBATCH -n 1
#SBATCH -t 0:10:00
module add anaconda/2023.03
module add r/4.2.2
source /nas/longleaf/home/demjus01/virtualEnv/bin/activate
 
#define variables
N=50000
numTrials=13
numTimePts=21
NSIM=3000

# polynomial models
scenario=1

analysis=naive
sbatch --output=/dev/null --error=/dev/null --mem=30g --time=0:45:00 --array=1-$NSIM --job-name=sc1_naive R CMD BATCH --no-save --no-restore "--args $scenario $analysis $N $numTrials $numTimePts $NSIM" simEstVE_poly.R simSc1naive_poly.Rout

analysis=proposed
sbatch --output=/dev/null --error=/dev/null --mem=30g --time=0:45:00 --array=1-$NSIM --job-name=sc1_prop  R CMD BATCH --no-save --no-restore "--args $scenario $analysis $N $numTrials $numTimePts $NSIM" simEstVE_poly.R simSc1prop_poly.Rout

scenario=2

analysis=naive
sbatch --output=/dev/null --error=/dev/null --mem=30g --time=0:45:00 --array=1-$NSIM --job-name=sc2_naive R CMD BATCH --no-save --no-restore "--args $scenario $analysis $N $numTrials $numTimePts $NSIM" simEstVE_poly.R simSc2naive_poly.Rout

analysis=proposed
sbatch --output=/dev/null --error=/dev/null --mem=30g --time=0:45:00 --array=1-$NSIM --job-name=sc2_prop  R CMD BATCH --no-save --no-restore "--args $scenario $analysis $N $numTrials $numTimePts $NSIM" simEstVE_poly.R simSc2prop_poly.Rout

scenario=3

analysis=naive
sbatch --output=/dev/null --error=/dev/null --mem=30g --time=0:45:00 --array=1-$NSIM --job-name=sc3_naive R CMD BATCH --no-save --no-restore "--args $scenario $analysis $N $numTrials $numTimePts $NSIM" simEstVE_poly.R simSc3naive_poly.Rout

analysis=proposed
sbatch --output=/dev/null --error=/dev/null --mem=30g --time=0:45:00 --array=1-$NSIM --job-name=sc3_prop  R CMD BATCH --no-save --no-restore "--args $scenario $analysis $N $numTrials $numTimePts $NSIM" simEstVE_poly.R simSc3prop_poly.Rout

# spline models
scenario=1
sbatch --output=/dev/null --error=/dev/null --mem=30g --time=1:15:00 --array=1-$NSIM --job-name=estimateVE_sc1_spline R CMD BATCH --no-save --no-restore "--args $scenario $N $numTrials $numTimePts $NSIM" simEstVE_spline.R simEstVE_spline_Sc1.Rout

scenario=2
sbatch --output=/dev/null --error=/dev/null --mem=30g --time=1:15:00 --array=1-$NSIM --job-name=estimateVE_sc2_spline R CMD BATCH --no-save --no-restore "--args $scenario $N $numTrials $numTimePts $NSIM" simEstVE_spline.R simEstVE_spline_Sc2.Rout

scenario=3
sbatch --output=/dev/null --error=/dev/null --mem=30g --time=1:15:00 --array=1-$NSIM --job-name=estimateVE_sc3_spline R CMD BATCH --no-save --no-restore "--args $scenario $N $numTrials $numTimePts $NSIM" simEstVE_spline.R simEstVE_spline_Sc3.Rout
