#!/bin/bash
 
#SBATCH -N 1
#SBATCH --mem=1g
#SBATCH -n 1
#SBATCH -t 3:00:00
#SBATCH --output=/dev/null
module add anaconda/2023.03
module add r/4.4.0
source /nas/longleaf/home/demjus01/virtualEnv/bin/activate
 
#define variables.  
numTrials=12
JOBID1=$(sbatch --output=/dev/null --parsable --mem=2g --time=0:10:00 --job-name=101process R CMD BATCH --no-save --no-restore "--args $numTrials" 101processItalyData.R 101processItalyData.Rout)

JOBID2=$(sbatch --output=/dev/null --parsable --dependency afterok:$JOBID1 --mem=16g --time=0:10:00 --job-name=102prep    R CMD BATCH --no-save --no-restore 102prepCohortData.R 102prepCohortData.Rout)

JOBID3=$(sbatch --output=/dev/null --parsable --dependency afterok:$JOBID2 --mem=23g --time=0:17:00 --job-name=103models  R CMD BATCH --no-save --no-restore 103tmtCensModels.R 103tmtCensModels.Rout)

JOBID4=$(sbatch --output=/dev/null --parsable --dependency afterok:$JOBID3 --mem=96g --time=1:30:00 --job-name=104prepSTE R CMD BATCH --no-save --no-restore 104prepSTEdata.R 104prepSTEdata.Rout)

#unstandardized analysis
estimator=unstandardized
JOBID5=$(sbatch --output=/dev/null --parsable --dependency afterok:$JOBID4 --mem=64g --time=0:23:00 --job-name=105estVE_un R CMD BATCH --no-save --no-restore "--args $estimator" 105estimateVE.R 105estimateVE_unstd.Rout)

testHypothesis=TRUE
sbatch --output=/dev/null --parsable --dependency afterok:$JOBID5 --mem=64g --time=2:00:00 --array=0-0 --job-name=estVar_hyp_un R CMD BATCH --no-save --no-restore "--args $estimator $testHypothesis" 106estimateVariance.R 106estimateVar_hyp_un.Rout

testHypothesis=FALSE
JOBID6=$(sbatch --output=/dev/null --parsable --dependency afterok:$JOBID5 --mem=64g --time=2:00:00 --array=0-0 --job-name=estVar_trl0_un R CMD BATCH --no-save --no-restore "--args $estimator $testHypothesis" 106estimateVariance.R 106estimateVar_trl0_un.Rout)

#define variables.
numTrialsMinus1=11
sbatch --output=/dev/null --dependency afterok:$JOBID6 --mem=64g --time=2:00:00 --array=1-$numTrialsMinus1 --job-name=estVar_all_un R CMD BATCH --no-save --no-restore "--args $estimator $testHypothesis" 106estimateVariance.R 106estimateVar_all_un.Rout
#end unstandardized analysis

#standardized analysis
estimator=standardized
JOBID7=$(sbatch --output=/dev/null --parsable --dependency afterok:$JOBID4 --mem=96g --time=2:00:00 --job-name=105estVE_st R CMD BATCH --no-save --no-restore "--args $estimator" 105estimateVE.R 105estimateVE_st.Rout)

testHypothesis=TRUE
sbatch --output=/dev/null --parsable --dependency afterok:$JOBID7 --mem=64g --time=5:00:00 --array=0-0 --job-name=estVar_hyp_st R CMD BATCH --no-save --no-restore "--args $estimator $testHypothesis" 106estimateVariance.R 106estimateVar_hyp_st.Rout

testHypothesis=FALSE
JOBID8=$(sbatch --output=/dev/null --parsable --dependency afterok:$JOBID7 --mem=64g --time=3:00:00 --array=0-0 --job-name=estVar_trl0_st R CMD BATCH --no-save --no-restore "--args $estimator $testHypothesis" 106estimateVariance.R 106estimateVar_trl0_st.Rout)

#define variables.
numTrialsMinus1=11
sbatch --output=/dev/null --dependency afterok:$JOBID8 --mem=64g --time=3:00:00 --array=1-$numTrialsMinus1 --job-name=estVar_all_st R CMD BATCH --no-save --no-restore "--args $estimator $testHypothesis" 106estimateVariance.R 106estimateVar_all_st.Rout
#end standardized analysis
