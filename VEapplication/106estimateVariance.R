# 106estimateVariance.R
# 
# args: estimator ("standardized" or "unstandardized")
# testHypothesis  (TRUE or FALSE)
#
# input: 
#
# 1. numTrials, a constant, number of trials to be emulated (i.e., J-1)
# /nas/longleaf/home/demjus01/VEapplication/data/numTrials.rds
# 
# 2. numTimePoints - a constant 
# /nas/longleaf/home/demjus01/VEapplication/data/numTimePts.rds
#
# 3. VUPM model object
# /nas/longleaf/home/demjus01/VEapplication/data/dFit.rds 
#
# 4. LTFU model object 
# /nas/longleaf/home/demjus01/VEapplication/data/hFit.rds 
#
# 5. estimator-specific coef ests for outcome model 
# /nas/longleaf/home/demjus01/VEapplication/data/[estimator]/initY_[estimator].rds 
#
# 6. estimator-specific logRR estimates 
# /nas/longleaf/home/demjus01/VEapplication/data/[estimator]/initLogRR_[estimator].rds
#
# 7. deliData 
# /nas/longleaf/home/demjus01/VEapplication/data/deliData.rds 
#
# 8. risk estimates for a \in {0,1} (standardized estimator only)
# /nas/longleaf/home/demjus01/VEapplication/data/standardized/riskaest.rds
#
# Output:  
#
# 1. N a constant, number of individuals in the analytic cohort
# /nas/longleaf/home/demjus01/VEapplication/data/[estimator]/N_[est].rds 
#
# 2. estimator-specific variance estimates
# /nas/longleaf/home/demjus01/VEapplication/data/[estimator]/deliVarAnalysis_[est].rds 
#
# 3. estimator-specific point estimates
# /nas/longleaf/home/demjus01/VEapplication/data/[estimator]/deliEstAnalysis_[est].rds 
#
# 4. estimator-specific CIs
# /nas/longleaf/home/demjus01/VEapplication/data/[estimator]/deliCIAnalysis_[est].rds 
library("tidyverse")
library("reticulate")
setwd('..') # move up one directory to level where virtualEnv lives
use_virtualenv('./virtualEnv/', required=TRUE)
setwd(paste0(getwd(), "/VEapplication")) # reset wd to directory where this script lives

# collect args passed in from SLURM job script
args           <- commandArgs(trailingOnly=TRUE)
estimator      <- as.character(args[1])
testHypothesis <- as.logical(args[2])
trialNo        <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# end collect args passed in from SLURM job script

dataDir        <- "data/"
# read in user-provided constants from upstream program
numTimePts     <- readRDS(file = paste0(dataDir, "numTimePts.rds"))
numTrials      <- readRDS(file = paste0(dataDir, "numTrials.rds"))

tau            <- numTimePts - 1
J              <- numTrials  - 1
K_J            <- tau-J

# read in model fit objects from upstream program
dFit        <- readRDS(file = paste0(dataDir, 'dFit.rds'))
hFit        <- readRDS(file = paste0(dataDir, 'hFit.rds'))

# read in matrix of trial-specific VE estimates
if (estimator == "unstandardized"){
  inDir     <- "data/unstandardized/"
  # read in unstandardized outcome model fit object
  oFit      <- readRDS(file = paste0(inDir, 'oFit_unstandardized.rds'))
  # read in permanent copy of unstandardized VE point estimates
  VEestMat  <- readRDS(file = paste0(inDir, 'VEest_unstandardized.rds'))

} else if (estimator == "standardized"){
  inDir     <- "data/standardized/"
  # read in standardized outcome model fit object
  oFit      <- readRDS(file = paste0(inDir, 'oFit_standardized.rds'))
  # read in permanent copy of standardized VE point estimates and risk estimates
  VEestMat  <- readRDS(file = paste0(inDir, 'VEest_standardized.rds'))
  risk0Mat  <- readRDS(file = paste0(inDir, 'risk0est.rds'))
  risk1Mat  <- readRDS(file = paste0(inDir, 'risk1est.rds'))
}
# read in wide ds formatted for use in delicatessen
deliData    <- readRDS(file = paste0(dataDir, 'deliData.rds'))

N           <- nrow(deliData)
# end boilerplate
##----------------------------------------------------------------
##  Obtain point estimates for beta
##----------------------------------------------------------------
# Calculate AUC_j for j=0, 1, ..., J
# AUC[j] contains \hat{AUC}_j
AUC           <- colSums(VEestMat[1:(tau-J), ])
jVec          <- 0:J
# regress \hat{AUC}_j on j 
betaData      <- data.frame(AUC, jVec)
betaFit       <- lm(AUC ~ jVec, data=betaData)
##----------------------------------------------------------------
##  Prepare point estimates for deli
##----------------------------------------------------------------
initD         <- list( coef(dFit) ) 
initH         <- list( coef(hFit) ) 
initO         <- list( coef(oFit) ) # estimates of alpha*
if ( estimator == "standardized" ){
  if ( testHypothesis ){
    initPi0     <- as.vector(risk0Mat[ 1:K_J, ])
    initPi1     <- as.vector(risk1Mat[ 1:K_J, ])
  } else {
    initPi0     <- risk0Mat[ 1:( tau-trialNo ), (trialNo + 1) ]
    initPi1     <- risk1Mat[ 1:( tau-trialNo ), (trialNo + 1) ]    
  }
}
initRR_j      <- 1 - VEestMat[ , (trialNo + 1) ]
initLogRR     <- log( initRR_j )[ 1:( tau-trialNo ) ]
initBeta      <- coef(betaFit)[2]

py_run_file( "106estimateVariance.py" )
##----------------------------------------------------------------
##  Save datasets containing delicatessen output
##----------------------------------------------------------------
if (estimator == "standardized"){

  outDir <- "data/standardized/"
  resDir <- "resTable/standardized/"

} else if (estimator == "unstandardized"){

  outDir <- "data/unstandardized/"
  resDir <- "resTable/unstandardized/"
}
if (testHypothesis){
  ##----------------------------------------------------------------
  ##  Conduct test of TEH assumption
  ##----------------------------------------------------------------
  hatBeta         <- tail(py$deliEst, 1)
  asVarBeta       <- tail(diag(py$deliVar), 1)
  seBeta          <- sqrt(asVarBeta/N)
  waldStat_beta   <- hatBeta/seBeta
  oneSidePval     <- pnorm(waldStat_beta)

  pVal            <- c(oneSidePval)
  print(pVal)

  if (estimator == "unstandardized"){
    hypTestFilename <- paste0(resDir, "hypTestRes_un.csv")
  } else if (estimator == "standardized"){
    hypTestFilename <- paste0(resDir, "hypTestRes_st.csv")
  }
  write.csv(pVal, hypTestFilename, row.names = FALSE)
} else{
  ##----------------------------------------------------------------------
  ##  Output point estimates and confidence intervals for trial trialNo
  ##----------------------------------------------------------------------
  if (estimator == "standardized"){

    saveRDS(py$N,       file=paste0(outDir, 'N_std.rds'))
    saveRDS(py$deliVar, file=paste0(outDir, 'deliVarAnalysis_std', trialNo, '.rds'))
    saveRDS(py$deliEst, file=paste0(outDir, 'deliEstAnalysis_std', trialNo, '.rds'))
    saveRDS(py$deliCI,  file=paste0(outDir, 'deliCIanalysis_std',  trialNo, '.rds'))

  } else if (estimator == "unstandardized"){

    saveRDS(py$N,       file=paste0(outDir, 'N_unstd.rds'))
    saveRDS(py$deliVar, file=paste0(outDir, 'deliVarAnalysis_unstd', trialNo, '.rds'))
    saveRDS(py$deliEst, file=paste0(outDir, 'deliEstAnalysis_unstd', trialNo, '.rds'))
    saveRDS(py$deliCI,  file=paste0(outDir, 'deliCIanalysis_unstd',  trialNo, '.rds'))
  }
}