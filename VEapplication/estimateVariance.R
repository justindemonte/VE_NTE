# Justin DeMonte  
# 240305
library("tidyverse")
library("reticulate")
## collect args passed in from SLURM job script
args            <- commandArgs(trailingOnly=TRUE)
testHypothesis  <- as.logical(args[1])
trialNo         <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# end collect args passed in from SLURM job script

# virtual Python environment contains required versions of py libraries
use_virtualenv('/nas/longleaf/home/demjus01/virtualEnv/', required=TRUE)
##----------------------------------------------------------------
## Read in permanent datasets 
##----------------------------------------------------------------
dataDir         <- "data/"
resDir          <- "results/"

numTimePts      <- readRDS(file=paste0(dataDir, "numTimePts.rds"))
numTrials       <- readRDS(file=paste0(dataDir, "numTrials.rds"))

d <- readRDS(file=paste0(dataDir, "d.rds")) %>% select(-starts_with("date")) %>% 
  select(-c("dropoutDate", "covidDeathDate"))

N <- nrow(d)

pFit            <- vector(mode="list", length = 1)
pFit[[1]]       <- readRDS(file=paste0(dataDir, "pFit_pooled.rds"))

cFit            <- vector(mode="list", length = 1)
cFit[[1]]       <- readRDS(file=paste0(dataDir, "cFit_pooled.rds"))

oFit            <- vector(mode="list", length = 1)
oFit[[1]]       <- readRDS(file=paste0(dataDir, 'oFit_calTime.rds'))

tFit            <- c(readRDS(file=paste0(dataDir, 'tFit_pooled.rds')))

logRRmat        <- readRDS(file=paste0(dataDir, 'logRRmat.rds'))
logRR_trialNo   <- logRRmat[, (trialNo + 1)]

logRRinit       <- logRR_trialNo[1:(numTimePts-trialNo)]  

tau             <- numTimePts - 1
J               <- numTrials  - 1

# estimate variance
py_run_file("ItalyVarianceEstimation.py")

if (testHypothesis){
  ##----------------------------------------------------------------
  ##  Conduct test of TEH assumption
  ##----------------------------------------------------------------
  numNuisParm     <- length(pFit[[1]]) + length(cFit[[1]]) + length(oFit[[1]])
  
  betaIndex       <- numNuisParm + 1
  hatBeta         <- py$deliEst[betaIndex]
  asVarBeta       <- py$deliVar[betaIndex, betaIndex]
  seBeta          <- sqrt(asVarBeta/N)
  waldStat_beta   <- hatBeta/seBeta
  oneSidePval     <- pnorm(waldStat_beta) 
  twoSidePval     <- 1 - pnorm(abs(waldStat_beta)) + pnorm(-abs(waldStat_beta))
  
  pVals           <- c(oneSidePval, twoSidePval)
  print(pVals)
  
  hypTestFilename <- paste0(resDir, "hypTestRes.csv")
  write.csv(pVals, hypTestFilename, row.names = FALSE)    
} else{
  ##----------------------------------------------------------------
  ##  Save datasets containing delicatessen output
  ##----------------------------------------------------------------
  saveRDS(N,          file=paste0(dataDir, 'N.rds'))
  saveRDS(py$deliVar, file=paste0(dataDir, 'deliVarAnalysis', trialNo, '.rds'))
  saveRDS(py$deliEst, file=paste0(dataDir, 'deliEstAnalysis', trialNo, '.rds'))
  saveRDS(py$deliCI,  file=paste0(dataDir, 'deliCIanalysis',  trialNo, '.rds'))
}
