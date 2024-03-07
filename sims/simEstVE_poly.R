# Justin DeMonte
# 240305
library("dplyr")
library("tibble")
library("tidyr")
library("purrr")
library("reticulate")
use_virtualenv('/nas/longleaf/home/demjus01/virtualEnv/', required=TRUE)
###########################################################################
### 
###   Generate one simulation 
### 
###########################################################################
## collect args passed in from SLURM job script
args       <- commandArgs(trailingOnly=TRUE)
scenario   <- as.numeric(args[1])
analysis   <- as.character(args[2])
N          <- as.numeric(args[3])
numTrials  <- as.numeric(args[4])
numTimePts <- as.numeric(args[5])
NSIM       <- as.numeric(args[6])
sim        <- Sys.getenv("SLURM_ARRAY_TASK_ID")
# # end collect args passed in from SLURM job script
tau        <- numTimePts - 1
J          <- numTrials  - 1

k_0        <- 5
j_0        <- 5
seTimes_j  <- c(0, 3,  6,  9, 12)
seTimes_k  <- c(1, 4,  8, 12, 15)

if (scenario==1){truParams <- c(-4,    0,     0, -2.5,   0,    0,  .02, .005)}
if (scenario==2){truParams <- c(-4, -.01, -.003, -2.5, .02,  .006,   0,    0)}
if (scenario==3){truParams <- c(-4, -.01, -.003, -2.5, .02,  .006, .02, .005)}

set.seed(sim)
source('simGenerate.R')

trialData <- bind_rows(trialData) # pool trial data
##----------------------------------------------------------------
##  Set the values for IPT weights
##----------------------------------------------------------------
pFit <- glm(Z ~ poly(trial, 2, raw = TRUE) + L1_c + L2 + L3,
          family = binomial(link = "logit"),
          data = trialData)

# Compute IPT weight and add corresponding column to analysis dataset
trialData$predPA_L  <- predict(pFit, newdata = trialData)
trialData           <- trialData %>%
   mutate(predPA_L  =  plogis(predPA_L)) %>%
   mutate(IPTweight =  if_else(Z==1, 1/predPA_L, 1/(1-predPA_L))) %>%
   select(-c(predPA_L))

# after obtaining predicted propensity scores from the pooled data
# separate data back into trial-specific datasets
tempList       <- vector(mode = "list", length = numTrials)
for (j in 0:(numTrials-1)){
  jIndex       <- j + 1
  tempList[[jIndex]] <- trialData %>% filter(trial==j)
}
rm(trialData)
trialData      <- tempList
rm(tempList)

trialData_long <- vector(mode = "list", length = numTrials)
for (trial in 0:J) {
  trialIndex   <- trial + 1
  ##----------------------------------------------------------------
  ##  Transpose trial data wide to long
  ##----------------------------------------------------------------
  trialData_long[[trialIndex]] <- trialData[[trialIndex]] %>%
     select(c("id", "trial", "L1_c", "L2", "L3", "Z", "IPTweight",
              num_range( "E_", 0:(numTimePts-trial-1) ),
              num_range( "C_", 0:(numTimePts-trial-1) ))) %>%
     pivot_longer(
                  cols = !c("IPTweight", "trial", "Z", "id", "L1_c", "L2", "L3"),
                  names_to = c(".value", "time"),
                  names_sep = "\\_"
                  ) %>%
     mutate(lagE = lag(E), lagC = lag(C))         %>%
     filter(time == 0 | (lagE == 0 & lagC == 0) ) %>%
     mutate(k    = as.numeric(time))              %>%
     mutate(oneMinusC = 1-C) %>% # create variable for being uncensored
     filter(k    != 0 )                           %>%
     mutate(kPlusTrial = k + trial)               %>%
     select(-c(lagE, lagC, time))
}

poolData_temp_long <- bind_rows(trialData_long)
rm(trialData_long)
trialData_long     <- poolData_temp_long %>%
  rename(idnum = id) %>% # create new id for each unique person-trial
  mutate(id = paste(idnum, "_t", trial, sep = ""))
rm(poolData_temp_long)

##----------------------------------------------------------------
##  Fit model to estimate probability of censoring
##----------------------------------------------------------------
tempData           <- vector(mode="list", length = 2)
tempData[[1]]      <- trialData_long %>% filter(Z==0)
tempData[[2]]      <- trialData_long %>% filter(Z==1)

cFit <- glm(oneMinusC ~ poly(trial, 2, raw = TRUE) + L1_c + L2 + L3,
           family  = binomial(link = "logit"),
           data    = tempData[[1]])

tempData[[1]]$prNotC <- predict(cFit, newdata = tempData[[1]])
tempData[[1]]       <- tempData[[1]] %>%
  mutate( prNotC = plogis(prNotC)) %>%
  group_by(id) %>% # unique id for each person-trial when pooling
  mutate( cumPrNotC = cumprod(prNotC) ) %>%
  ungroup %>%
  mutate( IPCweight = oneMinusC/cumPrNotC )

trialData_long      <- bind_rows( tempData[[1]], 
    tempData[[2]]   %>% mutate(IPCweight=1) ) %>%
  mutate(IPweight   = IPTweight * IPCweight) %>%
  mutate(oneMinusE  = 1 - E )

##----------------------------------------------------------------
##  Fit outcome model and compute plug-in estimates
##----------------------------------------------------------------
VEestMat            <- matrix(ncol=numTrials, nrow=tau)
colnames(VEestMat)  <- vector(length = numTrials)

if (analysis=="naive"){
  
  oFit <- glm(E ~ poly(kPlusTrial, 2, raw=TRUE) + Z + Z:poly(k, 2, raw=TRUE), 
              family  = binomial(link = "logit"),
              weights = IPweight,
              data    = trialData_long)
  
} else if (analysis=="proposed"){
  
  oFit <- glm(E ~ poly(kPlusTrial, 2, raw=TRUE) + Z +
                Z:poly(kPlusTrial, 2, raw=TRUE) +
                Z:poly(k, 2, raw=TRUE),
              family    = binomial(link = "logit"),
              weights   = IPweight,
              data      = trialData_long)
} 
##----------------------------------------------------------------
##  Calculate estimated VE
##----------------------------------------------------------------
estParams   <- coef(oFit)
if (analysis=="proposed"){    estHazMatrix <- getHazMatrix(estParams, naive=FALSE)
} else if (analysis=="naive"){estHazMatrix <- getHazMatrix(estParams, naive=TRUE)}
VEestMat    <- getVEmatrix(estHazMatrix)
##----------------------------------------------------------------
##  Obtain point estimate for beta
##----------------------------------------------------------------
# Calculate AUC_j for j=0, 1, ..., J
# AUC[j] contains \hat{AUC}_j
AUC         <- colSums(VEestMat[1:(tau-J), ])
jVec        <- 0:J
# regress j on \hat{AUC}_j
tData       <- data.frame(AUC, jVec)
tFit        <- lm(AUC ~ jVec, data=tData)
##----------------------------------------------------------------
##  Prepare point estimates for all nuisance parameters for deli
##----------------------------------------------------------------
initProp    <- list(coef(pFit))    # estimates of zeta
initCens    <- list(coef(cFit))    # estimates of kappa
initOutcome <- list(coef(oFit))    # estimates of alpha*
initBeta    <- c(coef(tFit)[2])    # estimate  of beta

initRR_j    <- 1 - VEestMat[k_0, (seTimes_j + 1)]
initRR_k    <- 1 - VEestMat[seTimes_k, (j_0 + 1)]
initLogRR   <- c(log(initRR_j), log(initRR_k))

py_run_file("simEstVar_poly.py")
##----------------------------------------------------------------
##  Conduct test of TEH assumption
##----------------------------------------------------------------
numNuisParm <- length(coef(pFit)) + length(coef(cFit)) + length(coef(oFit))

if (analysis=="proposed"){

  betaIndex       <- numNuisParm + 1
  hatBeta         <- py$deliEst[betaIndex]
  asVarBeta       <- py$deliVar[betaIndex, betaIndex]
  seBeta          <- sqrt(asVarBeta/N)
  waldStat_beta   <- hatBeta/seBeta
  oneSidePval     <- pnorm(waldStat_beta) 
  twoSidePval     <- 1 - pnorm(abs(waldStat_beta)) + pnorm(-abs(waldStat_beta))
  
  pVals           <- c(oneSidePval, twoSidePval)
}
##----------------------------------------------------------------
##  Prepare output for one replication
##----------------------------------------------------------------
logRRCI  <- tail(py$deliCI, 10)
VECI     <- 1 - exp(logRRCI)  # VE=1-RR
CIout    <- VECI[, rev(seq_len(ncol(VECI)))] # reverse cols so confidence lims are lower, upper
colnames(CIout) <- c("lowerCL","upperCL")

logRRres <- tail(py$deliEst, 10)
VEres    <- 1 - exp(logRRres) # VE=1-RR
trueVE   <- c(VEtrueMat[k_0, (seTimes_j+1)], VEtrueMat[seTimes_k, (j_0+1)])

# note these are for variance of log(RR)
VarRes   <- tail(diag(py$deliVar), 10)
SEres    <- sqrt(VarRes/N)
##----------------------------------------------------------------
##  Write output to csv files for single replication of sim
##----------------------------------------------------------------
directory <- paste0("results/poly/scenario", scenario, "/", analysis, "/")

saveRDS(scenario,   file=paste0(directory, 'scenario_sc',   scenario, '_', analysis, '.rds'))
saveRDS(analysis,   file=paste0(directory, 'analysis_sc',   scenario, '_', analysis, '.rds'))
saveRDS(N,          file=paste0(directory, 'N_sc',          scenario, '_', analysis, '.rds'))
saveRDS(numTrials,  file=paste0(directory, 'numTrials_sc',  scenario, '_', analysis, '.rds'))
saveRDS(numTimePts, file=paste0(directory, 'numTimePts_sc', scenario, '_', analysis, '.rds'))
saveRDS(NSIM,       file=paste0(directory, 'NSIM_sc',       scenario, '_', analysis, '.rds'))

saveRDS(seTimes_k,  file=paste0(directory, 'selectk_sc',    scenario, '_', analysis, '.rds'))
saveRDS(seTimes_j,  file=paste0(directory, 'selectj_sc',    scenario, '_', analysis, '.rds'))

if (analysis=="proposed"){
  hypTestFilename     <- paste0(directory, "hypTestRes_", sim, '_sc', scenario, '_', analysis, ".csv")
  write.csv(pVals, hypTestFilename, row.names = FALSE)
}
SEoutputFilename  <- paste0(directory, "SEresults_", sim, '_sc',  scenario, '_', analysis, ".csv")
write.csv(SEres,     SEoutputFilename, row.names = FALSE)
CIoutputFilename  <- paste0(directory, "CIresults_", sim, '_sc',  scenario, '_', analysis, ".csv")
write.csv(CIout,     CIoutputFilename, row.names = FALSE)
VEoutputFilename  <- paste0(directory, "VEresults_", sim, '_sc',  scenario, '_', analysis, ".csv")
write.csv(VEres,     VEoutputFilename, row.names = FALSE)
logRRoutFilename  <- paste0(directory, "logRRres_", sim, '_sc',  scenario, '_', analysis, ".csv")
write.csv(logRRres,  logRRoutFilename, row.names = FALSE)
VEtrueFilename    <- paste0(directory, "trueVE_sc", scenario, '_', analysis, ".csv")
write.csv(trueVE,    VEtrueFilename,   row.names = FALSE)