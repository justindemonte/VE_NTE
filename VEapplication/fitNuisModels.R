# Justin DeMonte  
# 240305

# fit nuisance models
library("tidyverse")
##----------------------------------------------------------------
## Read in permanent datasets 
##----------------------------------------------------------------
dataDir <- "data/"

numTimePts <- readRDS(file=paste0(dataDir, "numTimePts.rds"))
numTrials  <- readRDS(file=paste0(dataDir, "numTrials.rds"))

trialData  <- vector(mode="list", length = numTrials)
for (trl in 1:numTrials){
  trialData[[trl]] <- readRDS(file=paste0(dataDir, "trialData", trl, ".rds"))
}

tau        <- numTimePts - 1
J          <- numTrials  - 1
##----------------------------------------------------------------
##  Set the values for IPT weights
##----------------------------------------------------------------
trialData_pooled <- bind_rows(trialData)

pFit <- glm(Z ~ rms::rcs(trial, 4) + BLage_c + sex + comorbidity,
          family = binomial(link = "logit"),
          data   = trialData_pooled)

pSplineBasis <- unique(model.matrix(pFit)[, # obtain spline basis matrix
  c("rms::rcs(trial, 4)trial",
    "rms::rcs(trial, 4)trial'",
    "rms::rcs(trial, 4)trial''")])

# Compute IPT weight and add corresponding column to analysis dataset
trialData_pooled$predPA_L <- predict(pFit, newdata = trialData_pooled)
trialData_pooled <- trialData_pooled %>%
  mutate(predPA_L = plogis(predPA_L)) %>%
  mutate(IPTweight = if_else(Z==1, 1/predPA_L, 1/(1-predPA_L))) %>%
  select(-c(predPA_L))

# separate data back into trial-specific datasets
for (j in 0:(numTrials-1)){
  jIndex <- j + 1
  trialData[[jIndex]] <- trialData_pooled %>% filter(trial==j)
}

# exploratory analysis
for (j in 0:(numTrials-1)){
  jIndex <- j + 1
  print(cat("trial", j)) 
  print(fivenum(trialData[[jIndex]]$IPTweight))
  
  trialData[[jIndex]] %>% select(IPTweight) %>%
    summarise(mn = mean(IPTweight)) %>%
    print()
}
# End exploratory analysis
##----------------------------------------------------------------
##  Transpose trial data wide to long
##----------------------------------------------------------------
trialData_long <- vector(mode = "list", length = numTrials)
for ( t in 0:(numTrials-1) ) {
  tIndex <- t + 1
  trialData_long[[tIndex]] <- trialData[[tIndex]] %>%
    select(c("id", "trial", "BLage_c", "sex", "comorbidity", "Z", "IPTweight",
            num_range( "E_", 0:(numTimePts-t) ),
            num_range( "C_", 0:(numTimePts-t) ))) %>%
    pivot_longer(
      cols = !c("IPTweight", "trial", "Z", "id", "BLage_c", "sex", "comorbidity"),
      names_to = c(".value", "time"),
      names_sep = "\\_") %>%
    mutate(lagE = lag(E), lagC = lag(C)) %>%
    filter(time == 0 | (lagE == 0 & lagC == 0) ) %>%
    mutate(k = as.numeric(time)) %>%
    mutate(oneMinusC = 1 - C) %>% # create variable for being uncensored
    filter(k != 0 ) %>% # k=0 time point not used in estimation procedure
    mutate(kPlusTrial = k + trial) %>%
    select(-c(lagE, lagC, time))
}
##----------------------------------------------------------------
##  Fit model to estimate probability of censoring
##----------------------------------------------------------------
trialData_long_pooled <- bind_rows(trialData_long) %>%
  # create unique id for each person-trial
  rename(idnum = id) %>% 
  mutate(id    = paste0(idnum, "_t", trial))

cFit <- glm(oneMinusC ~ rms::rcs(kPlusTrial, 4) + Z + BLage_c + sex + comorbidity,
            family = binomial(link = "logit"),
            data   = trialData_long_pooled)  

cSplineBasis <- unique(model.matrix(cFit)[,
    c("rms::rcs(kPlusTrial, 4)kPlusTrial", 
      "rms::rcs(kPlusTrial, 4)kPlusTrial'", 
      "rms::rcs(kPlusTrial, 4)kPlusTrial''")])    

# Construct IPC and combined IP weights
trialData_long_pooled$prNotC <- predict(cFit, newdata = trialData_long_pooled)
trialData_long_pooled <- trialData_long_pooled %>%
  mutate( prNotC = plogis(prNotC)) %>%
  group_by(id) %>% # unique id for each person-trial 
  mutate( cumPrNotC = cumprod(prNotC) ) %>%
  ungroup %>%
  mutate( IPCweight = oneMinusC/cumPrNotC ) %>%
  mutate( IPweight = IPTweight * IPCweight) %>%
  mutate( oneMinusE = 1 - E )

for (trl in 0:(numTrials-1)){ # separate back into trial-specific long datasets
  
  trlIndex <- trl + 1
  trialData_long[[trlIndex]] <- trialData_long_pooled %>% filter(trial==trl)  
}
##----------------------------------------------------------------
##  Fit hazard model and compute plug-in estimates
##----------------------------------------------------------------
trialData_long_pooled <- bind_rows(trialData_long)
  
# initialize matrix of trial-specific VE estimates
VEestMat <- matrix(nrow = numTimePts, ncol = numTrials)

oFit_calTime <- glm(E ~ rms::rcs(kPlusTrial, 4) + Z
    + Z:(rms::rcs(k, 4)) + Z:(rms::rcs(kPlusTrial, 4)),
  family  = binomial(link = "logit"),
  weights = IPweight,
  data    = trialData_long_pooled)

oSplineBasis_kPlusTrial <- unique(model.matrix(oFit_calTime)[,
  c("rms::rcs(kPlusTrial, 4)kPlusTrial", 
    "rms::rcs(kPlusTrial, 4)kPlusTrial'", 
    "rms::rcs(kPlusTrial, 4)kPlusTrial''")])
oSplineBasis_k          <- unique(Hmisc::rcspline.eval(
  trialData_long_pooled$k,     nk=4, inclx = TRUE))

trialData_long_pooled$prEvent   <- predict(oFit_calTime, newdata = trialData_long_pooled)
trialData_long_pooled           <- trialData_long_pooled %>% 
  select(Z, trial, k, prEvent) %>%
  mutate(prNoEvent = 1-plogis(prEvent))

# obtain trial-specific predicted values for all possible combos of Z and k
for (trl in 0:(numTrials-1)){
  trlIndex      <- trl + 1
  outcomeScore  <- trialData_long_pooled %>% filter(trial==trl)
  outcomeScore  <- unique(outcomeScore[c("Z", "k", "prNoEvent")]) %>%  
    group_by(Z) %>%
    mutate( cumPrNoEvent = cumprod(prNoEvent) ) %>%
    ungroup
  outcomeScore0 <- outcomeScore %>% filter(Z==0) %>% rename(surv0 = cumPrNoEvent)
  outcomeScore1 <- outcomeScore %>% filter(Z==1) %>% rename(surv1 = cumPrNoEvent)
  
  scoreDsFinal  <- outcomeScore0 %>% full_join(outcomeScore1, by="k") %>%
    mutate( VE = 1 - ( (1-surv1)/(1-surv0) ) )
  
  VEtrl         <- c(scoreDsFinal$VE, rep(NA, trl))
  
  VEestMat[, trlIndex] <- VEtrl
}
##----------------------------------------------------------------
##  Obtain point estimates for beta
##----------------------------------------------------------------
# Calculate AUC_j for j=0, 1, ..., J
# AUC[j] contains \hat{AUC}_j
AUC         <- colSums(VEestMat[1:(tau-J), ])
jVec        <- 0:J
# regress \hat{AUC}_j on j 
tData       <- data.frame(AUC, jVec)
tFit        <- lm(AUC ~ jVec, data=tData)
##----------------------------------------------------------------
## Save permanent datasets for use downstream
##----------------------------------------------------------------
saveRDS(coef(pFit),          file=paste0(dataDir, 'pFit_pooled.rds'))
saveRDS(pSplineBasis,        file=paste0(dataDir, 'pSplineBasis_pooled.rds'))
saveRDS(coef(cFit),          file=paste0(dataDir, 'cFit_pooled.rds'))
saveRDS(cSplineBasis,        file=paste0(dataDir, 'cSplineBasis_pooled.rds'))
  
saveRDS(coef(oFit_calTime),      file=paste0(dataDir, 'oFit_calTime.rds'))
saveRDS(oSplineBasis_kPlusTrial, file=paste0(dataDir, 'oSplineBasis_kPlusTrial.rds'))
saveRDS(oSplineBasis_k,          file=paste0(dataDir, 'oSplineBasis_k.rds'))  

saveRDS(coef(tFit)[2],       file=paste0(dataDir, 'tFit_pooled.rds'))

initLogRR                    <- log(1-VEestMat)
saveRDS(initLogRR,           file=paste0(dataDir, 'logRRmat.rds'))