# 105estimateVE.R
# 
# args: estimator ("standardized" or "unstandardized")
#
# input: 
# 1. STE analysis dataset
# /nas/longleaf/home/demjus01/VEapplication/data/STEdataLong.rds 
#
# 2. numTrials, a constant, number of trials to be emulated (i.e., J-1)
# /nas/longleaf/home/demjus01/VEapplication/data/numTrials.rds
# 
# 3. numTimePoints - a constant 
# /nas/longleaf/home/demjus01/VEapplication/data/numTimePts.rds
#
# 4. numTrials - a constant 
# /nas/longleaf/home/demjus01/VEapplication/data/numTrials.rds
#
# 5. analytic cohort data (as described in Section 2) in wide format 
# /nas/longleaf/home/demjus01/VEapplication/data/elderdata.rds 
#
# 6. spline bases 
# /nas/longleaf/home/demjus01/VEapplication/data/l_Ybasis.rds.rds 
# /nas/longleaf/home/demjus01/VEapplication/data/k_Ybasis.rds.rds 
# /nas/longleaf/home/demjus01/VEapplication/data/X1_Ybasis.rds.rds 
#
# output:  
# 1. estimator-specific outcome model fit object
# /nas/longleaf/home/demjus01/VEapplication/data/[estimator]/oFit_[estimator].rds
# 
# 2. estimator-specific matrix of VE point estimates
# /nas/longleaf/home/demjus01/VEapplication/data/[estimator]/VEest_[estimator].rds
library("tidyverse")

## collect args passed in from SLURM job script
args           <- commandArgs(trailingOnly=TRUE)
estimator      <- as.character(args[1])

dataDir        <- "data/"
# read in user-provided constants from upstream program
numTimePts     <- readRDS(file = paste0(dataDir, "numTimePts.rds"))
numTrials      <- readRDS(file = paste0(dataDir, "numTrials.rds"))

tau            <- numTimePts - 1
J              <- numTrials  - 1

mathcalW <- 0 # |mathcal{W}| number of (j,k) time points in NTE study period
for (v in 0:J){
  mathcalW <- mathcalW + (tau-v)
}

STEdataLong    <- readRDS(file = paste0(dataDir, "STEdataLong.rds"))
cohortDataWide <- readRDS(file = paste0(dataDir, "elderData.rds"))

l_Ybasis       <- readRDS(file = paste0(dataDir, 'l_Ybasis.rds'))
k_Ybasis       <- readRDS(file = paste0(dataDir, 'k_Ybasis.rds'))
X1_Ybasis      <- readRDS(file = paste0(dataDir, 'X1_Ybasis.rds'))

# initialize matrix of trial-specific VE estimates
VEestMat       <- matrix(nrow = tau, ncol = numTrials)
##----------------------------------------------------------------
##  Fit 'o'utcome hazard model 
##----------------------------------------------------------------
if (estimator=="standardized"){
  
  oFit <- glm( Y ~ lspl1_Y + lspl2_Y   + lspl3_Y + 
                 X1_cspl1  + X1_cspl2  + X1_cspl3 +
                 X2 + X3 + X1X3 + A +
                 A:lspl1_Y + A:lspl2_Y + A:lspl3_Y +
                 A:kspl1_Y + A:kspl2_Y + A:kspl3_Y,
              family  = binomial( link = "logit" ),
              weights = Wjk, 
              data    = STEdataLong,
              subset  = ( ( R == 1 ) & ( k != 0 ) )
  )
} else if (estimator=="unstandardized"){
  
  oFit <- glm(  Y ~ lspl1_Y + lspl2_Y + lspl3_Y + A +
                  A:lspl1_Y + A:lspl2_Y + A:lspl3_Y +
                  A:kspl1_Y + A:kspl2_Y + A:kspl3_Y,
               family  = binomial( link = "logit" ),
               weights = Wjk,   
               data    = STEdataLong,
               subset  = ( ( R == 1 ) & ( k != 0 ) ))
}
print(oFit)
##-------------------------------------------------------------------------------
##  obtain trial-specific predicted values for all possible combos of A and k
##-------------------------------------------------------------------------------
if ( estimator=="unstandardized" ){
  STEdataLong$prEvent <- predict(oFit, newdata = STEdataLong)
  STEdataLong_t       <- STEdataLong %>% 
    filter( ( R==1 ) & ( k!=0 ) )    %>% 
    select(A, j, k, prEvent)         %>%
    mutate(prNoEvent = 1-plogis(prEvent)) %>% 
    select(-prEvent)
  
  for ( trl in 0:J ){
    trlIndex      <- trl + 1
    outcomeScore  <- STEdataLong_t %>% filter( j == trl )
    outcomeScore  <- unique(outcomeScore[c("A", "k", "prNoEvent")]) %>% 
      filter(k != 0) %>%
      group_by(A) %>%
      mutate( cumPrNoEvent = cumprod(prNoEvent) ) %>%
      ungroup()
    outcomeScore0 <- outcomeScore %>% filter(A==0) %>% rename(surv0 = cumPrNoEvent)
    outcomeScore1 <- outcomeScore %>% filter(A==1) %>% rename(surv1 = cumPrNoEvent)
    
    scoreDsFinal  <- outcomeScore0 %>% full_join(outcomeScore1, by="k") %>%
      mutate( VE = 1 - ( (1-surv1)/(1-surv0) ) )
    
    VEtrl         <- c(scoreDsFinal$VE, rep(NA, trl)) # trial j has tau-j time points
    
    VEestMat[, trlIndex] <- VEtrl
  }
} else if ( estimator=="standardized" ){
  
  cohortDataWide  <- cohortDataWide %>% mutate( X1X3 = X1_c * X3) 
  
  risk0mat        <- matrix(ncol = numTrials, nrow = tau)
  risk1mat        <- matrix(ncol = numTrials, nrow = tau)
  
  riskMatList     <- list(risk0mat, risk1mat)
  
  for (a in 0:1){
    # expand cohort data to long format with one row for each (j,k) in \mathcal{W}
    ds_a <- cohortDataWide %>% mutate( X1X3 = X1_c * X3) %>% 
      select(id, X1_c, X2, X3, X1X3) %>%
      slice(rep(1:n(), each = (J+1))) %>% # repeat each id J+1 times
      group_by(id)                    %>%
      mutate(j=row_number()-1)        %>% # add values for j \in {0,...,J}
      ungroup()                       %>%
      mutate(K_j = tau-j)             %>%
      tidyr::uncount(K_j)             %>% # expand each id, j combo K_j times 
      group_by(id, j)                 %>%
      mutate(k=row_number())          %>% # add values for k \in {1,...,K_j}
      ungroup()                       %>%
      mutate(l = j+k)
    
    ds_a$A <- a  # set treatment to a for all individuals
    
    ds_a <- ds_a %>% 
      left_join(X1_Ybasis, by=join_by(X1_c==X1_cspl1), relationship="many-to-one", keep=TRUE) %>%
      left_join(k_Ybasis,  by=join_by(k==kspl1_Y),     relationship="many-to-one", keep=TRUE) %>%
      left_join(l_Ybasis,  by=join_by(l==lspl1_Y),     relationship="many-to-one", keep=TRUE) 
    
    # check
    # each individual should have |\mathcal{W}|=tau+(tau-1)+...+(tau-J) entries. 
    # so nrow(ds_a) should be |\mathcal{W}| times nrow(cohortDataWide)
    # should be TRUE
    nrow(ds_a) == nrow(cohortDataWide)*mathcalW
    
    # check: ds_a cols X1_cspl1, X1_cspl2, X1_cspl3 should have no NAs
    sum(is.na(ds_a$X1_cspl1))
    sum(is.na(ds_a$X1_cspl2))
    sum(is.na(ds_a$X1_cspl3))
    
    # check: ds_a cols kspl1_Y, kspl2_Y, kspl3_Y should have no NAs
    sum(is.na(ds_a$kspl1_Y))
    sum(is.na(ds_a$kspl2_Y))
    sum(is.na(ds_a$kspl3_Y))
    
    # check: ds_a cols lspl1_Y, lspl2_Y, lspl3_Y should have no NAs
    sum(is.na(ds_a$lspl1_Y))
    sum(is.na(ds_a$lspl2_Y))
    sum(is.na(ds_a$lspl3_Y))
    
    ds_a$prEvent <- predict(oFit, newdata = ds_a, type="response")
    
    ds_a <- ds_a %>% mutate( oneMinusLambdaAhat = (1-prEvent) ) %>% 
      group_by( id, j ) %>%
      mutate( oneMinusRisk = cumprod(oneMinusLambdaAhat) ) %>%
      ungroup() %>%
      mutate( risk = 1-oneMinusRisk )
    
    for ( jind in 0:J ){
      for ( kind in 1:(tau-jind) ){
        riskMatList[[ (a+1) ]][kind, ( jind+1 )] <- 
          mean(
            ( ds_a %>% subset( j == jind & k == kind ) )$risk 
          )
      }
    }
  }
  VEestMat <- 1-( riskMatList[[ 2 ]]/riskMatList[[ 1 ]] )
}
print(VEestMat)
##----------------------------------------------------------------
##  Save permanent copies of estimator-specific output
##----------------------------------------------------------------
if (estimator == "unstandardized"){
  outDir <- "data/unstandardized/"
  # save unstandardized outcome model fit object
  saveRDS(oFit,     file = paste0(outDir, 'oFit_unstandardized.rds'))
  # save permanent copy of VE point estimates
  saveRDS(VEestMat, file = paste0(outDir, 'VEest_unstandardized.rds'))
  
} else if (estimator == "standardized"){
  outDir <- "data/standardized/"
  # save standardized outcome model fit object
  saveRDS(oFit,     file = paste0(outDir, 'oFit_standardized.rds'))
  # save permanent copy of VE point estimates
  saveRDS(VEestMat, file = paste0(outDir, 'VEest_standardized.rds'))
  # save permanent copy of risk point estimates
  saveRDS(riskMatList[[ 1 ]], file = paste0(outDir, 'risk0est.rds'))
  saveRDS(riskMatList[[ 2 ]], file = paste0(outDir, 'risk1est.rds'))
}