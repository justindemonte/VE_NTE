library("dplyr")
library("tibble")
library("tidyr")
library("purrr")
library("reticulate")
setwd('..') # move up one directory to level where virtualEnv lives
use_virtualenv('./virtualEnv/', required=TRUE)
setwd(paste0(getwd(), "/sims")) # reset wd to directory where this script lives
###----------------------------------------------------------------
###
###   Derive STE analysis data and analyze one replication of sim dataset
###
###----------------------------------------------------------------
# collect args passed in from SLURM job script
args        <- commandArgs(trailingOnly=TRUE)
scenario    <- as.numeric(args[1])
analysis    <- as.character(args[2])
N           <- as.numeric(args[3])
numTrials   <- as.numeric(args[4])
numTimePts  <- as.numeric(args[5])
NSIM        <- as.numeric(args[6])
ageProdTerm <- as.numeric(args[7])
simulation  <- as.character(args[8])
spl         <- as.logical(args[9])
sim         <- Sys.getenv("SLURM_ARRAY_TASK_ID")
# end collect args passed in from SLURM job script

tau         <- numTimePts - 1
J           <- numTrials  - 1

k_0        <- 5
j_0        <- 5
seTimes_j  <- c(0, 3, 6,  9, 12)
seTimes_k  <- c(1, 4, 8, 12, 15)

source("200helpers.R")
# end boilerplate

# generate simulated cohort
set.seed(sim)
source("201DGP.R")
cohortDataWide <- simd
rm(simd)
###----------------------------------------------------------------
###
###   Create long format data (one row per calendar time) 
###   for observational cohort
###
###----------------------------------------------------------------
for (l in 0:tau){ 
  cohortDataWide <- cohortDataWide %>% 
    mutate( "Y_{l}"     := if_else( Tstar >   l,     0, 1 ),
            "leadZ_{l}" := if_else( V_1   > ( l+1 ), 0, 1 ))  
}
cohortDataLong <- cohortDataWide %>% 
  select( c( "id", "X1_c", "X2", "X3", 
    num_range( "E_",     0:tau ),
    num_range( "Y_",     0:tau ),
    num_range( "leadZ_", 0:tau ))) %>%
  pivot_longer(
    cols        = starts_with( c( "Y_", "leadZ_") ), 
    names_to    = c( ".value", "time" ),
    names_sep   = "\\_" )  %>%
  mutate( l     = as.numeric(time) ) %>% group_by( id ) %>% 
  mutate( lagY  = dplyr::lag(Y,     default=0), 
          Z     = dplyr::lag(leadZ, default=0) )        %>%
  mutate( RsupLeadZ  = if_else( Y==0 & Z==0, 1, 0 )) %>%
          ungroup() %>%
  select( -time ) 
##----------------------------------------------------------------
##  Estimate IP of censoring weights
##----------------------------------------------------------------
# obtain vector of l among those records 'at risk' in the leadZ|X model
l_Zmodel <- (cohortDataLong %>% filter( RsupLeadZ==1 ))$l 
# compute default knot locations for 4 knots 
lKnots   <- Hmisc::rcspline.eval(l_Zmodel, nk=4, knots.only=TRUE)

# append spline basis to dataset
basisMat <- Hmisc::rcspline.eval(cohortDataLong$l, knots=lKnots, nk=4, inclx=TRUE)
colnames(basisMat) <- c("Z_lspl1", "Z_lspl2", "Z_lspl3")
cohortDataLong <- cohortDataLong %>% bind_cols(basisMat) 

if ( spl ){  
  pFit <- glm( leadZ ~ Z_lspl1 + Z_lspl2 + Z_lspl3 + X1_c + X2 + X3,
               family = binomial( link = "logit" ),
               data   = cohortDataLong,
               subset = ( RsupLeadZ==1 ) )  
} else{ # main sims, correctly specified vaccine uptake model (polyl fn of cal time)
  pFit <- glm( leadZ ~ poly( l, 2, raw = TRUE ) + X1_c + X2 + X3,
               family = binomial( link = "logit" ),
               data   = cohortDataLong,
               subset = ( RsupLeadZ==1 ) )  
}
# Compute IPT weight and add corresponding column to analysis dataset
cohortDataLong$predPZ_X <- predict( pFit, newdata = cohortDataLong, type="response")
cohortDataLong          <- cohortDataLong  %>%
   mutate( oneMinusLambdaCinv = case_when(
     Z     == 1 ~ 1,
     leadZ == 1 ~ 1 / predPZ_X,
     leadZ == 0 ~ 1 / ( 1-predPZ_X ),
     TRUE ~ NA_real_
     ))
zParmEst <- coef(pFit)
###----------------------------------------------------------------
###
###   Prepare wide dataset with 1 obs per person for delicatessen
###
###----------------------------------------------------------------
deliData <- cohortDataLong %>%
  select( "id", "Y", "lagY", "leadZ", "Z", "l", "RsupLeadZ",
      "X1_c", "X2", "X3", "Z_lspl1", "Z_lspl2", "Z_lspl3", 
      starts_with("E_") ) %>%
  pivot_wider(
    names_from  = l,
    values_from = c( Y, lagY, leadZ, Z, RsupLeadZ,
                     Z_lspl1, Z_lspl2, Z_lspl3 ),
    names_sep   = "_"
  ) %>% select(-l)
for ( l in 0:tau ){ # generate time variables in wide format
  deliData <- deliData %>% mutate( 'i_{l}' := l, 'iSqr_{l}' := l**2 )
}
###----------------------------------------------------------------
###
###   Create STE analysis dataset
###
###----------------------------------------------------------------
# for each trial j in 0 to J, create a flag indicating that
# a given row belongs in the trial-j dataset (BITJDS) defined as follows:
# if an individual is eligible for trial j,
# then rows with j <= l <= tau belong in the trial-j dataset (BITJDS)
for ( j in 0:J ){
  cohortDataLong[ paste0( "BITJDS_", j) ] <-
    cohortDataLong[ paste0( "E_", j) ] * ( cohortDataLong[ "l" ] >= j )
}
STEdataLong <- cohortDataLong %>%
  select( c( "id", "l", "Y", "lagY", "leadZ", "Z", "oneMinusLambdaCinv", 
             "X1_c", "X2", "X3", 
        num_range( "BITJDS_", 0:J ) )) %>%
  pivot_longer(
    cols      = c( starts_with( "BITJDS_" ) ),
    names_to  = c( ".value", "trialNo" ),
    names_sep = "\\_" )                %>%
  mutate( j   = as.numeric( trialNo )) %>%
  filter( BITJDS == 1 )                %>%
  select( -c( "trialNo", "BITJDS" ) )  %>%
  arrange( id, j, l )                  %>%
  mutate( k = l-j )                    %>% # derive trial time

  group_by( id, j )                    %>%
  mutate( A = leadZ[ k==0 ] )          %>% # treatment assignment for trial j
  mutate( R = case_when (                  # derive at-risk indicator
    ( Z != A ) | ( lagY == 1 ) ~ 0,
    TRUE ~ 1 ))                        %>%
  mutate(prod1minusLambdaCinv = cumprod(oneMinusLambdaCinv)) %>%
  mutate(Wjk = dplyr::lag(prod1minusLambdaCinv) )            %>%
  ungroup()
###----------------------------------------------------------------
###
###   Compute estimated VE
###
###----------------------------------------------------------------
if ( spl ){
  # obtain vector of l among those records at risk in the outcome model
  l_Ymodel <- (STEdataLong %>% filter(( R == 1 ) & ( k != 0 )))$l
  # compute default knot locations for 4 knots 
  lKnots   <- Hmisc::rcspline.eval(l_Ymodel, nk=4, knots.only=TRUE)
  
  # append spline basis to dataset
  basisMat <- Hmisc::rcspline.eval(STEdataLong$l, knots=lKnots, nk=4, inclx=TRUE)
  colnames(basisMat) <- c("lspl1_Y", "lspl2_Y", "lspl3_Y")
  STEdataLong <- STEdataLong %>% bind_cols(basisMat) 
  
  # Check, should be TRUE
  print(all.equal(STEdataLong$l, STEdataLong$lspl1_Y))
  
  # obtain vector of k among those records at risk in the outcome model
  k_Ymodel <- (STEdataLong %>% filter(( R == 1 ) & ( k != 0 )))$k
  # compute default knot locations for 4 knots 
  kKnots   <- Hmisc::rcspline.eval(k_Ymodel, nk=4, knots.only=TRUE)
  
  # append spline basis to dataset
  basisMat <- Hmisc::rcspline.eval(STEdataLong$k, knots=kKnots, nk=4, inclx=TRUE)
  colnames(basisMat) <- c("kspl1_Y", "kspl2_Y", "kspl3_Y")
  STEdataLong <- STEdataLong %>% bind_cols(basisMat) 
  
  # Check, should be TRUE
  print(all.equal(STEdataLong$k, STEdataLong$kspl1_Y))
}
##----------------------------------------------------------------
##  Fit outcome hazard model
##----------------------------------------------------------------
if ( simulation=="main" ){
  if ( analysis=="proposed" ){
    if ( spl ){
      oFit <- glm( Y ~ lspl1_Y + lspl2_Y + lspl3_Y + A +
                     A:lspl1_Y + A:lspl2_Y + A:lspl3_Y +
                     A:kspl1_Y + A:kspl2_Y + A:kspl3_Y,
                   family  = binomial( link = "logit" ),
                   weights = Wjk, 
                   data    = STEdataLong,
                   subset  = ( R == 1 ) & ( k != 0 ))  
    } else{
      oFit <- glm( Y ~ poly( l, 2, raw=TRUE ) + A +
                     A:poly( l, 2, raw=TRUE ) +
                     A:poly( k, 2, raw=TRUE ) ,
                   family  = binomial( link = "logit" ),
                   weights = Wjk, 
                   data    = STEdataLong,
                   subset  = ( R == 1 ) & ( k != 0 ))  
    }
  } else if ( analysis=="naive" ){
      oFit <- glm( Y ~ poly( l, 2, raw=TRUE ) + A +
                   A:poly( k, 2, raw=TRUE ) ,
                 family  = binomial( link = "logit" ),
                 weights = Wjk, 
                 data    = STEdataLong,
                 subset  = ( R == 1 ) & ( k != 0 ))
  } 
} else if ( simulation=="extension" ){

  if ( analysis=="proposed" ){
    oFit <- glm( Y ~ poly( l, 2, raw=TRUE ) + 
                   X1_c + X2 + X3 +
                   A +
                   A:poly( l, 2, raw=TRUE ) +
                   A:poly( k, 2, raw=TRUE ) +
                   A:X1_c,
                 family  = binomial( link = "logit" ),
                 weights = Wjk, 
                 data    = STEdataLong,
                 subset  = ( R == 1 ) & ( k != 0 ))
  } else if ( analysis=="naive" ){
    oFit <- glm( Y ~ poly( l, 2, raw=TRUE ) + A +
                   A:poly( l, 2, raw=TRUE ) +
                   A:poly( k, 2, raw=TRUE ) ,
                 family  = binomial( link = "logit" ),
                 weights = Wjk, 
                 data    = STEdataLong,
                 subset  = ( R == 1 ) & ( k != 0 ))
  }
}
estParams <- coef( oFit )
print(estParams)
##----------------------------------------------------------------
##  Main analysis: Compute hatVE as a function of alphaStar
##----------------------------------------------------------------
if ( simulation == "main" | ( simulation=="extension" & analysis=="naive" ) ){
  if ( spl ){
    
    VEestMat <- matrix(ncol = numTrials, nrow = tau)
    
    STEdataLong$prEvent <- predict(oFit, newdata = STEdataLong)
    STEdataLong_t       <- STEdataLong %>% 
      filter( ( R==1 ) & ( k!=0 ) )    %>% 
      select(A, j, k, prEvent)         %>%
      mutate(prNoEvent = 1-plogis(prEvent))
    
    for ( trl in 0:J ){
      trlIndex      <- trl + 1
      outcomeScore  <- STEdataLong_t %>% filter( j==trl )
      outcomeScore  <- 
        distinct(outcomeScore[c("A", "k", "prNoEvent")], A, k, .keep_all = TRUE) %>% 
        filter( k>0 )                               %>%
        group_by(A)                                 %>%
        mutate( cumPrNoEvent = cumprod(prNoEvent) ) %>%
        ungroup()
      outcomeScore0 <- outcomeScore %>% filter(A==0) %>% rename(surv0 = cumPrNoEvent)
      outcomeScore1 <- outcomeScore %>% filter(A==1) %>% rename(surv1 = cumPrNoEvent)
      
      scoreDsFinal  <- outcomeScore0 %>% full_join(outcomeScore1, by="k") %>%
        mutate( VE = 1 - ( (1-surv1)/(1-surv0) ) )
      
      VEtrl         <- c(scoreDsFinal$VE, rep(NA, trl))
      
      VEestMat[, trlIndex] <- VEtrl
    }
  } else{ # no spline, use helper functions to compute est VE as a function of est alpha parms
    if (simulation=="main"){
      if (analysis=="proposed"){    estHazMatrix <- getHazMatrix(estParams, naive=FALSE)
      } else if (analysis=="naive"){estHazMatrix <- getHazMatrix(estParams, naive=TRUE)}  
    } else{ # simulation="extension"
      estHazMatrix <- getHazMatrix(estParams, naive=FALSE)
    }
    VEestMat       <- getVEmatrix( estHazMatrix )
  }
  ##----------------------------------------------------------------
  ##  Obtain point estimate for beta
  ##----------------------------------------------------------------
  # Calculate AUC_j for j=0, 1, ..., J
  # AUC[j] contains \hat{AUC}_j
  AUC         <- colSums(VEestMat[1:(tau-J), ])
  jVec        <- 0:J
  # regress \hat{AUC}_j on j 
  betaData    <- data.frame(AUC, jVec)
  betaFit     <- lm(AUC ~ jVec, data=betaData)
  print(betaFit)  
}
##----------------------------------------------------------------
##  Extension analysis: Compute hatVE^s as a function of gammaStar
##----------------------------------------------------------------
if ( simulation=="extension" & analysis=="proposed" ){

  risk0mat <- matrix(ncol = numTrials, nrow = tau)
  risk1mat <- matrix(ncol = numTrials, nrow = tau)
  
  riskMatList <- list(risk0mat, risk1mat)
  
  for (a in 0:1){
    # expand cohort data to long format with one row for each (j,k) in \mathcal{W}
    ds_a <- cohortDataWide %>% select(id, X1_c, X2, X3) %>%
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
            ( ds_a %>% subset(j == jind & k == kind) )$risk 
          )
      }
    }
  }
  VEestMat <- 1-( riskMatList[[ 2 ]]/riskMatList[[ 1 ]] )
}
##----------------------------------------------------------------
##  Prepare point estimates for all parameters for deli
##----------------------------------------------------------------
initProp    <- list( zParmEst  ) 
initOutcome <- list( estParams ) 

initRR_j    <- 1 - VEestMat[ k_0, ( seTimes_j+1 ) ]
initRR_k    <- 1 - VEestMat[ seTimes_k, ( j_0+1 ) ]
if ( simulation=="main" ){
  initLogRR <- c( log(initRR_j), log(initRR_k) )  
  initBeta  <- coef(betaFit)[2]
} else{ # extension
  if ( analysis=="proposed" ){
    initPi0   <- riskMatList[[ 1 ]][ k_0, ( seTimes_j+1 ) ]
    initPi1   <- riskMatList[[ 2 ]][ k_0, ( seTimes_j+1 ) ] 
  }
  initLogRR <- log(initRR_j)  
}
##----------------------------------------------------------------
##  Prepare rcs basis matrices
##----------------------------------------------------------------
if ( spl ){
  oSpline_l <- STEdataLong %>% 
    filter(( R==1 ) & ( k!=0 )) %>%
    select(c("lspl1_Y", "lspl2_Y", "lspl3_Y")) %>%
    unique() %>% arrange(lspl1_Y)
  
  oSpline_k <- STEdataLong %>% 
    filter(( R==1 ) & ( k!=0 )) %>%
    select(c("kspl1_Y", "kspl2_Y", "kspl3_Y")) %>%
    unique() %>% arrange(kspl1_Y)
  
  for ( m in 1:tau ){
    deliData <- deliData %>% mutate( # add columns to the deli ds for spline basis matrices
      "oSpline_l_{m}"   := oSpline_l[[m, 1]],
      "oSpline_l'_{m}"  := oSpline_l[[m, 2]],
      "oSpline_l''_{m}" := oSpline_l[[m, 3]],
      "oSpline_k_{m}"   := oSpline_k[[m, 1]],
      "oSpline_k'_{m}"  := oSpline_k[[m, 2]],
      "oSpline_k''_{m}" := oSpline_k[[m, 3]])
  }
}
py_run_file( "203estVar.py" )
##----------------------------------------------------------------
##  Conduct test of TEH assumption
##----------------------------------------------------------------
if ( simulation=="main" & analysis=="proposed" ){ 
  numNuisParm     <- length(zParmEst) + length(estParams)  
  betaIndex       <- numNuisParm + 1
  hatBeta         <- py$deliEst[ betaIndex ]
  asVarBeta       <- py$deliVar[ betaIndex, betaIndex ]
  seBeta          <- sqrt( asVarBeta/N ) 
  waldStat_beta   <- hatBeta/seBeta
  oneSidePval     <- pnorm( waldStat_beta ) 
  pVal            <- c( oneSidePval )
}
##----------------------------------------------------------------
##  Prepare output for one replication
##----------------------------------------------------------------
if ( simulation=="main" ){
  logRRCI  <- tail(py$deliCI, 10)
  VECI     <- 1 - exp(logRRCI)  # VE=1-RR
  CIout    <- VECI[, rev(seq_len(ncol(VECI)))] # reverse cols so confidence lims are lower, upper
  colnames(CIout) <- c("lowerCL","upperCL")
  
  logRRres <- tail(py$deliEst, 10)
  VEres    <- 1 - exp(logRRres) # VE=1-RR
  
  # note these are for variance of log(RR)
  VarRes   <- tail(diag(py$deliVar), 10)
  SEres    <- sqrt(VarRes/N) 
  
} else{ # simulation="extension"
  logRRCI  <- tail( py$deliCI, 5 )
  VECI     <- 1 - exp( logRRCI )  # VE=1-RR
  CIout    <- VECI[, rev(seq_len(ncol(VECI)))] # reverse cols so confidence lims are lower, upper
  colnames(CIout) <- c( "lowerCL","upperCL" )
  
  logRRres <- tail( py$deliEst, 5 )
  VEres    <- 1 - exp( logRRres ) # VE=1-RR
  
  # note these are for variance of log(RR)
  VarRes   <- tail( diag(py$deliVar), 5 )
  SEres    <- sqrt( VarRes/N ) 
}
##----------------------------------------------------------------
##  Write output for single replication of sim
##----------------------------------------------------------------
simCode <- ifelse( simulation=="main", "m", "e" )
anlCode <- ifelse( analysis=="naive",  "n", "p" )
splCode <- ifelse( spl,                "s", "p" )
settingCode <- paste0( simCode, anlCode, splCode, scenario )

if ( simulation=="main" ){
  if ( spl ){
    directory <- paste0( "results/spline/scenario",  scenario, "/", analysis, "/")
  } else{
    directory <- paste0( "results/poly/scenario",    scenario, "/", analysis, "/")
  }
}else{ #simulation=="extension"
  directory <- paste0( "results/extension/scenario", scenario, "/", analysis, "/")
}
saveRDS(scenario,   file=paste0(directory, 'scenario_',   settingCode, '.rds'))
saveRDS(numTrials,  file=paste0(directory, 'numTrials_',  settingCode, '.rds'))
saveRDS(numTimePts, file=paste0(directory, 'numTimePts_', settingCode, '.rds'))
saveRDS(NSIM,       file=paste0(directory, 'NSIM_',       settingCode, '.rds'))
saveRDS(seTimes_j,  file=paste0(directory, 'selectj_',    settingCode, '.rds'))
saveRDS(seTimes_k,  file=paste0(directory, 'selectk_',    settingCode, '.rds'))

SEoutFilename    <- paste0(directory, "SEres_", settingCode, '_', sim, ".csv")
write.csv(SEres,    SEoutFilename, row.names = FALSE)

CIoutFilename    <- paste0(directory, "CIres_", settingCode, '_', sim, ".csv")
write.csv(CIout,    CIoutFilename, row.names = FALSE)

VEoutFilename    <- paste0(directory, "VEres_", settingCode, '_', sim, ".csv")
write.csv(VEres,    VEoutFilename, row.names = FALSE)

if ( simulation=="main" & analysis=="proposed" ){
  hypTestFilename <- paste0(directory, "hypTestRes_", settingCode, '_', sim, '.csv')
  write.csv(pVal, hypTestFilename,  row.names = FALSE)
}