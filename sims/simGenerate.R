# Justin DeMonte 
# 240305
source('00_balancing_numerical_solver_source_code.R')

# Read in parameter values for distributions of L, Z_j, E_j.
probL2givenL1   <- read.csv("probL2givenL1.csv")
probL3givenL1L2 <- read.csv("probL3givenL1L2.csv")
propensityFit   <- read.csv("propensityFit.csv")
outcomeFit      <- read.csv("outcomeFit.csv")
betaVec         <- c(outcomeFit$x[outcomeFit$X=="age_c"],
                     outcomeFit$x[outcomeFit$X=="sex"],
                     outcomeFit$x[outcomeFit$X=="comorbidity"])

###########################################################################
### 
###  True params and VE
### 
###########################################################################
##----------------------------------------------------------------
##  Helper functions
##----------------------------------------------------------------
# discrete hazard at a single time point as a function of regression params
getHaz       <- function(Z, vaccTime, k, betaVec) {
  plogis( betaVec[1] + ( betaVec[2] * (vaccTime+k) ) + ( betaVec[3] * ((vaccTime+k)^2) ) +
         (betaVec[4] * Z) + (betaVec[5] * (vaccTime+k) * Z) + (betaVec[6] * ((vaccTime+k)^2) * Z) +
         (betaVec[7] * k * Z) + (betaVec[8] * (k^2) * Z) )
}

# takes parameter vector of length 6
# Use this version to compute estimated haz from naive model parameter estimates
getHaz_naive <- function(Z, vaccTime, k, betaVec) { 
     plogis( betaVec[1] +     ( betaVec[2] * (vaccTime+k) ) +    ( betaVec[3] * ((vaccTime+k)^2) ) + 
            (betaVec[4] * Z) + (betaVec[5] * k * Z)  + (betaVec[6] * (k^2) * Z))
}

# obtain trial-specific VE over k.  Takes a vector, returns a vector,
getVE        <- function(haz1, haz0){
  1 - ( (1-cumprod(1-haz1)) / (1-cumprod(1-haz0)) )
}

# Obtain all hazards for Z \in {0, 1} 
# col j+1 of hazMatz holds haz_j^z(k) for k=1 to tau-j
getHazMatrix <- function(params, naive){
  hazMat1    <- matrix(nrow=tau, ncol=tau)
  hazMat0    <- matrix(nrow=tau, ncol=tau)
  for (j in 0:(tau-1)){
    for (k in 1:(tau-j)){
      if (naive==FALSE){
        hazMat1[k, (j+1)] <- getHaz(Z=1,       vaccTime=j, k=k,  betaVec=params)
        hazMat0[k, (j+1)] <- getHaz(Z=0,       vaccTime=j, k=k,  betaVec=params)   
      } else{ # for computing estimated haz from naive model parameter estimates
        hazMat1[k, (j+1)] <- getHaz_naive(Z=1, vaccTime=j, k=k,  betaVec=params)
        hazMat0[k, (j+1)] <- getHaz_naive(Z=0, vaccTime=j, k=k,  betaVec=params)  
      }   
    }
  }
  return(list(hazMat1, hazMat0))
}

# Obtain matrix of true VE across all j and k 
# hazMatList must have hazard matrix for z=1 as first element, 
# hazard matrix for z=0 as second element.
getVEmatrix   <- function(hazMatList){
  VEtrueMat   <- matrix(nrow=tau, ncol=numTrials)
  for (j in 0:J){
    VEtrueMat[1:(tau-j), (j+1)] <- 
      getVE(hazMatList[[1]][1:(tau-j), j+1], hazMatList[[2]][1:(tau-j), j+1])
  }
  return(VEtrueMat)
} 
##----------------------------------------------------------------
##  Calculate true values
##----------------------------------------------------------------
trueHazMatrix <- getHazMatrix(truParams, naive = FALSE)
VEtrueMat     <- getVEmatrix(trueHazMatrix)

# generate vector of true counterfactual hazards across calendar time 
# separately for each possible vacc uptake time (incl. never vacced)
# cfHazard[[1]]   is for those who never uptake treatment
# cfHazard[[j+2]] is for those who uptake treatment at calendar time j
cfHazard      <- vector(mode = "list", length = numTimePts+1)  
cfHazard[[1]] <- trueHazMatrix[[2]][1:tau, 1] # never vacced
cfHazard[[2]] <- trueHazMatrix[[1]][1:tau, 1] # vacc at j=0
for (j in 1:(tau-1)){ 
  cfHazard[[j+2]] <- c(trueHazMatrix[[2]][1:j, 1], trueHazMatrix[[1]][1:(tau-j), j+1])
}
# vacced at last time point => same haz as never vacced:
cfHazard[[numTimePts+1]] <- trueHazMatrix[[2]][1:tau, 1] 
###########################################################################
###
###  Build the simulated dataset
###
###########################################################################
### Covariates
simd <- tibble(id = c(1:N)) %>%
  add_column(intercept = rep(1,N)) %>%
  add_column(L1 = VGAM::rfoldnorm(n = nrow(.), mean = 0, sd = 7, a1 = 1, a2 = 1) + 80) %>%
  # mutate(meanL1 = mean(L1)) %>%
  # ungroup() %>%
  mutate(L1_c = L1 - 86.2) %>% # centered on population mean
  mutate(prL2 = plogis(probL2givenL1$x[probL2givenL1$X=="(Intercept)"] +
      (probL2givenL1$x[probL2givenL1$X=="age_c"] * L1_c))) %>%
  add_column(runifL2 = runif(nrow(.))) %>%
  mutate(L2 = if_else( (runifL2 < prL2), 1, 0 )) %>%
  mutate(prL3 = plogis(
    probL3givenL1L2$x[probL3givenL1L2$X=="(Intercept)"] +
      (probL3givenL1L2$x[probL3givenL1L2$X=="age_c"] * L1_c) +
      (probL3givenL1L2$x[probL3givenL1L2$X=="sex"] * L2))) %>%
  add_column(runifL3 = runif(nrow(.))) %>%
  mutate(L3 = if_else( (runifL3 < prL3), 1, 0))

### Vaccine uptake
# 1st element in list is unvacced, denoted by vaccTime=NA
# 2nd element in vaccTime 0, 3rd vaccTime 1 and so on
vaccUptakeCohorts      <- vector(mode="list", length = numTimePts+1)
vaccUptakeCohorts[[1]] <- simd %>% add_column(runifA = rep(NA, nrow(.)))
for (i in 0:(numTimePts-1)){
  vaccUptakeTemp <- vaccUptakeCohorts[[1]] %>%
    mutate(prA = plogis(propensityFit$x[propensityFit$X=="(Intercept)"] +
        (propensityFit$x[propensityFit$X=="poly(k, 2, raw = TRUE)1"] * i ) +
        (propensityFit$x[propensityFit$X=="poly(k, 2, raw = TRUE)2"] * (i^2) ) +
        (propensityFit$x[propensityFit$X=="age_c"] * L1_c) +
        (propensityFit$x[propensityFit$X=="sex"] * L2) +
        (propensityFit$x[propensityFit$X=="comorbidity"] * L3))) %>%
    mutate(runifA = runif(nrow(.))) %>%
    mutate(A = if_else( (runifA < prA), 1, 0 ),
          eventTime = Inf ) %>% # initialize eventTime
    mutate(vaccTime = if_else(A==1, i, NA_integer_))

  vaccUptakeCohorts[[1]]     <- vaccUptakeTemp %>% filter(A==0)
  vaccUptakeCohorts[[(i+2)]] <- vaccUptakeTemp %>% filter(A==1)
}

### Events
# Unvacced
for (t in 1:tau){      # t indexes calendar time
  intres <- determine_intercept(
    tolerance = .000000001,
    marg      = cfHazard[[1]][t], # P(E_t=1|E_{t-1}=0, vacctime=NA)
    intercept = NULL,
    beta_vec  = betaVec,
    lower_bound = -30, upper_bound = 30,
    external_dataset = subset(vaccUptakeCohorts[[1]],
        select = c(intercept, L1_c, L2, L3)),
    printout = FALSE)
  
  vaccUptakeCohorts[[1]] <- vaccUptakeCohorts[[1]] %>%
    add_column(intres = rep(intres$intercept, nrow(.)) ) %>%
    mutate(condProb = plogis(
      intres +
        (outcomeFit$x[outcomeFit$X=="age_c"] * L1_c) +
        (outcomeFit$x[outcomeFit$X=="sex"] * L2)  +
        (outcomeFit$x[outcomeFit$X=="comorbidity"] * L3))) %>%
    mutate(runifE = runif(nrow(.))) %>%
    mutate(E = if_else( runifE < condProb , 1, 0 )) %>%
    mutate(eventTime = case_when(
             eventTime < Inf ~ eventTime, # already had event prior to t
             E == 1 ~ as.numeric(t),      # event at time t
             TRUE ~ Inf))                 # no event yet as of t
           
  vaccUptakeCohorts[[1]] <- vaccUptakeCohorts[[1]] %>%
    select(-c("intres", "E", "runifE", "condProb"))
}
for (j in 0:tau){ # j is calendar time of vacc uptake
  for (t in 1:tau){ # t also indexes calendar time
    intres <- determine_intercept(
      tolerance = .0000000000001,
      marg = cfHazard[[j+2]][t], # P(E_t=1|E_{t-1}=0, vacctime=j)
      intercept = NULL,
      beta_vec = betaVec,
      lower_bound = -30, upper_bound = 30,
      external_dataset = subset(vaccUptakeCohorts[[j+2]],
                                select = c(intercept, L1_c, L2, L3)),
      printout = FALSE)
    
    vaccUptakeCohorts[[j+2]] <- vaccUptakeCohorts[[j+2]] %>%
      add_column(intres = rep(intres$intercept, nrow(.)) ) %>%
      mutate(condProb = plogis(
        intres +
          (outcomeFit$x[outcomeFit$X=="age_c"] * L1_c) +
          (outcomeFit$x[outcomeFit$X=="sex"] * L2)  +
          (outcomeFit$x[outcomeFit$X=="comorbidity"] * L3)))  %>%
      
      mutate(runifE = runif(nrow(.))) %>%
      mutate(E = if_else( runifE < condProb , 1, 0 )) %>%
      mutate(eventTime = case_when(
               eventTime < Inf ~ eventTime, # already had event prior to t
               E == 1 ~ as.numeric(t),      # event at time t
               TRUE ~ Inf))                 # no event yet as of t

    vaccUptakeCohorts[[j+2]] <- vaccUptakeCohorts[[j+2]] %>%
      select(-c("intres", "E", "runifE", "condProb"))
  }
}

trialData <- vector(mode="list", length = numTrials)
for (tl in 0:J){ # tl is trial number
  # Trial cohort is those who did not yet get vacc and did not yet have event

  # index 1 is unvacced, remaining indices for vacced at start of `tl` or later
  vTimeIndeces_t <- c(tl:tau) + 2
  vTimeIndeces   <- c(1, vTimeIndeces_t)

  # did not yet have an event
  trialCohort <- lapply(vaccUptakeCohorts[vTimeIndeces],
                        function(x) filter(x, eventTime > tl))

  trialData[[tl+1]] <- bind_rows(trialCohort) %>%
    mutate(Z = case_when(vaccTime == tl  ~ 1,
                    is.na(vaccTime) ~ 0,
                    TRUE ~ 0),
      # time of censoring and event on trial time scale
      # Assuming only form of censoring is going off regimen
      # Note `trialCohort` contains no one with vaccTime < tl, by definition above
      cTime = case_when(is.na(vaccTime) ~ Inf,
        vaccTime == tl   ~ Inf,
        vaccTime > tl    ~ as.numeric(vaccTime-tl)) ,
      eTime = (eventTime - tl),
      trial = tl)

  for (time in tl:numTimePts){
    k <- time - tl # Shift to trial time scale
    trialData[[(tl+1)]]     <- trialData[[(tl+1)]] %>%
      mutate("E_{k}"        := if_else(eTime > k, 0, 1),
            "C_{k}"         := if_else(cTime > k, 0, 1),
            "oneMinusE_{k}" := if_else(eTime > k, 1, 0),
            "oneMinusC_{k}" := if_else(cTime > k, 1, 0),
            "Y_{k}" := if_else( ( (cTime >= k) & (eTime >= k) ), 1, 0), # @risk indicator
            "k_{k}"         := k,
            "kSqr_{k}"      := k^2,
            "Z_k_{k}"       := Z*k,
            "Z_kSqr_{k}"    := Z*k^2)
  }
}
###########################################################################
###
###   Prepare a version of nested trial-specific data for delicatessen
###
###########################################################################
d <- bind_rows(vaccUptakeCohorts) %>% arrange(id)
for (j in 0:J){
  
  d <- d %>% mutate(
    "j_{j}"    := j,
    "jSqr_{j}" := j**2,
    "G_{j}" := case_when( ( (vaccTime >= j) & (eventTime > j) ) ~ 1,
                          ( (is.na(vaccTime)) & (eventTime > j) ) ~ 1,
                          ( (is.na(vaccTime)) & (eventTime <= j) ) ~ 0,
                          TRUE ~ 0),
    "Z_{j}" := case_when(vaccTime <= j   ~ 1,
                         is.na(vaccTime) ~ 0,
                         TRUE ~ 0),
    # time of censoring and event on trial time scale
    # Assuming only form of censoring is going off regimen
    cTime = case_when(is.na(vaccTime) ~ Inf,
                      vaccTime == j   ~ Inf,
                      vaccTime  > j   ~ as.numeric(vaccTime-j),
                      vaccTime  < j   ~ 0,
                      TRUE ~ 0),
    eTime = (eventTime - j))

  for (time in j:numTimePts){
    k <- time - j # Shift to trial time scale
    d <- d %>%
      mutate("E_{j}_{k}"    := if_else(eTime > k, 0, 1),
        "C_{j}_{k}"         := if_else(cTime > k, 0, 1),
        "oneMinusE_{j}_{k}" := if_else(eTime > k, 1, 0),
        "oneMinusC_{j}_{k}" := if_else(cTime > k, 1, 0),
        "Y_{j}_{k}" 
          := if_else( ( (cTime >= k) & (eTime >= k) ), 1, 0), # @risk indicator
        "k_{j}_{k}"         := k,
        "kSqr_{j}_{k}"      := k^2,
        "j_{j}_{k}"         := time,
        "jSqr_{j}_{k}"      := time^2,
        "kPlusTrial_{j}_{k}"    := j + k,
        "kPlusTrialSqr_{j}_{k}" := (j + k)^2,
        "Z_k_{j}_{k}"       := eval(parse(text=paste("Z_", j, sep = ""))) *k,
        "Z_kSqr_{j}_{k}"    := eval(parse(text=paste("Z_", j, sep = ""))) *k^2,
        "Z_j_{j}_{k}"       := eval(parse(text=paste("Z_", j, sep = ""))) *j,
        "Z_jSqr_{j}_{k}"    := eval(parse(text=paste("Z_", j, sep = ""))) *(j^2),
        "Z_kPlusTrial_{j}_{k}"    := eval(parse(text=paste("Z_", j, sep = ""))) *(j+k),
        "Z_kPlusTrialSqr_{j}_{k}" := eval(parse(text=paste("Z_", j, sep = ""))) *((j+k)^2)
        )
  }
}
