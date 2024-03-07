# Justin DeMonte  
# 240305
#
library("tidyverse")
##----------------------------------------------------------------
## Read in permanent datasets 
##----------------------------------------------------------------
dataDir <- "data/"

numTimePts  <- readRDS(file=paste0(dataDir, "numTimePts.rds"))
numTrials   <- readRDS(file=paste0(dataDir, "numTrials.rds"))

vaccUptakeCohorts <- vector(mode="list", length = numTimePts+1)
for (trl in 1:(numTimePts+1)){
  vaccUptakeCohorts[[trl]] <- readRDS(file=paste0(dataDir, "vaccUptakeCohorts", trl, ".rds"))
}
##----------------------------------------------------------------
## Prepare a version of nested-trial-specific data for delicatessen
##----------------------------------------------------------------
d <- bind_rows(vaccUptakeCohorts) %>% arrange(id) %>% 
  haven::zap_labels() # remove variable formatting from Stata dataset
for (j in 0:(numTimePts-1)){
    
  d <- d %>% mutate(
    "G_{j}" := case_when( # eligibility indicator
      ( (vaccTime >= j)   & (covidTime > j & deathTime > j & 
          turn80time <= j & infectionTime > j) ) ~ 1,
      ( (is.na(vaccTime)) & (covidTime > j & deathTime > j & 
          turn80time <= j & infectionTime > j) ) ~ 1,
      ( (is.na(vaccTime)) & (covidTime <= j | deathTime <= j | 
          turn80time > j | infectionTime <= j) ) ~ 0,
      TRUE ~ 0), 
    
    "Z_{j}" := case_when(
      vaccTime <= j   ~ 1,
      is.na(vaccTime) ~ 0,
      TRUE            ~ 0),
    
    # time of censoring and event on trial time scale 
    cTime1 = case_when( # going off regimen 0 by getting a first dose
      is.na(vaccTime) ~ Inf,
      vaccTime == j   ~ Inf, # received regimen 1 
      vaccTime >  j   ~ as.numeric(vaccTime - j), # goes off regimen 0 at vaccTime-j
      vaccTime <  j   ~ 0,   # not eligible for trial j
      TRUE ~ 0),
    # censor those who get enrolled in the treatment arm of this trial
    # and fail to get a second dose within recommended window
    # i.e., within 6 weeks for Moderna and Phizer
    # within 12 weeks for Astra Zeneca
    # If JandJ (i.e., vaccinetypedose1==4) never censor because no need for second dose.
    cTime2 = case_when( 
      ( vaccTime == j & vaccinetypedose1 != 4 & vaccinetypedose1 != 3 & (datedose2-datedose1 > 42) ) ~ 6, 
      ( vaccTime == j & vaccinetypedose1 == 3 & (datedose2-datedose1 > 84) ) ~ 12,
      TRUE ~ Inf),
    cTime3 = (dose3time - j),   # get a third dose
    doTime = (dropoutTime - j), # non-COVID death 
    eTime  = (eventTime - j))
  
  d <- d %>% mutate(cTime = pmin(cTime1, cTime2, cTime3, doTime))   
  
  for (time in j:numTimePts){
    k <- time - j # Shift to trial time scale
    d <- d %>%
      mutate(
        "E_{j}_{k}"             := if_else(eTime > k, 0, 1),
        "C_{j}_{k}"             := if_else(cTime > k, 0, 1),
        "oneMinusE_{j}_{k}"     := if_else(eTime > k, 1, 0),
        "oneMinusC_{j}_{k}"     := if_else(cTime > k, 1, 0),
        "Y_{j}_{k}"             := 
          if_else( ( (cTime >= k) & (eTime >= k) ), 1, 0), # @risk indicator
        "k_{j}_{k}"             := k ,
        "kSqr_{j}_{k}"          := k^2,
        "j_{j}_{k}"             := time,
        "jSqr_{j}_{k}"          := time^2,
        "kPlusTrial_{j}_{k}"    := j + k,
        "kPlusTrialSqr_{j}_{k}" := (j + k)^2,
        "Z_k_{j}_{k}"     := eval(parse(text=paste("Z_", j, sep = ""))) *k,
        "Z_kSqr_{j}_{k}"  := eval(parse(text=paste("Z_", j, sep = ""))) *k^2,
        "Z_j_{j}_{k}"     := eval(parse(text=paste("Z_", j, sep = ""))) *j,
        "Z_jSqr_{j}_{k}"  := eval(parse(text=paste("Z_", j, sep = ""))) *(j^2)
      )
  }
}
for (j in 0:(numTrials-1)){
  d <- d %>% mutate(
        "trial_{j}"        := j,
        "trialSquared_{j}" := j**2
    )
}
##----------------------------------------------------------------
## Save permanent dataset for use downstream
##----------------------------------------------------------------
# Remove variable labels read in from STATA format dataset
d <- d %>% select(-c(birthdate, starts_with("date"))) 
saveRDS(d, file=paste0(dataDir, 'd_t.rds'))