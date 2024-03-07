# Justin DeMonte  
# 240305
#
library("tidyverse")

## collect args passed in from SLURM job script
tl          <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# end collect args passed in from SLURM job script

##----------------------------------------------------------------
## Read in permanent datasets 
##----------------------------------------------------------------
dataDir <- "data/"

# read in user-provided constants from upstream program
numTimePts  <- readRDS(file=paste0(dataDir, "numTimePts.rds"))
numTrials   <- readRDS(file=paste0(dataDir, "numTrials.rds"))

vaccUptakeCohorts <- vector(mode="list", length = (numTimePts+1))
for (trl in 1:(numTimePts+1)){
  vaccUptakeCohorts[[trl]] <- readRDS(file=paste0(dataDir, 
    "vaccUptakeCohorts", trl, ".rds"))
}  
##----------------------------------------------------------------
## Build trial cohort dataset
##----------------------------------------------------------------

# Trial cohort consists of those who meet inclusion criteria 
# and did not yet get vacc.
# Recall vaccUptakeCohorts[[1]] holds those who remained unvaccinated 
# throughout entire observation period (through Dec. 2021)
# vaccUptakeCohorts[[j+2]] holds those who initiated treatment at time j. 
vTimeIndeces_t <- c(tl:(numTimePts-1)) + 2
vTimeIndeces   <- c(1, vTimeIndeces_t)

# Did not yet get vacc. and meet inclusion criteria
trialCohort <- lapply(vaccUptakeCohorts[vTimeIndeces], function(x){filter(x, 
  covidTime  >  tl & 
  turn80time <= tl & 
  deathTime  >  tl &
  infectionTime > tl)})

trialData <- bind_rows(trialCohort) %>%
  haven::zap_labels() %>% 
  mutate(
    Z = case_when(
      vaccTime == tl  ~ 1,
      is.na(vaccTime) ~ 0,
      TRUE ~ 0),
    
    # time of censoring and event on trial time scale 
    cTime1 = case_when(
      is.na(vaccTime) ~ Inf,
      vaccTime == tl  ~ Inf,
      vaccTime >  tl  ~ as.numeric(vaccTime - tl)), # goes off regimen at vaccTime-tl
    
    # censor those who get enrolled in the treatment arm of this trial
    # and fail to get a second dose within recommended window
    # i.e., within 6 weeks for Moderna and Phizer
    # within 12 weeks for Astrazeneca
    cTime2 = case_when( 
      ( vaccTime == tl & vaccinetypedose1 != 4 & vaccinetypedose1 != 3 & (datedose2-datedose1 > 42) ) ~ 6, 
      ( vaccTime == tl & vaccinetypedose1 == 3 & (datedose2-datedose1 > 84) ) ~ 12,
       TRUE ~ Inf),
    cTime3 = (dose3time   - tl), # censor at time of 3rd dose
    doTime = (dropoutTime - tl), # trial time of death from non-COVID cause 
    eTime  = (eventTime   - tl),
    trial  = tl) 
  
  trialData <- trialData %>% mutate(cTime = pmin(cTime1, cTime2, cTime3, doTime))

for (time in tl:numTimePts){
  k <- time - tl # Shift to trial time scale
  trialData <- trialData %>%
    mutate(
      "E_{k}"         := if_else(eTime > k, 0, 1),
      "C_{k}"         := if_else(cTime > k, 0, 1),
      "oneMinusE_{k}" := if_else(eTime > k, 1, 0),
      "oneMinusC_{k}" := if_else(cTime > k, 1, 0),
      "Y_{k}" := if_else( ( (cTime >= k) & (eTime >= k) ), 1, 0), # @risk indicator
      "k_{k}"         := k)
}
##----------------------------------------------------------------
## Save permanent dataset for use downstream
##----------------------------------------------------------------
saveRDS(trialData, file = paste0(dataDir, 'trialData', (tl+1), '.rds'))
##----------------------------------------------------------------
## Create tables of trial-specific descriptive stats
##----------------------------------------------------------------
library("xtable")
outputDir <- "output/"

Ztable <- trialData %>%
  group_by(Z) %>% 
  count()

print(xtable(Ztable), type="html", file=paste0(outputDir, "Ztable", (tl+1), ".html"))
