# Justin DeMonte  
# 240305
#
# Initial data preparation
library("tidyverse")

## collect args passed in from SLURM job script
args        <- commandArgs(trailingOnly=TRUE)
numTrials   <- as.numeric(args[1])
# end collect args passed in from SLURM job script
readin <- haven::read_dta(file="/nas/longleaf/home/demjus01/ItalyData/Database clean.dta")

### end boilerplate
elderData <- readin 

elderData <- elderData %>%
  mutate(comorbidity = if_else( 
    (cancer + copd + cvd + diabetes + hypertension + kidneydisease) >= 1, 1, 0), 
      datedose1 = datedose1plus14days - 14,
      datedose2 = datedose2plus14days - 14,
      datedose3 = datedose3plus14days - 14,
      baselineAge = as.numeric((parse_date("2021-02-08")-birthdate)/365.25)) %>% 
    filter(baselineAge > 79) %>% # handled more accurately on per-trial basis below
    mutate(meanBLage = mean(baselineAge)) %>%
    ungroup() %>%
    mutate(BLage_c = baselineAge - meanBLage)
##----------------------------------------------------------------
## Define study period
##----------------------------------------------------------------
startDate <- as.Date("2021-02-15") 

# admin censoring date for our analysis
# Defined as the time where we no longer have new vaccine uptake information
endDate   <- as.Date("2021-12-18")

# Determine numTimePts as a function of startDate and endDate
numTimePts   <- 0
start        <- startDate
while (start <= endDate) {
  numTimePts <- numTimePts + 1
  start      <- start + 7  
}
# Restrict dataset to those who are eligible for first trial
# i.e., trial initiated at `startDate`.
elderData <- elderData %>% 
  filter(is.na(dateinfection) | dateinfection > startDate) %>% 
  filter(is.na(datedeath)     | datedeath >  startDate) %>% 
  filter(is.na(datecovid)     | datecovid >  startDate) %>% 
  filter(is.na(datedose1)     | datedose1 >= startDate) %>%
  filter(is.na(datedose2)     | datedose2 >= startDate) %>%
  filter(is.na(datedose3)     | datedose3 >= startDate) %>%
  
  # cannot identify time at which 1,079 (~1%) individuals initiated treatment
  # excluding them 
  # as an aside, none of these individuals have severe covid 
  # or covid-related death during the study period.
  filter((is.na(datedose1) & dosenumber==0) | !is.na(datedose1))

##----------------------------------------------------------------
## Derive event and dropout times
##----------------------------------------------------------------
# initial values for derived variables
elderData <- elderData %>% mutate(
  
  eventTime      = Inf,
  dropoutTime    = Inf,
  covidDeathTime = Inf,
  infectionTime  = Inf,
  covidTime      = Inf,
  deathTime      = Inf,
  turn80time     = Inf,
  dose3time      = Inf,
  dropoutDate    = if_else((!is.na(datedeath) & death==0), datedeath, as.Date(NA_character_)),
  covidDeathDate = if_else((!is.na(datedeath) & death==1), datedeath, as.Date(NA_character_))
)

start <- startDate
j     <- 0
while (start <= endDate){
  end <- start + 6
    
  elderData <- elderData %>% mutate(turn80time = case_when( 
    turn80time < Inf ~ turn80time, # already turned 80 prior to j
    ((start - birthdate)/365.25 > 80) ~ as.numeric(j), # turned 80 at time j
    TRUE ~ Inf)) %>% # did not yet turn 80 as of time j
    
    # for comparing strategies "2 doses to no doses" 
    # if someone gets a 3rd dose, censor them 
    # following Acuti Martellucci definition of 3rd dose by vacc type
    mutate(dose3time = case_when(
      dose3time < Inf ~ dose3time,   # already had dose3 prior to j
      (start <= datedose3 & datedose3 <= end) ~ as.numeric(j+1), # dose3 prior to j+1
      TRUE ~ Inf)) %>%
    
    # outcomes 
    mutate(infectionTime = case_when(
      infectionTime < Inf ~ infectionTime, # already had infection prior to j
      (start <= dateinfection & dateinfection <= end) ~ as.numeric(j+1), # infection at time j+1
      TRUE ~ Inf)) %>%             # no infection yet as of time j 
    
    mutate(covidTime = case_when(
      covidTime < Inf ~ covidTime, # already had severe covid prior to j
      (start <= datecovid & datecovid <= end) ~ as.numeric(j+1), # severe covid at time j+1
      TRUE ~ Inf)) %>%             # no severe covid yet as of time j 
    
    # This `deathTime` variable is all-cause
    mutate(deathTime = case_when(
      deathTime < Inf ~ deathTime, # already had death prior to j+1
      (start <= datedeath & datedeath <= end) ~ as.numeric(j+1), # death at time j+1
      TRUE ~ Inf)) %>%             # no death yet as of time j+1 
    
    mutate(covidDeathTime = case_when(
      covidDeathTime < Inf ~ covidDeathTime, # already had death from COVID prior to j
      (start <= covidDeathDate & covidDeathDate <= end) ~ as.numeric(j+1), # death from COVID at time j+1
      TRUE ~ Inf)) %>%             # death from COVID at time j+1 yet as of time j
  
    # dropout due to death from nonCOVID cause 
    mutate(dropoutTime = case_when(
      dropoutTime < Inf ~ dropoutTime, # already dropped out prior to j
      (start <= dropoutDate & dropoutDate <= end) ~ as.numeric(j+1), # dropout at time j+1
      TRUE ~ Inf))                 # did not yet dropout as of time j
      
  elderData   <- elderData %>% mutate(eventTime = case_when(
    eventTime < Inf ~ eventTime, # already had event prior to j
    ( (start <= datecovid & datecovid <= end) | # event at time j+1
    start <= covidDeathDate & covidDeathDate <= end ) ~ as.numeric(j+1),
    TRUE ~ Inf))               # no event yet as of time j
      
  j     <- j + 1
  start <- start + 7   
}
##----------------------------------------------------------------
## Create vaccination-time-specific cohorts 
##----------------------------------------------------------------
# vaccUptakeCohorts[[1]] holds those who remained unvaccinated at end of study
# vaccUptakeCohorts[[j+2]] holds those who initiated vaccination at time j
vaccUptakeCohorts      <- vector(mode="list", length = numTimePts+1)
vaccUptakeCohorts[[1]] <- elderData
start                  <- startDate
j                      <- 0
while (start <= endDate) {
  end                  <- start + 6
  vaccUptakeTemp       <- vaccUptakeCohorts[[1]] %>% 
    mutate( A = case_when ( # A indicates vaccine uptake during jth interval
      ( start <= datedose1 & datedose1 <= end ) ~ 1, 
      is.na(datedose1) ~ 0, 
      TRUE ~ 0 )) %>%
    mutate(vaccTime = if_else(A==1, j, NA_integer_)) 
  
  vaccUptakeCohorts[[1]]   <- vaccUptakeTemp %>% filter(A==0) %>% select(-c("A"))
  vaccUptakeCohorts[[j+2]] <- vaccUptakeTemp %>% filter(A==1) %>% select(-c("A"))
  j     <- j + 1
  start <- start + 7    
}
##----------------------------------------------------------------
## Calculate distribution of vaccine type in analytic cohort
##----------------------------------------------------------------
numVaccUptkCohorts <- length(vaccUptakeCohorts)
trialJindex        <- numTrials + 1

for (df in 2:trialJindex){ # index 1 corresponds to never uptake vacc
  vaccUptakeCohorts[[df]] <- vaccUptakeCohorts[[df]] %>%
    mutate(completeReg = case_when( 
      is.na(datedose2) ~ 0,
      vaccinetypedose1 == 4 ~ 1,
      ( vaccinetypedose1 != 3 & (datedose2-datedose1 <= 42) ) ~ 1, 
      ( vaccinetypedose1 == 3 & (datedose2-datedose1 <= 84) ) ~ 1,
      TRUE ~ 0) ) %>%
    # vType 5 denotes mixed 1st and 2nd doses
    mutate(vType = case_when(
      completeReg == 0 ~ vaccinetypedose1,
      ( (completeReg == 1) & (vaccinetypedose1 != 4) & (vaccinetypedose1 != vaccinetypedose2) ) ~ 5,
      ( (completeReg == 1) & vaccinetypedose1 == 4) ~ 4,
      TRUE ~ vaccinetypedose1)) 
  print(vaccUptakeCohorts[[df]] %>% group_by(completeReg, vType) %>% count())
}

vaccUptakeTempList <- vaccUptakeCohorts[-1]
vaccUptakeTemp <- bind_rows(vaccUptakeTempList)
print(vaccUptakeTemp %>% group_by(completeReg, vType) %>% count())
vaccUptakeTemp %>% filter(completeReg==1 & vType==4) %>%
  select(datedose1, datedose2, datedose3, vaccinetypedose1, vaccinetypedose2, vaccinetypedose3) %>% 
  print(width=Inf)
vaccUptakeTemp %>% filter(completeReg==0 & vType==4) %>%
  select(datedose1, datedose2, datedose3, vaccinetypedose1, vaccinetypedose2, vaccinetypedose3) %>%
  print(width=Inf)
##----------------------------------------------------------------
## Save permanent datasets for use downstream
##----------------------------------------------------------------
dataDir <- "data/"

dfNames <- vector(length = numVaccUptkCohorts)
for (df in 1:numVaccUptkCohorts){
  dfNames[df] <- paste0(dataDir, "vaccUptakeCohorts", df)  
}
names(vaccUptakeCohorts) <- dfNames

# save a separate dataset for each vacc uptake cohort
lapply(names(vaccUptakeCohorts), function(x) {
  x1 <- vaccUptakeCohorts[[x]]
  saveRDS(x1, file=paste0(x, '.rds'))
})
# pass constants from SLURM script on to downstream Rscripts.
saveRDS(numTimePts,  file = paste0(dataDir, 'numTimePts.rds'))
saveRDS(numTrials,   file = paste0(dataDir, "numTrials.rds"))