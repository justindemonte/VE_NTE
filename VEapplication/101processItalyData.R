# 101processItalyData.R
#
# input: 'raw' Abruzzo data 
# /nas/longleaf/home/demjus01/ItalyData/Database clean.dta
#
# args: numTrials = number of trials to be emulated (i.e., J+1)
# 
# output: 1. analytic cohort data (as described in Section 2) in wide format 
# elderdata.rds 
# 
# 2. numTrials - a constant (arg passed to this program)
# /nas/longleaf/home/demjus01/VEapplication/data/numTrials.rds
# 
# 3. numTimePoints - a constant 
# /nas/longleaf/home/demjus01/VEapplication/data/numTimePts.rds

library("tidyverse")
# collect args passed in from SLURM job script
args      <- commandArgs(trailingOnly=TRUE)
numTrials <- as.numeric(args[1])
# end collect args passed in from SLURM job script

J         <- numTrials - 1

readin    <- haven::read_dta(file="/nas/longleaf/home/demjus01/ItalyData/Database clean.dta")

minAge    <- 80

pseudoInf <- 999

# define constants 
kminPfizer  <- 3
kminModerna <- 4
kminAstraZ  <- 8

kmaxPfizer  <- 6
kmaxModerna <- 6
kmaxAstraZ  <- 12
### end boilerplate
elderData <- readin 

elderData <- elderData %>% mutate(X3 = if_else( 
    (cancer + copd + cvd + diabetes + hypertension + kidneydisease) >= 1, 1, 0), 
      X2        = sex, 
      datedose1 = datedose1plus14days - 14,
      datedose2 = datedose2plus14days - 14,
      datedose3 = datedose3plus14days - 14,
      X1 = as.numeric((parse_date("2021-02-15")-birthdate)/365.25)) %>% 
    filter(X1 >= minAge)      %>% 
    mutate(meanX1 = mean(X1)) %>%
    mutate(X1_c = X1 - meanX1)
##----------------------------------------------------------------
## Define study period
##----------------------------------------------------------------
startDate <- as.Date("2021-02-15") 

# admin censoring date for Abruzzo data analysis
# Defined as the time where we no longer have new vaccine uptake information
endDate   <- as.Date("2021-12-18")

# Determine numTimePts as a function of startDate and endDate
numTimePts   <- 1
start        <- startDate
while (start <= endDate) {
  numTimePts <- numTimePts + 1
  start      <- start + 7  
}
tau          <- numTimePts-1

print(paste0("numTimePts: ", numTimePts))
print(paste0("tau: ", tau))

# Narrow dataset to those who meet certain eligibility crit for first trial
elderData <- elderData %>% 
  filter(is.na(dateinfection) | dateinfection >= startDate) %>% 
  filter(is.na(datedeath)     | datedeath >= startDate)     %>% 
  filter(is.na(datecovid)     | datecovid >= startDate)     %>% 
  filter(is.na(datedose1)     | datedose1 >= startDate)     %>%
  filter(is.na(datedose2)     | datedose2 >= startDate)     %>%
  filter(is.na(datedose3)     | datedose3 >= startDate) 
##----------------------------------------------------------------
## Data cleaning
##----------------------------------------------------------------
# cannot identify time at which 972 (<1%) individuals initiated treatment
# excluding them 
# as an aside, none of these individuals have severe covid 
# or covid-related death during the study period.
print("number of rows before dropping those with unknown date of first dose")
elderData %>% count()

elderData <- elderData %>% 
  filter( ( is.na(datedose1) & dosenumber==0 ) | !is.na(datedose1) )

print("number of rows after dropping those with unknown date of first dose")
print("Should differ by 972.")
elderData %>% count()

# In 'raw' data, people with 1 dose of Janssen have datedose1==datedose2
# Recode these so that datedose2 is NA to comport with coding scheme below
print("number of rows with datedose2=NA before recoding Janssen")
elderData %>% filter(is.na(datedose2)) %>% count() 

elderData <- elderData %>% mutate(datedose2 = 
  if_else( ( datedose1==datedose2 & vaccinetypedose1==4 ), NA, datedose2))

print("number of rows with datedose2=NA after recoding Janssen")
print("Should differ by 167.")
elderData %>% filter(is.na(datedose2)) %>% count() 

# In 'raw' data, vaccinetypedose2 is occasionally NA when datedose2 != NA.
# All of these cases involve a first dose in Spring 2021 and a second dose in fall 2021
# Hence in our analysis, these will all result in censoring.  
# Need to set vaccinetypedose2 to an arbitrary non-NA value to comport with coding scheme below
print("number of rows with datedoseX not missing and vaccinetypedoseX missing before recoding")
elderData %>% filter(
  ( is.na(vaccinetypedose1) & !is.na(datedose1) ) | ( is.na(vaccinetypedose2) & !is.na(datedose2) ) | ( is.na(vaccinetypedose3) & !is.na(datedose3) )
  ) %>% count() 

elderData <- elderData %>% mutate(vaccinetypedose2 = 
  if_else( ( is.na(vaccinetypedose2) & !is.na(datedose2) ), vaccinetypedose3, vaccinetypedose2))

print("number of rows with datedoseX not missing and vaccinetypedoseX missing after recoding")
print("Should be zero.")
elderData %>% filter(
  ( is.na(vaccinetypedose1) & !is.na(datedose1) ) | ( is.na(vaccinetypedose2) & !is.na(datedose2) ) | ( is.na(vaccinetypedose3) & !is.na(datedose3) )
) %>% count() 

print("Identify records with misordered vaccination dates")
elderData %>% filter(
  (datedose2 <= datedose1) | (datedose3 <= datedose1) | (datedose3 <= datedose2)
) %>% print()

print("number of rows before removing person with misordered vaccination dates")
elderData %>% count()
elderData <- elderData %>% filter(id!=895895) 
print("number of rows after removing person with misordered vaccination dates")
print("should differ by 1")
elderData %>% count()
##----------------------------------------------------------------
## Derive event and dropout times
##----------------------------------------------------------------
# initial values for derived variables
elderData <- elderData %>% mutate(
  
  T_             = pseudoInf, # calendar time of event, i.e., T in manuscript
  U              = pseudoInf,
  V_1            = pseudoInf, # time of first vaccine dose
  V_2            = pseudoInf, # time of second vaccine dose
  
  dropoutDate    = if_else((!is.na(datedeath) & death==0), datedeath, as.Date(NA_character_)),
  covidDeathDate = if_else((!is.na(datedeath) & death==1), datedeath, as.Date(NA_character_))
)
elderData  <- elderData %>% mutate(Z_0 = 0) 

start <- startDate
l     <- 0
while (start <= endDate){
  end        <- start + 6
  lPlus1     <- l + 1 
  print(start)
  print(end)
  print(l)
  print(lPlus1)
  
  print(paste0("numRows for which datedoseX and datedoseY both occur within interval ", l))
  print(
    elderData %>% filter(
      ( (start <= datedose1 & datedose1 <= end) & (start <= datedose2 & datedose2 <= end) ) |  
      ( (start <= datedose2 & datedose2 <= end) & (start <= datedose3 & datedose3 <= end) ) |  
      ( (start <= datedose1 & datedose1 <= end) & (start <= datedose3 & datedose3 <= end) )  
    ) %>% count()
  )
  
  elderData  <- elderData %>% mutate(U  = case_when(
  
    # dropout due to death from nonCOVID cause 
    U < pseudoInf ~ U,          # already dropped out prior to l
    (start <= dropoutDate & dropoutDate <= end) ~ as.numeric(l+1), # dropout at time l+1
    TRUE ~ pseudoInf)) %>%      # did not yet dropout as of time l
      
    # event time
    mutate(T_ = case_when(
      T_ < pseudoInf ~ T_,        # already had event prior to l
      ( (start <= datecovid & datecovid <= end) | # event at time l+1
      ( start <= covidDeathDate & covidDeathDate <= end ) ) ~ as.numeric(l+1),
      TRUE ~ pseudoInf)) %>%      # no event yet as of time l
    
    mutate(Tstar = pmin(T_, U),
      Delta = if_else(T_ < U, 1, 0)) %>% 
    
    mutate( "B_{lPlus1}" := case_when(
      (start <= datedose1 & datedose1 <= end) ~ vaccinetypedose1,   # dose1 at time l+1
      (start <= datedose2 & datedose2 <= end) ~ vaccinetypedose2,   # dose2 at time l+1
      (start <= datedose3 & datedose3 <= end) ~ vaccinetypedose3,   # dose3 at time l+1
      TRUE ~ 0)) %>%
    
    mutate(V_1 = case_when(
      V_1    <  pseudoInf ~ V_1,  # already had dose1 prior to l
      (start <= datedose1 & datedose1 <= end) ~ as.numeric( l+1 ), # received dose1 at time l+1
      TRUE ~ pseudoInf)) %>%
    
    mutate(V_2 = case_when(
      V_2    <  pseudoInf ~ V_2,  # already had dose2 prior to l
      (start <= datedose2 & datedose2 <= end) ~ as.numeric( l+1 ), # received dose2 at time l+1
      TRUE ~ pseudoInf)) 
  
  B_lplus1.nm <- paste0("B_", lPlus1)
  Z_l.nm      <- paste0("Z_", l     )
  
  # derive Z_l
  elderData <- elderData %>% mutate( "Z_{lPlus1}" := case_when(
    ( get({{B_lplus1.nm}}) == 0) & ( get({{Z_l.nm}}) == 0) & (Tstar > l) ~ 0,  
    
    ( get({{B_lplus1.nm}}) == 0) & ( get({{Z_l.nm}}) == 1) & (Tstar > l) ~ 1,
    ( get({{B_lplus1.nm}}) %in% c(1, 2, 3, 4) ) & ( get({{Z_l.nm}}) == 0) & (Tstar > l) ~ 1,
    
    ( get({{B_lplus1.nm}}) == 0) & ( get({{Z_l.nm}}) == 2) & (Tstar > l) ~ 2,
    ( get({{B_lplus1.nm}}) %in% c(1, 2, 3, 4) ) & ( get({{Z_l.nm}}) == 1) & (Tstar > l) ~ 2,
    
    ( get({{Z_l.nm}}) == 3) & (Tstar > l) ~ 3,
    ( get({{B_lplus1.nm}}) %in% c(1, 2, 3, 4) ) & ( get({{Z_l.nm}}) == 2) & (Tstar > l) ~ 3,
    
    ( Tstar <= l) ~ get({{Z_l.nm}}), # let previous value ride to avoid NAs in the dataset
    TRUE ~ pseudoInf))  
  
  # derive eligibility indicators for each calendar time point
  elderData[ paste0("E_", l) ] <- ( elderData["V_1"] > l ) *  # aka Z_l=0
    ( elderData["Tstar"] > l ) # free of event and LTFU at time l

  l     <- l + 1
  start <- start + 7   
}
print("Number of individuals in the analysis dataset")
elderData %>% count()

print("table of U")
print(table(elderData$U))
print("table of T")
print(table(elderData$T_))
print("table of Tstar")
print(table(elderData$Tstar))
print("table of Delta")
print(table(elderData$Delta))
print("table of V_1")
print(table(elderData$V_1))
##----------------------------------------------------------------
##  Descriptive analyses for first paragraph of data analysis results
##----------------------------------------------------------------
print("Final N for analytic cohort - number of indv'ls who meet elg crit for trial zero")
elderData %>% count()

print("Number of people who received a first COVID-19 vaccine dose by May 9")
elderData %>% filter(V_1 <= 12) %>% count()

print("Number of people who completed Phizer regimen as directed")
elderData %>% mutate(sojourn = V_2-V_1) %>%
  filter( (V_1 <= 12) & vaccinetypedose1==1 & 
    sojourn>=kminPfizer & sojourn<=kmaxPfizer) %>% count()

print("Number of people who completed Moderna regimen as directed")
elderData %>% mutate(sojourn = V_2-V_1) %>%
  filter( (V_1 <= 12) & vaccinetypedose1==2 & 
    sojourn>=kminModerna & sojourn<=kmaxModerna) %>% count()

print("Number of people who completed AstraZ regimen as directed")
elderData %>% mutate(sojourn = V_2-V_1) %>%
  filter( (V_1 <= 12) & vaccinetypedose1==3 & 
    sojourn>=kminAstraZ & sojourn<=kmaxAstraZ) %>% count()

print("Number of people who completed JandJ regimen as directed")
elderData %>% mutate(sojourn = V_2-V_1) %>%
  filter( (V_1 <= 12) & vaccinetypedose1==4) %>% count()

elderData <- elderData %>% select(starts_with("E_"),
                        starts_with("Z_"),
                        id, X1_c, X2, X3,
                        V_1, vaccinetypedose1,
                        Tstar, Delta) %>% haven::zap_labels()

dataDir <- "data/"

# save a permanent copy of analytic cohort data in wide format
saveRDS(elderData,  file = paste0(dataDir, 'elderData.rds'))

# pass constants from SLURM script on to downstream Rscripts.
saveRDS(numTimePts, file = paste0(dataDir, 'numTimePts.rds'))
saveRDS(numTrials,  file = paste0(dataDir, "numTrials.rds"))