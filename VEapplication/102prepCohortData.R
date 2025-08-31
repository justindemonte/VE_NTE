# 102prepCohortData.R
#
# input: 
# 1. analytic cohort data (as described in Section 2) in wide format 
# /nas/longleaf/home/demjus01/VEapplication/data/elderdata.rds 
#
# 2. numTrials, a constant, number of trials to be emulated(i.e., J-1)
# /nas/longleaf/home/demjus01/VEapplication/data/numTrials.rds
# 
# 3. numTimePoints - a constant 
# /nas/longleaf/home/demjus01/VEapplication/data/numTimePts.rds
#
# output: 
# 1. analytic cohort data (as described in Section 2) in long format 
# /nas/longleaf/home/demjus01/VEapplication/data/cohortDataLong.rds

library("tidyverse")

dataDir    <- "data/"
# read in user-provided constants from upstream program
numTimePts <- readRDS(file=paste0(dataDir, "numTimePts.rds"))
numTrials  <- readRDS(file=paste0(dataDir, "numTrials.rds"))

tau        <- numTimePts - 1
J          <- numTrials  - 1

cohortDataWide <- readRDS(file=paste0(dataDir, "elderData.rds"))

# define constants, open and close of brand-specific grace periods 
kminPfizer  <- 3
kminModerna <- 4
kminAstraZ  <- 8

kmaxPfizer  <- 6
kmaxModerna <- 6
kmaxAstraZ  <- 12
# end boilerplate
###########################################################################
###
###   Create long format data (one row per calendar time) for observational cohort
###   
###########################################################################
for ( l in 0:tau ){
  cohortDataWide <- cohortDataWide %>% 
    mutate( "Y_{l}" := if_else( (Tstar <= l & Delta == 1), 1, 0 ),
            "H_{l}" := if_else( (Tstar <= l & Delta == 0), 1, 0 )) # censor due to LTFU
}
cohortDataLong <- cohortDataWide %>% 
  mutate( X1X3 = X1_c * X3)      %>% # mean-centered baseline age by comorbid interaction
  # censor due to LTFU
  select( c( "id", "X1_c", "X2", "X3", "X1X3", 
             "V_1", "vaccinetypedose1",
             num_range( "E_",    0:tau ),
             num_range( "Y_",    0:tau ),
             num_range( "H_",    0:tau ),
             num_range( "Z_",    0:tau ))) %>%
  pivot_longer(
    cols       = starts_with( c( "Y_", "H_", "Z_") ),
    names_to   = c( ".value", "time" ),
    names_sep  = "\\_" )  %>%
  mutate( l    = as.numeric(time) ) %>%
  group_by( id ) %>% 
  mutate( lagY = dplyr::lag(Y, default = 0), 
          lagH = dplyr::lag(H, default = 0),
          lagZ = dplyr::lag(Z, default = 0) ) %>%
  ungroup() %>%
  mutate( D = Z - lagZ) %>%   # increment to vaccine uptake process
            
  # derive dummy variables for Z to be used in LTFU model downstream
  mutate(Zsup1  = if_else( Z=="1", 1, 0 ),  # dummy variable for Z=1
         Zsup2  = if_else( Z=="2", 1, 0 ),  # dummy variable for Z=2
         Zsup3  = if_else( Z=="3", 1, 0 )   # dummy variable for Z=3
         ) %>%
  mutate( # for use in linear pred of vacc model
    IlagZeq0            = if_else( lagZ == 0, 1, 0 ), 
    IlagZeq1            = if_else( lagZ == 1, 1, 0 ), 
    IlagZeq2            = if_else( lagZ == 2, 1, 0 ), 
    
    # Indicators for pre-grace period time points by brand
    IpreGrace_Pfizer    = if_else( (vaccinetypedose1==1) & (l > V_1) & (l < (V_1 + kminPfizer)),  1, 0),
    IpreGrace_Moderna   = if_else( (vaccinetypedose1==2) & (l > V_1) & (l < (V_1 + kminModerna)), 1, 0),
    IpreGrace_AstraZ    = if_else( (vaccinetypedose1==3) & (l > V_1) & (l < (V_1 + kmaxAstraZ)),  1, 0),
    
    # Indicators for time points inside grace period by brand
    Igrace_mRNA    = if_else( ( (vaccinetypedose1==1) & (l >= (V_1 + kminPfizer))  & (l < (V_1 + kmaxPfizer))  ) |
                              ( (vaccinetypedose1==2) & (l >= (V_1 + kminModerna)) & (l < (V_1 + kmaxModerna)) ), 
                              1, 0),
    
    Igrace_AstraZ  = if_else( (vaccinetypedose1==3) & (l >= (V_1 + kminAstraZ)) & (l < (V_1 + kmaxAstraZ)), 1, 0),
    
    # Generic indicators for time points inside grace period
    IstateOfGrace  = ifelse( Igrace_mRNA==1 | Igrace_AstraZ==1, 1, 0),
    
    # Indicator for time point at which grace period closes
    IfallFromGrace = if_else( ( (vaccinetypedose1==1) & (l == (V_1 + kmaxPfizer)) ) | 
                              ( (vaccinetypedose1==2) & (l == (V_1 + kmaxModerna))) | 
                              ( (vaccinetypedose1==3) & (l == (V_1 + kmaxAstraZ)) ), 
                              1, 0), 
    ##----------------------------------------------------------------
    ## build 'at-risk' indicator for vaccine uptake process model (VUPM)
    ##----------------------------------------------------------------
    ### * Mark the first row where censoring occurs due to early dose with a 1
    secondDose2early = if_else( (Z==2) & 
                         (IpreGrace_Pfizer==1 | 
                         IpreGrace_Moderna==1 | 
                         IpreGrace_AstraZ==1), 
                         1, 0),
    
    secondDose2late  = if_else( ( (vaccinetypedose1==1) & (Z==1) & (l >= (V_1 + kmaxPfizer))) | 
                         ( (vaccinetypedose1==2) & (Z==1) & (l >= (V_1 + kmaxModerna) )) | 
                         ( (vaccinetypedose1==3) & (Z==1) & (l >= (V_1 + kmaxAstraZ) ) ),
                         1, 0))

cohortDataLong <- cohortDataLong %>% 
  group_by(id) %>%
  
  ### * Now fill in 1's for the censoring indicator for all time points after censoring occurs
  mutate( secondDose2early = if_else(cumsum(secondDose2early) > 0, 1, 0 )) %>%
  mutate( secondDose2late  = if_else(cumsum(secondDose2late)  > 0, 1, 0 )) %>%
  ungroup() %>%
  mutate(
    ### three components of 'at-risk' indicator for VUPM:
    
    # 1. can't get booster before booster rollout
    # Note l=32 aligns with booster dose rollout for the Abruzzo data
    IlLtl_23.FV         = if_else( 
      ( (vaccinetypedose1 %in% c(1,2,3)) & (lagZ == 2) & (l < 32)) |
      ( (vaccinetypedose1==4) & (lagZ == 1) & (l < 32)), 
      1, 0 ), 
    # 2. on-study 
    ITstarGtlminus1     = if_else( ( lagH == 0 ) & ( lagY == 0 ), 1, 0 ),
    # 3. I(nonadherence to A_{V_1})=0
    # composite indicator for nonadherence to A_{V_1}
    InonAdherentToA_v1  = if_else( (vaccinetypedose1==4 & Z==2) | # J and J
                                   Z==3 | 
                                   secondDose2early | 
                                   secondDose2late, 
                                   1, 0),
    
    # Composite at-risk indicator for VUPM (referred to as R^Z in supplement): 
    RsupD  = if_else( ITstarGtlminus1==1 & 
                      IlLtl_23.FV==0     & 
                      dplyr::lag(InonAdherentToA_v1, default=0)==0, # aka C_{l-2}
                      1, 0),
    
    # Composite at-risk indicator for L2FU model: 
    RsupH  = if_else( ITstarGtlminus1==1 & 
                      InonAdherentToA_v1==0, # aka C_{l-1}
                      1, 0)) %>%
  select( -time ) 
##----------------------------------------------------------------
##  Prep the data for fitting VUPM
##----------------------------------------------------------------
# obtain vector of l among those records 'at risk' in the VUPM
l_VUPM   <- (cohortDataLong %>% filter( ( l>0 ) & (RsupD == 1) ))$l
# compute default knot locations for 4 knots 
lKnots   <- Hmisc::rcspline.eval(l_VUPM, nk=4, knots.only=TRUE)

# append spline basis to dataset
basisMat <- Hmisc::rcspline.eval(cohortDataLong$l, knots=lKnots, nk=4, inclx=TRUE)
colnames(basisMat) <- c("lspl1_VUPM", "lspl2_VUPM", "lspl3_VUPM")
cohortDataLong     <- cohortDataLong %>% bind_cols(basisMat) 

# Check, should be TRUE
print(all.equal(cohortDataLong$l, cohortDataLong$lspl1_VUPM))
##----------------------------------------------------------------
##  Prep the data for fitting LTFU model
##----------------------------------------------------------------
# obtain vector of l among those records 'at risk' in the L2FU model
l_L2FU   <- (cohortDataLong %>% filter( ( l>0 ) & (RsupH == 1) ))$l
# compute default knot locations for 4 knots 
lKnots   <- Hmisc::rcspline.eval(l_L2FU, nk=4, knots.only=TRUE)

# append spline basis to dataset
basisMat <- Hmisc::rcspline.eval(cohortDataLong$l, knots=lKnots, nk=4, inclx=TRUE)
colnames(basisMat) <- c("lspl1_L2FU", "lspl2_L2FU", "lspl3_L2FU")
cohortDataLong     <- cohortDataLong %>% bind_cols(basisMat) 

# Check, should be TRUE
print(all.equal(cohortDataLong$l, cohortDataLong$lspl1_L2FU))
###########################################################################
###
###              Output 
###   
###########################################################################
# save a permanent copy of analytic cohort data in long format
saveRDS(cohortDataLong, file = paste0(dataDir, 'cohortDataLong.rds'))