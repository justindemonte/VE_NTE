# 104prepSTEdata.R
#
# input: 
# 1. analytic cohort data (as described in Section 2) in long format 
# /nas/longleaf/home/demjus01/VEapplication/data/cohortDataLong.rds 
#
# 2. numTrials, a constant, number of trials to be emulated(i.e., J-1)
# /nas/longleaf/home/demjus01/VEapplication/data/numTrials.rds
# 
# 3. numTimePoints - a constant 
# /nas/longleaf/home/demjus01/VEapplication/data/numTimePts.rds
#
# 4. model fit object for VUPM 
# /nas/longleaf/home/demjus01/VEapplication/data/dFit.rds 
#
# 5. model fit object for L2FU model
# /nas/longleaf/home/demjus01/VEapplication/data/hFit.rds 
#
# output:
# 1. STE analysis dataset
# /nas/longleaf/home/demjus01/VEapplication/data/STEdataLong.rds 
#
# 2. base deli dataset
# /nas/longleaf/home/demjus01/VEapplication/data/deliData.rds 
library("tidyverse")

dataDir    <- "data/"
# read in user-provided constants from upstream program
numTimePts <- readRDS(file=paste0(dataDir, "numTimePts.rds"))
numTrials  <- readRDS(file=paste0(dataDir, "numTrials.rds"))

tau        <- numTimePts - 1
J          <- numTrials  - 1

# read in model fit objects from upstream program
dFit       <- readRDS(file = paste0(dataDir, 'dFit.rds'))
hFit       <- readRDS(file = paste0(dataDir, 'hFit.rds'))

cohortDataLong <- readRDS(file=paste0(dataDir, "cohortDataLong.rds"))
# end boilerplate
##----------------------------------------------------------------
## Obtain predicted values from fitted nuisance models
##----------------------------------------------------------------
# Compute predicted probs of LTFU 
cohortDataLong$predLambdaH   <- predict( hFit, newdata = cohortDataLong, type="response" )

# Compute predicted VU process transition probs 
cohortDataLong$predTransProb <- predict( dFit, newdata = cohortDataLong, type="response" )
##----------------------------------------------------------------
##  Calculate estimated IPT weight for each time point
##----------------------------------------------------------------
cohortDataLong <- cohortDataLong %>%
  
  mutate(lagLambdaC = case_when(
    RsupD == 0                                     ~ 0, # referred to as R^Z in supplement
    lagZ  == 0 & Z==1                              ~ 1-predTransProb, # this is 1-p^{0,1}   
    lagZ  == 0                                     ~ predTransProb, 
    (lagZ == 1) & (IpreGrace_Pfizer==1 | IpreGrace_Moderna==1 | IpreGrace_AstraZ==1) ~ predTransProb, 
    (lagZ == 1) & (IstateOfGrace==1)               ~ 0, 
    (lagZ == 1) & (IfallFromGrace==1)              ~ 1-predTransProb, 
    IlLtl_23.FV == 1                               ~ 0, 
    (IlLtl_23.FV==0) & (lagZ == 1) & (vaccinetypedose1==4)            ~ predTransProb, 
    (IlLtl_23.FV==0) & (lagZ == 2) & (vaccinetypedose1 %in% c(1,2,3)) ~ predTransProb, 
    TRUE                                           ~ NA_real_)) %>%
  mutate( IPTwt  = 1/(1-lagLambdaC), 
          IPCwt  = 1/(1-predLambdaH)) %>%
  mutate( IPTCwt = IPTwt * IPCwt)     
###########################################################################
###
###   Prepare wide dataset with 1 obs per person for delicatessen
###
###########################################################################
deliData <- cohortDataLong    %>%
  select( "id", "H", "Y", "Z", "D", "l", "X1_c", "X2", "X3", "X1X3", 
          "lspl1_L2FU", "lspl2_L2FU", "lspl3_L2FU",
          "lspl1_VUPM", "lspl2_VUPM", "lspl3_VUPM",
          "IstateOfGrace", "IfallFromGrace",
          "IpreGrace_Pfizer",
          "IpreGrace_Moderna",
          "IpreGrace_AstraZ",
          "Igrace_mRNA",
          "Igrace_AstraZ",
          "IlLtl_23.FV",
          "Zsup1", "Zsup2", 
          "IlagZeq0", "IlagZeq1", "IlagZeq2",             
          "InonAdherentToA_v1",               
          "RsupD", "RsupH", 
          starts_with("E_") ) %>%
  pivot_wider(
    names_from  = l,
    values_from = c( H, Y, Z, D, lspl1_L2FU, lspl2_L2FU, lspl3_L2FU,
                     lspl1_VUPM, lspl2_VUPM, lspl3_VUPM,
                     IstateOfGrace, IfallFromGrace, IlLtl_23.FV,
                     IpreGrace_Pfizer,
                     IpreGrace_Moderna,
                     IpreGrace_AstraZ,
                     Igrace_mRNA,
                     Igrace_AstraZ,
                     Zsup1, Zsup2, 
                     IlagZeq0, IlagZeq1, IlagZeq2, 
                     InonAdherentToA_v1,
                     RsupD, RsupH),
    names_sep   = "_") 
###########################################################################
###
###   Create STE analysis dataset
###
###########################################################################
# for each trial j in 0 to J, create a flag indicating that
# a given row belongs in the trial-j dataset (BITJDS)
# defined as follows:
# if an individual is eligible for trial j,
# then rows with j <= l <= tau belong in the trial-j dataset (BITJDS)
for ( j in 0:J ){
  cohortDataLong[ paste0( "BITJDS_", j) ] <-
    cohortDataLong[ paste0( "E_", j) ] * ( cohortDataLong[ "l" ] >= j )
}
STEdataLong <- cohortDataLong %>%
  select( c( "id", "l", "H", "Y", "Z", "lagH", "lagY", "lagZ",
             "X1_c", "X2", "X3", "X1X3", 
             "IPTCwt", "InonAdherentToA_v1",
             num_range( "BITJDS_", 0:J )
  )) %>%
  pivot_longer(
    cols      = c( starts_with( "BITJDS_" ) ),
    names_to  = c( ".value", "trialNo" ),
    names_sep = "\\_" )  %>%
  mutate( j   = as.numeric( trialNo ))         %>%
  filter( BITJDS == 1 )                        %>%
  select( -c( "trialNo", "BITJDS" ) )          %>%
  arrange( id, j, l )                          %>%
  mutate( k = l-j )                            %>% # derive trial time
  # ensure that cumulative product of IPTC weights below is taken starting at k=1
  mutate( IPTCwt = if_else( k==0, 1, IPTCwt) ) %>% 
  group_by( id, j )                            %>%
    # trial-j treatment assignment based on vacc uptake during 1st week of trial j
    mutate( Z_keq1 = Z[ k==1 ] )               %>%
    mutate( A     = if_else(Z_keq1 > 0, 1, 0)) %>%
    # 'compute trial-specific censoring indicator
    mutate( C = case_when(
      ( (A == 0) & (Z > 0) )                   ~ 1,
      ( (A == 1) & (InonAdherentToA_v1 == 1) ) ~ 1,
      TRUE ~ 0))                               %>%
    mutate( R = case_when (
      ( C == 1 | lagY == 1 | H == 1 )       ~ 0,
      TRUE ~ 1 ))                              %>%
    mutate(Wjk = cumprod(IPTCwt))              %>%
  ungroup()                                    %>% 
  select(-c(Z_keq1))                           
##----------------------------------------------------------------
##  Construct spline basis matrices for outcome model
##----------------------------------------------------------------
# obtain vector of l among those records at risk in the outcome model
l_Ymodel <- (STEdataLong %>% filter(( R == 1 ) & ( k != 0 )))$l
# compute default knot locations for 4 knots 
lKnots   <- Hmisc::rcspline.eval(l_Ymodel, nk=4, knots.only=TRUE)

# append spline basis to dataset
l_Ybasis <- Hmisc::rcspline.eval(STEdataLong$l, knots=lKnots, nk=4, inclx=TRUE)
colnames(l_Ybasis) <- c("lspl1_Y", "lspl2_Y", "lspl3_Y")
STEdataLong <- STEdataLong %>% bind_cols(l_Ybasis) 

# Check, should be TRUE
print(all.equal(STEdataLong$l, STEdataLong$lspl1_Y))

# obtain vector of k among those records at risk in the outcome model
k_Ymodel <- (STEdataLong %>% filter(( R == 1 ) & ( k != 0 )))$k
# compute default knot locations for 4 knots 
kKnots   <- Hmisc::rcspline.eval(k_Ymodel, nk=4, knots.only=TRUE)

# append spline basis to dataset
k_Ybasis <- Hmisc::rcspline.eval(STEdataLong$k, knots=kKnots, nk=4, inclx=TRUE)
colnames(k_Ybasis) <- c("kspl1_Y", "kspl2_Y", "kspl3_Y")
STEdataLong <- STEdataLong %>% bind_cols(k_Ybasis) 

# Check, should be TRUE
print(all.equal(STEdataLong$k, STEdataLong$kspl1_Y))

# obtain vector of X1_c among those records at risk in the outcome model
X1_Ymodel <- (STEdataLong %>% filter(( R == 1 ) & ( k != 0 )))$X1_c
# compute default knot locations for 4 knots 
X1knots   <- Hmisc::rcspline.eval(X1_Ymodel, nk=4, knots.only=TRUE)

# append spline basis to dataset
X1_Ybasis <- Hmisc::rcspline.eval(STEdataLong$X1_c, knots=X1knots, nk=4, inclx=TRUE)
colnames(X1_Ybasis) <- c("X1_cspl1", "X1_cspl2", "X1_cspl3")
STEdataLong <- STEdataLong %>% bind_cols(X1_Ybasis) 

# Check, should be TRUE
print(all.equal(STEdataLong$X1_c, STEdataLong$X1_cspl1))
##----------------------------------------------------------------
##  Prepare rcs basis matrices for 'o'utcome model
##----------------------------------------------------------------
l_Ybasis  <- as.data.frame(l_Ybasis)  %>% distinct() %>% arrange(lspl1_Y)
k_Ybasis  <- as.data.frame(k_Ybasis)  %>% distinct() %>% arrange(kspl1_Y)
X1_Ybasis <- as.data.frame(X1_Ybasis) %>% distinct()
# Check, should be TRUE
print( nrow(X1_Ybasis)==length(unique(deliData$X1_c)) )

for ( m in 1:tau ){
  deliData <- deliData %>% mutate( # add columns to the deli ds for spline basis matrices
    "oSpline_l_{m}"   := l_Ybasis[[(m+1), 1]],
    "oSpline_l'_{m}"  := l_Ybasis[[(m+1), 2]],
    "oSpline_l''_{m}" := l_Ybasis[[(m+1), 3]],
    "oSpline_k_{m}"   := k_Ybasis[[(m+1), 1]],
    "oSpline_k'_{m}"  := k_Ybasis[[(m+1), 2]],
    "oSpline_k''_{m}" := k_Ybasis[[(m+1), 3]])
}
deliData <- left_join(deliData, X1_Ybasis, by=join_by(X1_c==X1_cspl1), relationship="many-to-one", keep=TRUE)

#check: deliData cols X1_cspl1, X1_cspl2, X1_cspl3 should have no NAs
sum(is.na(deliData$X1_cspl1))
sum(is.na(deliData$X1_cspl2))
sum(is.na(deliData$X1_cspl3))
##----------------------------------------------------------------
##  Save permanent copies of analysis data sets
##----------------------------------------------------------------
# save a permanent copy of STE analysis data in long format
saveRDS(STEdataLong, file = paste0(dataDir, 'STEdataLong.rds'))
# save dataset for delicatessen
saveRDS(deliData,    file = paste0(dataDir, 'deliData.rds'))
# save spline basis matrices 
saveRDS(l_Ybasis,    file = paste0(dataDir, 'l_Ybasis.rds'))
saveRDS(k_Ybasis,    file = paste0(dataDir, 'k_Ybasis.rds'))
saveRDS(X1_Ybasis,   file = paste0(dataDir, 'X1_Ybasis.rds'))
