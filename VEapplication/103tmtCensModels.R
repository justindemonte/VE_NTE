# 103tmtCensModels.R
#
# input: 
# 1. analytic cohort data (as described in Section 2) in long format 
# /nas/longleaf/home/demjus01/VEapplication/data/cohortDataLong.rds 
#
# 2. numTrials, a constant from upstream program
# /nas/longleaf/home/demjus01/VEapplication/data/numTrials.rds 
#
# output:
#
# 1. model fit object for VUPM 
# /nas/longleaf/home/demjus01/VEapplication/data/dFit.rds 
#
# 2. model fit object for L2FU model
# /nas/longleaf/home/demjus01/VEapplication/data/hFit.rds 
library("tidyverse")
dataDir    <- "data/"

# read in user-provided constants from upstream program
numTrials  <- readRDS(file=paste0(dataDir, "numTrials.rds"))
J          <- numTrials - 1

# end boilerplate
cohortDataLong <- readRDS(file = paste0(dataDir, "cohortDataLong.rds")) 
##----------------------------------------------------------------
##  Fit the VUPM model
##----------------------------------------------------------------
dFit <- glm( D ~ lspl1_VUPM + lspl2_VUPM + lspl3_VUPM +
                         X1_c + X2 + X3      + 
                         IlagZeq1 + IlagZeq2 +
                         IpreGrace_Pfizer    + 
                         IpreGrace_Moderna   +
                         IpreGrace_AstraZ    +
                         Igrace_mRNA         +
                         Igrace_AstraZ,
             family = binomial( link = "logit" ),
             data   = cohortDataLong,
             subset = ( l>0 ) & (RsupD == 1) ) # reffered to as R^Z in supplement
print(dFit)
##----------------------------------------------------------------
##  Fit the LTFU model 
##----------------------------------------------------------------
hFit <- glm( H ~ lspl1_L2FU + lspl2_L2FU + lspl3_L2FU +
             Zsup1 + Zsup2 + 
             X1_c + X2 + X3 + X1X3,
             family = binomial( link = "logit" ),
             data   = cohortDataLong,
             subset = ( l>0 ) & (RsupH == 1) )
print(hFit)
##----------------------------------------------------------------
##  Save a permanent copy of model fit objects
##----------------------------------------------------------------
saveRDS(dFit,  file = paste0(dataDir, 'dFit.rds'))
saveRDS(hFit,  file = paste0(dataDir, 'hFit.rds'))