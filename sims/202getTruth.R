library("dplyr")
library("tibble")
library("tidyr")
library("purrr")
###########################################################################
###
###   Generate superpopulation and calculate true VE empirically
###
###########################################################################
## collect args passed in from SLURM job script
args          <- commandArgs(trailingOnly=TRUE)
N             <- as.numeric(args[1])
numTrials     <- as.numeric(args[2])
numTimePts    <- as.numeric(args[3])
ageProdTerm   <- as.numeric(args[4])
simulation    <- as.character(args[5])
scenario      <- Sys.getenv("SLURM_ARRAY_TASK_ID")
# end collect args passed in from SLURM job script

tau           <- numTimePts - 1
J             <- numTrials  - 1
# generate superpopulation 
set.seed(3001)
source("201DGP.R")

k_0           <- 5
j_0           <- 5
seTimes_j     <- c( 0, 3, 6,  9, 12 )
seTimes_k     <- c( 1, 4, 8, 12, 15 )
###########################################################################
###
###  Calculate VE_j(k) and VE^s_j(k) in the superpopulation
###
###########################################################################
## VE_j(k) = 1 - { E(Y^1_{j+k} | E_j=1) / E(Y^0_{j+k} | E_j=1) }
if (simulation=="main"){
  trueVEjk     <- vector( mode = "numeric", length = 10 )
  trueVEjkSups <- vector( mode = "numeric", length = 10 )  
} else {
  trueVEjk     <- vector( mode = "numeric", length = 5 )
  trueVEjkSups <- vector( mode = "numeric", length = 5 ) 
}
counter <- 1
k       <- 5
for ( j in seTimes_j ){  
  
  ###########################################################################
  ##  compute VE_j(5) for select j
  ###########################################################################
  # sum_i ( Y^{v_1=j+1}_{j+k} E_j=1 ) / sum_i ( E_j )
  risk1_jk <- sum( ( simd[ paste0( "TsupV1eq", ( j+1 ) )] <= ( j+k ) )
                   * simd[ paste0( "E_", j )] ) 
  
  # sum_i ( Y^{v_1>tau}_{j+k} E_j=1 ) / sum_i ( E_j )
  risk0_jk <- sum( 
    ( simd["TsupV1gtTau"] <= ( j+k ) ) * simd[paste0( "E_", j )] ) 
  
  # denominators in above comments will cancel in ratio
  trueVEjk[counter] <- 1 - ( risk1_jk / risk0_jk )
  ###########################################################################
  ##  compute VE^s_j(5) for select j
  ###########################################################################
  unvaccHazValues <- simd %>% select(starts_with("oneMinuslambda0"))
  vaccHazValues   <- simd %>% select(starts_with(paste0("oneMinuslambdaSup", (j+1), "_")))
  
  # (1/n) \sum_{i=1}^n E[Y^0_j(k)\mid E_j=1, X_i]
  unvaccMat    <- as.matrix(unvaccHazValues)
  risk0_jkSups <- sum(
    1-matrixStats::rowCumprods(unvaccMat[,(j+1):ncol(unvaccMat)])[,k]
  )
  # (1/n) \sum_{i=1}^n E[Y^1_j(k)\mid E_j=1, X_i]
  risk1_jkSups <- sum(
    1-matrixStats::rowCumprods(as.matrix(vaccHazValues))[,k]
  )
  # (1/n) in above comments will cancel in ratio
  trueVEjkSups[counter]  <- 1 - ( risk1_jkSups / risk0_jkSups )
  
  counter                <- counter + 1
}
if (simulation=="main"){
  j <- 5
  for ( k in seTimes_k ){  # compute VE_5(k) for select k
    
    ## In the following, denominators in comments will cancel in ratio
    # sum_i ( Y^{v_1=j+1}_{j+k} E_j=1 ) / sum_i ( E_j )
    risk1_jk <- sum( ( simd[ paste0( "TsupV1eq", ( j+1 ) )] <= ( j+k ) )
                     * simd[ paste0( "E_", j )] )
    
    # sum_i ( Y^{v_1>tau}_{j+k} E_j=1 ) / sum_i ( E_j )
    risk0_jk <- sum( ( simd["TsupV1gtTau"] <= ( j+k ) ) * simd[paste0( "E_", j )] )
    
    trueVEjk[counter]      <- 1 - ( risk1_jk / risk0_jk )
    counter                <- counter + 1
  }  
}
if ( simulation=="main" ){
  directory <- paste0( "results/true/scenario", scenario, "/unstandardized/" )
  saveRDS( trueVEjk,      file=paste0( directory, 'trueValues_m', scenario, '.rds' ))
} else {
  directory <- paste0( "results/true/scenario", scenario, "/standardized/" )
  saveRDS( trueVEjkSups, file=paste0( directory, 'trueValues_e', scenario, '.rds' ))
}