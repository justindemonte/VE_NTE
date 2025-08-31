##----------------------------------------------------------------
##  Helper functions
##----------------------------------------------------------------
# Helper functions in this program employ a binary treatment variable Z.
# 'Z' in these programs corresponds to A (treatment assignment) in the manuscript.  

# discrete hazard at a single time point as a function of regression params
# takes parameter vector of length 8, returns a scalar (hazard).
getHaz       <- function(Z, vaccTime, k, betaVec) {
  plogis( betaVec[1] + 
            ( betaVec[2] * (vaccTime+k) ) + 
            ( betaVec[3] * ((vaccTime+k)^2) ) +
            ( betaVec[4] * Z ) + 
            ( betaVec[5] * (vaccTime+k) * Z ) + 
            ( betaVec[6] * ((vaccTime+k)^2) * Z) +
            ( betaVec[7] * k * Z ) + 
            ( betaVec[8] * (k^2) * Z) )
}
# takes parameter vector of length 6, returns a scalar (hazard).
# Use this version to compute estimated haz from naive model parameter estimates
getHaz_naive <- function(Z, vaccTime, k, betaVec) { 
  plogis( betaVec[1] +     
            ( betaVec[2] * (vaccTime+k) ) +    
            ( betaVec[3] * ((vaccTime+k)^2) ) + 
            ( betaVec[4] * Z ) + 
            ( betaVec[5] * k * Z ) + 
            ( betaVec[6] * (k^2) * Z) )
}
# obtain trial-specific VE over k.  Takes a vector, returns a vector.
getVE        <- function(haz1, haz0){
  1 - ( (1-cumprod(1-haz1)) / (1-cumprod(1-haz0)) )
}
# Obtain all hazards for A \in {0, 1} 
# col j+1 of hazMatz holds lambda_j^z(k) for k=1 to tau-j
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
# Obtain matrix of VE across all j and k 
# hazMatList must have hazard matrix for z=1 as first element, 
# hazard matrix for z=0 as second element.
getVEmatrix   <- function(hazMatList){
  VEmat   <- matrix(nrow=tau, ncol=numTrials)
  for (j in 0:J){
    VEmat[1:(tau-j), (j+1)] <- 
      getVE(hazMatList[[1]][1:(tau-j), j+1], hazMatList[[2]][1:(tau-j), j+1])
  }
  return(VEmat)
} 