# Justin DeMonte  
# 240305

# append spline basis matrices to deli dataset
library("tidyverse")
##----------------------------------------------------------------
##  Read in permanent datasets
##----------------------------------------------------------------
dataDir      <- "data/"
# deli dataset
d            <- readRDS(file=paste0(dataDir, "d_t.rds"))

# constants args provided to upstream programs
numTimePts   <- readRDS(file=paste0(dataDir, "numTimePts.rds"))
numTrials    <- readRDS(file=paste0(dataDir, "numTrials.rds"))

pSplineBasis <- readRDS(file=paste0(dataDir, 'pSplineBasis_pooled.rds'))
cSplineBasis <- readRDS(file=paste0(dataDir, 'cSplineBasis_pooled.rds'))  

oSplineBasis_kPlusTrial <- readRDS(file=paste0(dataDir, 'oSplineBasis_kPlusTrial.rds'))
oSplineBasis_k          <- readRDS(file=paste0(dataDir, 'oSplineBasis_k.rds'))  
##----------------------------------------------------------------
##  Add spline basis columns to deli dataset
##----------------------------------------------------------------
for (j in 0:(numTrials-1)){
  
  jindex <- j + 1
  d <- d %>% mutate( # pooled propensity score model
    "prop_j_{j}"        := pSplineBasis[jindex, 1],
    "prop_j'_{j}"       := pSplineBasis[jindex, 2],
    "prop_j''_{j}"      := pSplineBasis[jindex, 3])
  
  for (time in jindex:numTimePts){ 
    k <- time - j # Shift to trial time scale
    d <- d %>% mutate( 
      "cens_kPlusJ_{j}_{k}"   := cSplineBasis[time, 1],
      "cens_kPlusJ'_{j}_{k}"  := cSplineBasis[time, 2],
      "cens_kPlusJ''_{j}_{k}" := cSplineBasis[time, 3],
      "oc_kPlusJ_{j}_{k}"     := oSplineBasis_kPlusTrial[time, 1],
      "oc_kPlusJ'_{j}_{k}"    := oSplineBasis_kPlusTrial[time, 2],
      "oc_kPlusJ''_{j}_{k}"   := oSplineBasis_kPlusTrial[time, 3],
      "oc_Z_kPlusJ_{j}_{k}"   := eval(parse(text=paste0("Z_", j))) * oSplineBasis_kPlusTrial[time, 1],
      "oc_Z_kPlusJ'_{j}_{k}"  := eval(parse(text=paste0("Z_", j))) * oSplineBasis_kPlusTrial[time, 2],
      "oc_Z_kPlusJ''_{j}_{k}" := eval(parse(text=paste0("Z_", j))) * oSplineBasis_kPlusTrial[time, 3])
  }
}
for (k in 1:numTimePts){
  for (j in 0:(numTrials-1)){
    d <- d %>% mutate(
      "oc_k_{k}"         := oSplineBasis_k[k, 1],
      "oc_k'_{k}"        := oSplineBasis_k[k, 2],
      "oc_k''_{k}"       := oSplineBasis_k[k, 3],
      "oc_Z_k_{j}_{k}"   := eval(parse(text=paste0("Z_", j))) * oSplineBasis_k[k, 1],
      "oc_Z_k'_{j}_{k}"  := eval(parse(text=paste0("Z_", j))) * oSplineBasis_k[k, 2],
      "oc_Z_k''_{j}_{k}" := eval(parse(text=paste0("Z_", j))) * oSplineBasis_k[k, 3])
  }
}
##----------------------------------------------------------------
## Save permanent dataset for use downstream
##----------------------------------------------------------------
saveRDS(d, file=paste0(dataDir, 'd.rds'))