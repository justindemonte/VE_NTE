# Justin DeMonte
# 240305
library("tidyverse")

outDir_poly     <- "results/poly/"
outDir_spline   <- "results/spline/"

finalRes_poly   <- list(list(), list(), list(), list(), list(), list())
finalRes_spline <- list()

for (scenario in 1:3){
  
  for (analysis in c("naive", "proposed")){
    for (timeFunc in c("poly", "spline")){
      
      if (timeFunc=="poly"){
        readDir  <- paste0("results/", timeFunc, "/scenario", scenario, "/", analysis, "/")
      } else if (timeFunc=="spline"){
        readDir  <- paste0("results/", timeFunc, "/scenario", scenario, "/")
      }
      writeDir   <- paste0(readDir, "results/")
      
      # Read in constants (passed as args to simulation shell script)
      if (timeFunc=="poly"){
        fileSuffix <- paste0(scenario, '_', analysis)
      } else if (timeFunc=="spline"){
        fileSuffix <- paste0(scenario, '_', timeFunc)
      }
      
      N          <- readRDS(file=paste0(readDir, 'N_sc',          fileSuffix, '.rds'))
      numTrials  <- readRDS(file=paste0(readDir, 'numTrials_sc',  fileSuffix, '.rds'))
      numTimePts <- readRDS(file=paste0(readDir, 'numTimePts_sc', fileSuffix, '.rds'))
      num.sims   <- readRDS(file=paste0(readDir, 'NSIM_sc',       fileSuffix, '.rds'))
      
      seTimes_j  <- readRDS(file=paste0(readDir, 'selectj_sc',    fileSuffix, '.rds'))
      seTimes_k  <- readRDS(file=paste0(readDir, 'selectk_sc',    fileSuffix, '.rds'))
      
      trueVE     <- read.csv(file=paste0(readDir, "trueVE_sc",    fileSuffix, '.csv')) 
      
      ##----------------------------------------------------------------
      ##  Put results into lists and evaluate CI coverage
      ##----------------------------------------------------------------
      VEresultsList     <- vector(mode = "list", length = num.sims)  
      logRRresultsList  <- vector(mode = "list", length = num.sims)  
      SEresultsList     <- vector(mode = "list", length = num.sims)
      CIresultsList     <- vector(mode = "list", length = num.sims)
      hypTestResList    <- vector(mode = "list", length = num.sims)
      
      # if a replication failed, and print
      # the seed value to investigate further
      badSeeds        <- 0 
      
      numVEests       <- length(seTimes_j) + length(seTimes_k) # constant
      CIcover         <- rep(0, numVEests) # initialize CI coverage counter
      if (analysis=="proposed" | timeFunc=="spline"){
        numTEHrejected  <- 0
        numWaneRejected <- 0
      }
      for(sim in 1:num.sims){
        ##----------------------------------------------------------------
        ##  Put results into lists 
        ##----------------------------------------------------------------
        filenm <- paste0(readDir, "VEresults_", sim, '_sc', fileSuffix, ".csv")
        if (file.exists(filenm)){
          VEresultsList[[sim]]    <- read.csv(file = filenm) 
          logRRresultsList[[sim]] <- log(1-VEresultsList[[sim]])
        } else {
          print(paste0(filenm, " does not exist."))
      	  badSeeds <- badSeeds + 1	
        }
        filenm <- paste0(readDir, "SEresults_", sim, '_sc', fileSuffix, ".csv")
        if (file.exists(filenm)){
          SEresultsList[[sim]] <- read.csv(file = filenm) 
        } else {
          print(paste0(filenm, " does not exist."))
        }
        ##----------------------------------------------------------------
        ##  Evaluate CI coverage
        ##----------------------------------------------------------------
        filenm <- paste0(readDir, "CIresults_", sim, '_sc', fileSuffix, ".csv")
        if (file.exists(filenm)){
          CIresults_i <- read.csv(file = filenm)
          for (point in 1:numVEests){
            if ((CIresults_i$lowerCL[point] < trueVE$x[point]) & 
                (CIresults_i$upperCL[point] > trueVE$x[point])) {
                    CIcover[point] <- CIcover[point] + 1 # increment CI coverage counter
            }  
          }
        } else {
          print(paste0(filenm, " does not exist."))
        }
        ##----------------------------------------------------------------
        ##  evaluate hypothesis test results
        ##----------------------------------------------------------------
        if (analysis=="proposed" | timeFunc=="spline"){
          filenm <- paste0(readDir, "hypTestRes_", sim, '_sc', fileSuffix, ".csv")
          if (file.exists(filenm)){
            HTresults_i   <- read.csv(file = filenm)
            TEHtestRes_i  <- HTresults_i$x[1]
            waneTestRes_i <- HTresults_i$x[2]
            if (TEHtestRes_i  < .05){
              numTEHrejected  <- numTEHrejected + 1
            }
            if (waneTestRes_i < .05){
              numWaneRejected <- numWaneRejected + 1
            }
          } else { # file does not exist
            print(paste0(filenm, " does not exist."))
          }
        }
      }
        
      ##----------------------------------------------------------------
      ##  Summarize results
      ##----------------------------------------------------------------
      succSims 	  <- num.sims - badSeeds
      print(paste("There were ", badSeeds, " bad seeds.", sep=""))
      
      # results vectors, each of length numVEests:
      empAveVE    <- apply(do.call(cbind,VEresultsList), 1, mean)
      bias        <- empAveVE-trueVE$x
      
      empSE       <- apply(do.call(cbind, logRRresultsList), 1, sd)
      aveOfSEests <- apply(do.call(cbind, SEresultsList),    1, mean)
      CIcoverRes  <- (CIcover/succSims)*100 # CI coverage as percentage of successful sims
      
      if (analysis=="proposed" | timeFunc=="spline"){
        TEHhypTestRes <- (numTEHrejected/succSims)*100  
        waneHypTstRes <- (numWaneRejected/succSims)*100
        write.csv(TEHhypTestRes, paste0(writeDir, 
          "TEHhypTestRes",  "_", timeFunc, "_sc", fileSuffix, ".csv"))
        write.csv(waneHypTstRes, paste0(writeDir, 
          "waneHypTestRes", "_", timeFunc, "_sc", fileSuffix, ".csv"))
      }
      
      # prepare LaTeX results output
      
      # 1 for naive, 2 for proposed
      anIndex <- as.integer(analysis=="naive") + 
        (2*as.integer(analysis=="proposed"))
      
      if (timeFunc=="poly"){
        finalRes_poly[[scenario]][[anIndex]] <- 
          data.frame((trueVE$x*100), (bias*100), (empSE*100), (aveOfSEests*100), CIcoverRes)  
      } else if (timeFunc=="spline"){
        finalRes_spline[[scenario]] <- 
          data.frame((trueVE$x*100), (bias*100), (empSE*100), (aveOfSEests*100), CIcoverRes)  
      }
    }
  }
}
##----------------------------------------------------------------
## final results in format ready for LaTeX table 
##----------------------------------------------------------------
# placement of horizontal lines in LaTeX table
hlines <- seq(5, 25, by=5)

# add informational columns
scenCol <- c(1, rep(NA, 9),
             2, rep(NA, 9), 
             3, rep(NA, 9))  

estimandCol <- c("$VE_0(5)$",
                 "$VE_3(5)$",
                 "$VE_6(5)$",
                 "$VE_{9}(5)$",
                 "$VE_{12}(5)$",
                 "$VE_5(1)$",
                 "$VE_5(4)$",
                 "$VE_5(8)$",
                 "$VE_5(12)$",
                 "$VE_5(15)$")

##----------------------------------------------------------------
## poly sim results
##----------------------------------------------------------------
poly1  <- bind_cols(estimandCol, 
            finalRes_poly[[1]][[1]], finalRes_poly[[1]][[2]][,-1])
poly2  <- bind_cols(estimandCol, 
            finalRes_poly[[2]][[1]], finalRes_poly[[2]][[2]][,-1])
poly3  <- bind_cols(estimandCol, 
            finalRes_poly[[3]][[1]], finalRes_poly[[3]][[2]][,-1])
resFileName_poly <- "simRes_poly.txt"


final_poly_t    <- bind_rows(poly1, poly2, poly3)
final_poly      <- bind_cols(scenCol, final_poly_t)

finalRes_poly_x <- xtable::xtable(final_poly, digits=c(0,0,0,1,1,1,1,0,1,1,1,0))
print(finalRes_poly_x, 
      include.rownames=FALSE, 
      hline.after=hlines,
      file=paste0(outDir_poly, resFileName_poly), 
      sanitize.text.function=identity)

##----------------------------------------------------------------
## spline sim results
##----------------------------------------------------------------
scn1Final_spline    <- bind_cols(estimandCol, finalRes_spline[[1]])
scn2Final_spline    <- bind_cols(estimandCol, finalRes_spline[[2]])
scn3Final_spline    <- bind_cols(estimandCol, finalRes_spline[[3]])
resFileName_spline  <- "simRes_spline.txt"

final_spline_t <- bind_rows(scn1Final_spline, scn2Final_spline, scn3Final_spline)
final_spline   <- bind_cols(scenCol, final_spline_t)

finalRes_spline_x <- xtable::xtable(final_spline, digits=c(0,0,0,1,1,1,1,0))
print(finalRes_spline_x,
      include.rownames=FALSE,
      hline.after=hlines,
      file=paste0(outDir_spline, resFileName_spline),
      sanitize.text.function=identity)