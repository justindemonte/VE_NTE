library("tidyverse")

# output directories
outDir_poly     <- "results/poly/"
outDir_spline   <- "results/spline/"
outDir_ext      <- "results/extension/"

# prepare data structures for results tables. 
finalRes_poly   <- list( list(), list(), list())
finalRes_spline <- list()
finalRes_ext    <- list( list(), list(), list())

for ( s in 1:3 ){
  for ( simulation in c("main", "extension") ){
    for ( analysis in c("naive", "proposed") ){
      for ( timeFunc in c("poly", "spline") ){
        if ( timeFunc=="spline" & (simulation=="extension" | analysis=="naive") ){break}
        
        # used in file names
        simCode       <- ifelse( simulation=="main", "m", "e" )
        anlCode       <- ifelse( analysis=="naive",  "n", "p" )
        splCode       <- ifelse( timeFunc=="poly",   "p", "s" )
        settingCode   <- paste0( simCode, anlCode, splCode, s )
        
        if ( simulation=="main" ){
          readDirTrue <- paste0( "results/true/scenario", s, "/unstandardized/" )
          readDirEst  <- paste0( "results/", timeFunc, "/scenario", s, "/", analysis, "/" )
          writeDir    <- paste0( "results/", timeFunc, "/scenario", s, "/" )
          trueVE      <- readRDS(file=paste0(readDirTrue,"trueValues_m", s, '.rds')) 
        } else { # simulation = "extension"
          readDirTrue <- paste0( "results/true/scenario",      s, "/standardized/" )
          readDirEst  <- paste0( "results/extension/scenario", s, "/", analysis, "/" )
          writeDir    <- paste0( "results/extension/scenario", s, "/" )
          trueVE      <- readRDS(file=paste0(readDirTrue,"trueValues_e", s, '.rds')) 
        }
        # Read in constants (passed as args to simulation shell script)
        num.sims   <- readRDS(file=paste0(readDirEst, 'NSIM_',    settingCode, '.rds'))
        seTimes_j  <- readRDS(file=paste0(readDirEst, 'selectj_', settingCode, '.rds'))
        seTimes_k  <- readRDS(file=paste0(readDirEst, 'selectk_', settingCode, '.rds'))
        ##----------------------------------------------------------------
        ##  Put results into lists and evaluate CI coverage
        ##----------------------------------------------------------------
        VEresList      <- vector(mode = "list", length = num.sims)  
        logRRresList   <- vector(mode = "list", length = num.sims)  
        SEresList      <- vector(mode = "list", length = num.sims)
        CIresList      <- vector(mode = "list", length = num.sims)
        hypTestResList <- vector(mode = "list", length = num.sims)
        
        # if a replication failed, and print
        # the seed value to investigate further
        badSeeds         <- 0 
        
        if (simulation=="main"){
          numVEests      <- length(seTimes_j) + length(seTimes_k) # constant  
        } else{
          numVEests      <- length(seTimes_j) # constant  
        }
        CIcover          <- rep(0, numVEests) # initialize CI coverage counter
        if ( analysis=="proposed" ){ numTEHrejected <- 0 }
        for( sim in 1:num.sims ){
          ##----------------------------------------------------------------
          ##  Put results into lists 
          ##----------------------------------------------------------------
          filenm <- paste0(readDirEst, "VEres_", settingCode, '_', sim, ".csv")
          if (file.exists(filenm)){
            VEresList[[sim]]    <- read.csv(file = filenm) 
            logRRresList[[sim]] <- log(1-VEresList[[sim]])
          } else {
            print(paste0(filenm, " does not exist."))
        	  badSeeds <- badSeeds + 1	
          }
          filenm <- paste0(readDirEst, "SEres_", settingCode, '_', sim, ".csv")
          if (file.exists(filenm)){
            SEresList[[sim]]    <- read.csv(file = filenm) 
          } else {print(paste0(filenm, " does not exist."))}
          ##----------------------------------------------------------------
          ##  Evaluate CI coverage
          ##----------------------------------------------------------------
          filenm <- paste0(readDirEst, "CIres_", settingCode, '_', sim, ".csv")
          if ( file.exists(filenm) ){
            CIres_i <- read.csv(file = filenm)
            for ( point in 1:numVEests ){
              if ( (CIres_i$lowerCL[point] < trueVE[point]) & 
                  (CIres_i$upperCL[point] > trueVE[point]) ){
                      CIcover[point] <- CIcover[point] + 1 # increment CI coverage counter
              }  
            }
          } else {print(paste0(filenm, " does not exist."))}
          ##----------------------------------------------------------------
          ##  evaluate hypothesis test results
          ##----------------------------------------------------------------
          if (analysis=="proposed" & simulation=="main"){
            filenm <- paste0(readDirEst, "hypTestRes_", settingCode, '_', sim, ".csv")
            if (file.exists(filenm)){
              HTres_i      <- read.csv(file = filenm)
              TEHtestRes_i <- HTres_i$x[1]
              if (TEHtestRes_i  < .05){
                numTEHrejected  <- numTEHrejected + 1
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
        empAveVE    <- apply(do.call(cbind,VEresList), 1, mean)
        bias        <- empAveVE-trueVE
        
        empSE       <- apply(do.call(cbind, logRRresList), 1, sd)
        aveOfSEests <- apply(do.call(cbind, SEresList),    1, mean)
        CIcoverRes  <- (CIcover/succSims)*100 # CI coverage as percentage of successful sims
        
        if (analysis=="proposed" & simulation=="main"){
          TEHhypTestRes <- (numTEHrejected/succSims)*100  
          write.csv(TEHhypTestRes, paste0(writeDir, 
            "TEHhypTestRes",  "_", settingCode, ".csv"))
        }
        ##----------------------------------------------------------------
        ##  prepare LaTeX results output
        ##----------------------------------------------------------------
        # 1 for naive, 2 for proposed
        anIndex <- as.integer(analysis=="naive") + (2*as.integer(analysis=="proposed"))
        
        if (simulation=="main"){
          if (timeFunc=="poly"){
            finalRes_poly[[s]][[anIndex]] <- 
              data.frame((trueVE*100), (bias*100), (empSE*100), (aveOfSEests*100), CIcoverRes)  
          } else if (timeFunc=="spline"){
            finalRes_spline[[s]] <- 
              data.frame((trueVE*100), (bias*100), (empSE*100), (aveOfSEests*100), CIcoverRes)  
          }
        } else { # simulation="extension"
          finalRes_ext[[s]][[anIndex]] <- 
            data.frame((trueVE*100),   (bias*100), (empSE*100), (aveOfSEests*100), CIcoverRes)    
        }
      }
    }
  }
}
###----------------------------------------------------------------
### 
### final main sim results in format ready for LaTeX table 
###
###----------------------------------------------------------------
# placement of horizontal lines in LaTeX table
hlines <- c(-1, 0, seq(5, 25, by=5))

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

final_poly_t     <- bind_rows(poly1, poly2, poly3)
final_poly       <- bind_cols(scenCol, final_poly_t)

addtorow         <- list()
addtorow$pos     <- list(0,0,0,0,30, 30, 30, 30, 30)
addtorow$command <- c("& & & \\multicolumn{4}{c}{Model (\\ref{main-timeVaryingIntercept})} &\\multicolumn{4}{c}{Model (\\ref{main-calTimeTerms})} \\\\\n ",
                      "&  & True & & &  & 95\\% CI &  &  &  & 95\\% CI\\\\\n",
                      "&  & value &  & ESE$^\\ddagger$ & ASE$^\\ddagger$ & coverage &  & ESE$^\\ddagger$ & ASE$^\\ddagger$ & coverage \\\\\n",
                      "Scn. & Estimand &(\\%)  & Bias & $\\times10^2$ & $\\times10^2$ & (\\%) & Bias & $\\times10^2$ & $\\times10^2$& (\\%) \\\\\n",
                      "\\hline",
                      "\\multicolumn{11}{l}{Abbreviations: ASE, average estimated standard error; ESE, empirical standard error; } \\\\\n",
                      "\\multicolumn{11}{l}{CI, confidence interval}\\\\\n",
                      "\\multicolumn{11}{l}{$^\\ddagger$Empirical standard error and average estimated standard error were computed on the}\\\\\n",
                      "\\multicolumn{11}{l}{log risk ratio scale.  All other quantities were computed on the VE scale.}"
                      )

finalRes_poly_x  <- xtable::xtable(final_poly, 
      digits=c(0,0,0,1,1,1,1,0,1,1,1,0), 
      caption = "Simulation study results by scenario across $3{,}000$ replications, $\\tau=20$ time points of follow-up, and $n=50{,}000$.",
      align = "cccccccccccc"
      )
print(finalRes_poly_x, 
      include.rownames=FALSE, 
      hline.after=hlines,
      file=paste0(outDir_poly, resFileName_poly), 
      sanitize.text.function=identity, 
      caption.placement = "top", 
      add.to.row = addtorow, 
      include.colnames = FALSE
      )
##----------------------------------------------------------------
## spline sim results
##----------------------------------------------------------------
scn1Final_spline   <- bind_cols(estimandCol, finalRes_spline[[1]])
scn2Final_spline   <- bind_cols(estimandCol, finalRes_spline[[2]])
scn3Final_spline   <- bind_cols(estimandCol, finalRes_spline[[3]])
resFileName_spline <- "simRes_spline.txt"

addtorow           <- list()
addtorow$pos       <- list(0,0,0,30, 30, 30, 30, 30, 30)
addtorow$command   <- c("&&True &  & &&  95\\% CI    \\\\\n ",
                        "&&value  & &  ESE$^\\ddagger$ &ASE$^\\ddagger$ &  coverage  \\\\\n",
                        "Scenario &Estimand& (\\%)&Bias& $\\times 10^2$  & $\\times 10^2$ & (\\%) \\\\\n",
                        "\\hline",
                        "\\multicolumn{7}{l}{Abbreviations: ASE, average estimated standard error; ESE,} \\\\\n",
                        "\\multicolumn{7}{l}{empirical standard error; CI, confidence interval} \\\\\n",
                        "\\multicolumn{7}{l}{$^\\ddagger$Empirical standard error and average estimated standard error}\\\\\n",
                        "\\multicolumn{7}{l}{were computed on the log risk ratio scale.  All other quantities}\\\\\n",
                        "\\multicolumn{7}{l}{were computed on the VE scale.}"
                        )
final_spline_t <- bind_rows(scn1Final_spline, scn2Final_spline, scn3Final_spline)
final_spline   <- bind_cols(scenCol, final_spline_t)

finalRes_spline_x <- xtable::xtable(final_spline, 
      digits=c(0,0,0,1,1,1,1,0), 
      caption = "Additional simulation study results by scenario across $3{,}000$ replications, $\\tau=20$ time points of follow-up, and $n=50{,}000$.  All time functions in the analytical models were specified using restricted cubic splines.",
      align = "cccccccc"
      )
print(finalRes_spline_x,
      include.rownames=FALSE,
      include.colnames = FALSE, 
      hline.after=hlines,
      file=paste0(outDir_spline, resFileName_spline),
      sanitize.text.function=identity,
      caption.placement = "top", 
      add.to.row = addtorow
      )
###----------------------------------------------------------------
### 
### final standardization sim results in format ready for LaTeX table 
###
###----------------------------------------------------------------
# placement of horizontal lines in LaTeX table
hlines  <- c(-1, 0, seq(5, 10, by=5))

# add informational columns
scenCol <- c(1, rep(NA, 4),
             2, rep(NA, 4), 
             3, rep(NA, 4))  

estimandCol <- c("$VE^s_0(5)$",
                 "$VE^s_3(5)$",
                 "$VE^s_6(5)$",
                 "$VE^s_{9}(5)$",
                 "$VE^s_{12}(5)$")
##----------------------------------------------------------------
## standardization sim results
##----------------------------------------------------------------
res1  <- bind_cols(estimandCol, finalRes_ext[[1]][[1]], finalRes_ext[[1]][[2]][,-1])
res2  <- bind_cols(estimandCol, finalRes_ext[[2]][[1]], finalRes_ext[[2]][[2]][,-1])
res3  <- bind_cols(estimandCol, finalRes_ext[[3]][[1]], finalRes_ext[[3]][[2]][,-1])
resFileName_ext <- "simRes_ext.txt"

final_t <- bind_rows(res1, res2, res3)
final   <- bind_cols(scenCol, final_t)

addtorow           <- list()
addtorow$pos       <- list(0,0,0,0,15, 15, 15, 15, 15)
addtorow$command   <- c("& & & \\multicolumn{4}{c}{Unstandardized} & \\multicolumn{4}{c}{Standardized}  \\\\\n ",
                        "\\cmidrule(lr){4-7}  \\cmidrule(lr){8-11} &  & True & & & & 95\\% CI && & & 95\\% CI \\\\\n",
                        "&  & Value & & ESE$^\\ddagger$ & ASE$^\\ddagger$ & coverage  &  &   ESE$^\\ddagger$ & ASE$^\\ddagger$   & coverage  \\\\\n",
                        "Scn. & Estimand & (\\%) & Bias & $\\times 10^2$ & $\\times 10^2$ & (\\%) & Bias& $\\times 10^2$ & $\\times 10^2$ & (\\%)  \\\\\n",
                        "\\hline",
                        "\\multicolumn{11}{l}{Abbreviations: ASE, average estimated standard error; ESE, empirical standard error;} \\\\\n",
                        "\\multicolumn{11}{l}{CI, confidence interval} \\\\\n",
                        "\\multicolumn{11}{l}{$^\\ddagger$Empirical standard error and average estimated standard error were computed on the log risk}\\\\\n",
                        "\\multicolumn{11}{l}{ratio scale.  All other quantities were computed on the VE scale.}"
                        )
finalRes_x <- xtable::xtable(final, 
      digits=c(0,0,0,1,1,1,1,0,1,1,1,0),
      caption = "Simulation study results by scenario for standardization methods.  Results are based on $3{,}000$ replications of a simulated cohort with $\\tau=20$ time points of follow-up for $n=50{,}000$ individuals.",
      align = "rrlrrrrrrrrr"
      )
print(finalRes_x, 
      include.rownames=FALSE, 
      include.colnames=FALSE, 
      hline.after=hlines,
      file=paste0(outDir_ext, resFileName_ext), 
      sanitize.text.function=identity,
      caption.placement = "top", 
      add.to.row = addtorow)