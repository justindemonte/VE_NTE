# Justin DeMonte  
# 240305
#
# display results of data analysis
rm(list=ls())
cat("\014")
library("tidyverse")
library("plotly")
library("tidyquant")
library("reshape2")
library("mgcv")
##----------------------------------------------------------------
##   Helper functions
##----------------------------------------------------------------
coarsePlot      <- function(resDir, numTrials, numTimes){
  
  pFit          <- readRDS(file = paste0(resDir, "pFit_pooled.rds"))
  cFit          <- readRDS(file = paste0(resDir, "cFit_pooled.rds"))
  oFit          <- readRDS(file = paste0(resDir, "oFit_calTime.rds"))
  
  modelParms    <- length(pFit) + length(cFit) + length(oFit)
  
  resList       <- vector(mode  = "list",  length = 3)
  resMatVE      <- matrix(nrow  = numTrials, ncol = numTimes)
  resMatCIlow   <- matrix(nrow  = numTrials, ncol = numTimes)
  resMatCIhi    <- matrix(nrow  = numTrials, ncol = numTimes)
  for (j in c(0:(numTrials-1))){
    logRRpoint  <- readRDS(file = paste0(resDir, "deliEstAnalysis", j, ".rds"))
    logRRCI     <- readRDS(file = paste0(resDir, "deliCIAnalysis",  j, ".rds"))
  
    totalParams <- length(logRRpoint)

    logRRpoint  <- logRRpoint[(modelParms+1):(totalParams)]
    
    logRRCI     <- logRRCI[(modelParms+1):(totalParams),]
    
    VE          <- 1-exp(logRRpoint) # convert from logRR to VE 
    
    VECI        <- 1-exp(logRRCI)    # convert from logRR to VE confidence limits
    CIout       <- VECI[, rev(seq_len(ncol(VECI)))] # reverse cols so confidence lims are lower, upper
    numTimePts  <- length(logRRpoint)
    k           <- c(1:numTimePts)
    
    resMatVE[j+1, k]    <- VE
    resMatCIlow[j+1, k] <- CIout[,1]
    resMatCIhi[j+1, k]  <- CIout[,2]
  }

  resList[[1]]  <- resMatVE
  resList[[2]]  <- resMatCIlow
  resList[[3]]  <- resMatCIhi
  return(resList)
}

getResTable <- function(VEres, trials, times){
  resTable  <- matrix(nrow = length(times), 
                      ncol = length(trials)+1)
  
  colnames(resTable) <- paste0("col", 1:ncol(resTable))
  
  resTable[,1] <- times
  
  jIndex <- 1
  for (j in trials){
    jIndex <- jIndex + 1
    kIndex <- 0
    for (k in times){
      kIndex <- kIndex + 1
      # j+1 because trial 0 appears in row 1 and so on
      resTable[kIndex, jIndex] <- 
        paste0(sprintf("%.0f", (VEres[[1]][(j+1), (k)] *100)),  " (", 
               sprintf("%.0f", (VEres[[2]][(j+1), (k)] *100)),  ", ",
               sprintf("%.0f", (VEres[[3]][(j+1), (k)] *100)),  ")")    
    }
  }  
  return(resTable)
}
##----------------------------------------------------------------
##   Prepare axis labels for 3-d contour plots
##----------------------------------------------------------------
yAxisLabs <- seq.Date( dmy("15/2/2021"), length=12, by='1 week' )
yAxisLabs <- str_replace_all(yAxisLabs, "2021-", "")

yAxisLabs <- seq(ymd("2021-02-15"), ymd("2021-05-03"), by = "7 days")
yAxisLabs <- as.Date(yAxisLabs)
yAxisLabs

yAxisLabs2 <- c("Feb. 15", 
                "Feb. 22",
                "Mar. 01",
                "Mar. 08",
                "Mar. 15",
                "Mar. 22",
                "Mar. 29",
                "Apr. 05",
                "Apr. 12",
                "Apr. 19",
                "Apr. 26",
                "May  03")

xAxisLabs <- c(sprintf('%#.1g',c(0,5))|> stringr::str_remove("[.]$"),
               sprintf('%#.2g',c(10,15,20,25,30,35,40))|> stringr::str_remove("[.]$"))

axx <- list(title = "Time Since Vaccination (Weeks)", 
            tickvals = c(0,50,100,150,200,250,300,350,400),
            ticktext = xAxisLabs)
axz <- list(title = "Vaccine Effectiveness", autotick = F,  
            tickvals=c(-1.8,-1.6,-1.4, -1.2, -1,-.8, -.6,-.4, -.2, 0, .2, .4, .6, .8, 1), 
            range = c(-2,1))
axy <- list(title = list(text="Calendar Date of Vaccination"),
            ticktext = yAxisLabs2, 
            tickvals = seq(1,111, by=10),
            ticketmode = "array", 
            range = c(1,119)
            , standoff=20
            )

##----------------------------------------------------------------
##   Process results
##----------------------------------------------------------------
dirCD_12_s4  <- "C:/Users/Justin/OneDrive - University of North Carolina at Chapel Hill/Documents/UNCspring2024/VEapplication/data/"
VE_CD        <- coarsePlot(dirCD_12_s4, 12, 44)

# point estimate and confidence interval for VE_0(44) for use in discussion
VE_CD[[1]][1,44]          
VE_CD[[2]][1,44]          
VE_CD[[3]][1,44]          

jSeq <- 1:120
kSeq <- 1:440

coarseMelted_CD <- vector(mode="list", length = 3)
predicted_CD    <- vector(mode="list", length = 3)
gamMod_CD       <- vector(mode="list", length = 3)
plot_CD         <- vector(mode="list", length = 3)

for (p in 1:3){
  coarseMelted_CD[[p]] <- melt(VE_CD[[p]]) %>% rename(jPlus1 = Var1) %>%
    rename(k = Var2) %>%
    rename(z = value)

  gamMod_CD[[p]] <- gam(z ~ te(jPlus1) + te(k) + ti(jPlus1, k), data=coarseMelted_CD[[p]])

  predicted_CD[[p]]                     <- matrix(data=NA, nrow=120, ncol=440)
  jIndex                                <- 0
  for (j in seq(1, 12.9, by = 0.1) ){
    jIndex                              <- jIndex + 1
    kIndex                              <- 0
    for (k in seq(1, (44-j), by = 0.1) ){
      kIndex                            <- kIndex + 1
      newdat                            <- data.frame(jPlus1=j, k=k)
      predicted_CD[[p]][jIndex, kIndex] <- predict(gamMod_CD[[p]], newdata=newdat)
    }
  }
  plot_CD[[p]] <- cbind(rep(NA, 120), predicted_CD[[p]])
}

axz <- list(title = "Vaccine Effectiveness", autotick = F,  
            tickvals=c(-1,-.8, -.6,-.4, -.2, 0, .2, .4, .6, .8, 1), 
            range = c(-1.2,1))
# 3d contour plots
fig_CD <- plot_ly()
fig_CD <- fig_CD %>% add_surface(z = ~plot_CD[[2]], opacity=.25, showscale=FALSE, cmin = -.6,cmax = 1)
fig_CD <- fig_CD %>% add_surface(z = ~plot_CD[[3]], opacity=.25, showscale=FALSE, cmin = -.6,cmax = 1)
fig_CD <- fig_CD %>% add_surface(z = ~plot_CD[[1]], opacity=1,
                                 
                           showscale=FALSE,
                           # colorbar=list(title='VE', y=.89, x=.7),
                           
                           cmin = -.6,
                           cmax = 1
                           )
fig_CD <- fig_CD %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
fig_CD

##----------------------------------------------------------------
##   Create output table for point est and 95% CI.
##----------------------------------------------------------------
outTableDir <- "C:/Users/Justin/OneDrive - University of North Carolina at Chapel Hill/Documents/VEapplication/resTable/"
resTable_CD <- getResTable(VE_CD,
  trials = c(0, 3, 6,  9),
  times  = c(1, 7, 14, 21, 28, 34))

print(xtable::xtable(resTable_CD), file=paste0(outTableDir, "table.txt"))
