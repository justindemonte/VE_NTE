###########################################################################
###
###                         Generate a cohort 
###
###########################################################################
# arbitrary large value for time-to-event variables to stand in for Inf
pseudoInf     <- 99 
minAge        <- 80 
###########################################################################
### 
###   True parameter values (constants)
### 
###########################################################################

# birth sex
sModel_int    <- -0.423166262
sModel_age    <- -0.047068132

# comorbidity 
cModel_int    <- 0.442877583 
cModel_age    <- 0.008693094
cModel_sex    <- 0.3696876
##----------------------------------------------------------------
##  Vacc uptake 
##----------------------------------------------------------------
vModel_l      <- 0.2459143
vModel_lSqr   <- -0.021850663
vModel_int    <- -2.644752947
vModel_age    <- -0.052336868
vModel_sex    <- 0.029531031
vModel_com    <- -0.04762281
##----------------------------------------------------------------
##  Events
##----------------------------------------------------------------
# time-varying intercept for treated and untreated in outcome haz model 
if ( scenario==1 ){ tvi <- c( -4,    0,     0, -2.5,   0,    0,  .02, .005 )}
if ( scenario==2 ){ tvi <- c( -4, -.01, -.003, -2.5, .02,  .006,   0,    0 )}
if ( scenario==3 ){ tvi <- c( -4, -.01, -.003, -2.5, .02,  .006, .02, .005 )}

alphaInt      <- tvi[1]
alpha_l       <- tvi[2]
alpha_lSqr    <- tvi[3]  
# z=1
alphaZint     <- tvi[4]  
alphaZ_l      <- tvi[5]  
alphaZ_lSqr   <- tvi[6]   
alphaZ_k      <- tvi[7]  
alphaZ_kSqr   <- tvi[8]  

# covariate main effects for outcome hazard model
alphaAge      <- -0.012623418
alphaSex      <- -0.254764979
alphaCom      <- 0.424559332

# covariate interaction term for outcome hazard model
alphaZage     <- ageProdTerm
###########################################################################
###
###  Build simulated cohort
###
###########################################################################
### Covariates
simd <- tibble( id = c( 1:N ) ) %>%
    
  # age 
  add_column( X1 = VGAM::rfoldnorm(n = nrow(.),
      mean = 0, sd = 7, a1 = 1, a2 = 1) + minAge )  %>%
  mutate( X1_c = X1 - ( minAge + (7*sqrt(2/pi)) ) ) %>%  # centered on population mean
    
  # birth sex
  mutate( prX2 = plogis( sModel_int + ( sModel_age*X1_c ) ) )  %>%
  add_column( runifX2 = runif(nrow(.)) ) %>%
  mutate( X2   = if_else( (runifX2 < prX2), 1, 0) ) %>%
    
  # comorbidity 
  mutate( prX3 = plogis( cModel_int + ( cModel_age*X1_c ) + ( cModel_sex*X2 ) )) %>%
  add_column( runifX3 = runif(nrow(.)) ) %>%
  mutate( X3   = if_else( (runifX3 < prX3), 1, 0) ) 
##----------------------------------------------------------------
##  Generate cfls for the cohort
##----------------------------------------------------------------
# initialize event time to large number (for practical purposes Infty)
simd <- simd %>% mutate( TsupV1gtTau = pseudoInf ) 

# counterfactual event time T^{v_1 > \tau} for individuals in the cohort 
for ( l in 1:tau ){ # l indexes calendar time
  
  simd <- simd %>% mutate(condProb = plogis( # logit^{-1}(linear predictor)
    alphaInt       
    + ( alpha_l    * l )     
    + ( alpha_lSqr * ( l**2 ) )        
    + ( alphaAge   * X1_c )         
    + ( alphaSex   * X2 )    
    + ( alphaCom   * X3 )    
  )) %>%
    # event at time l ? 
    add_column( runifY = runif(nrow(.)) ) %>%
    mutate( Y = if_else( runifY < condProb, 1, 0 ) ) %>%
    
    # update TsupV1gtTau 
    mutate( TsupV1gtTau = case_when(
      TsupV1gtTau < l   ~ TsupV1gtTau,   # already had event prior to l
      Y == 1            ~ as.numeric(l), # event at time l
      TRUE ~ pseudoInf)                  # no event yet as of l
    )  
  simd <- simd %>%
    mutate( "oneMinuslambda0_{l}" := 1-condProb ) %>% 
    select( -c("Y", "runifY") )
}
# counterfactual event time T^{v_1} for individuals in the cohort 
for (v in 1:( J+1 ) ){       
  
  # initialize event time under no anticipation effect assumption
  # counterfactuals under strategy "get vacc at v_1" before v_1
  # involve same hazard as counterfactuals under "never get vaccinated" 
  simd <- simd %>% mutate( eventTime := TsupV1gtTau ) 
  
  for ( l in v:tau ){ # l indexes calendar time
    
    simd <- simd %>% mutate( condProb = plogis( # logit^{-1}(linear predictor)
      alphaInt       
      + ( alpha_l     * l )     
      + ( alpha_lSqr  * (l**2) )        
      + ( alphaAge    * X1_c )  
      + ( alphaSex    * X2 )    
      + ( alphaCom    * X3 )    
      
      + alphaZint 
      + ( alphaZ_l    * l )                 
      + ( alphaZ_lSqr * (l**2) )      
      + ( alphaZ_k    * ( l-v+1 ) )         
      + ( alphaZ_kSqr * ( ( l-v+1 )**2 ) )     
      
      + ( alphaZage   * X1_c )
    )) %>%
      
      # event at time l ? 
      add_column( runifY = runif(nrow(.)) ) %>%
      mutate( Y = if_else( runifY < condProb, 1, 0 ) ) %>%
      
      # update T^{v_1 = l}  
      mutate( eventTime = case_when(
        eventTime < l   ~ eventTime,     # already had event prior to l
        Y == 1          ~ as.numeric(l), # event at time l
        TRUE            ~ pseudoInf)     # no event yet as of l
      )
    simd <- simd %>%
      mutate( "oneMinuslambdaSup{v}_{l}" := 1-condProb ) %>% 
      select( -c("Y", "runifY") )
  }
  simd <- simd %>% rename("TsupV1eq{v}" := eventTime)
}
##----------------------------------------------------------------
##  Vaccine uptake 
##----------------------------------------------------------------
# initialize vacc uptake time to large number (for practical purposes Infty)
simd <- simd %>% mutate(V_1 = pseudoInf)
for (l in 0 : tau){ # l indexes calendar time
  
  simd <- simd %>%
    mutate(condProb = plogis( # logit^{-1}(linear predictor)
      vModel_int 
      + ( vModel_l     * l )         
      + ( vModel_lSqr  * (l**2) )         
      + ( vModel_age   * X1_c ) 
      + ( vModel_sex   * X2 )
      + ( vModel_com   * X3 )
    )) %>%  
    
    # vacc uptake at calendar time l+1 ? 
    add_column( runifA = runif(nrow(.)) ) %>%
    mutate( A = if_else( runifA < condProb, 1, 0 ) ) %>%
    
    # update V_1
    mutate( V_1     = case_when(
      V_1 < ( l+1 ) ~ V_1,               # already got vacced prior to l+1
      A == 1        ~ as.numeric( l+1 ), # vacced at time l+1
      TRUE          ~ pseudoInf)         # no vacc yet as of l+1
    )
  
  simd <- simd %>% select(-c("A", "runifA", "condProb"))
}
##----------------------------------------------------------------
##  Observed events and eligibility indicators
##----------------------------------------------------------------
# derive observed outcome from vacc uptake time and corresponding counterfactual
simd["Tstar"] <- simd["TsupV1gtTau"] * ( simd["V_1"] > ( J+1 ) )
for ( l in 1:( J+1 ) ){
  
  simd["Tstar"] <- simd["Tstar"] + 
    ( simd[ paste0("TsupV1eq", l) ]  * ( simd["V_1"] == l ) ) 
}
# derive eligibility indicators for each calendar time point
for ( l in 0 : tau ){
  
  simd[ paste0("E_", l) ] <- ( simd["V_1"]   > l ) *  # aka Z_l=0 
    ( simd["Tstar"] > l )               # free of event at time l
}