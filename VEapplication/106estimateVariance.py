import numpy as np
from delicatessen import MEstimator
from delicatessen.utilities import inverse_logit
from delicatessen.utilities import logit

d              = r.deliData # wide data set, n rows

# Constants and option switches
trialNo        = int(r.trialNo)
numTpts        = int(r.numTimePts)
numTrials      = int(r.numTrials) 
estimator      = r.estimator
testHypothesis = r.testHypothesis

J              = numTrials - 1
tau            = numTpts - 1
Jpls1K_j       = numTrials * int(r.K_J)
meanj          = J/2 # for setj={0, ..., J} mean(setj) = {J(J+1)/2}/(J+1) = J/2

# prepare data structures
d['one']       = 1
d['zero']      = 0

N              = d.shape[0]

H              = [[] for l in range(numTpts)]
Y              = [[] for l in range(numTpts)]
D              = [[] for l in range(numTpts)]
# lpQ is the linear predictor for model with outcome variable Q
lpD            = [[] for l in range(numTpts)] 
lpH            = [[] for l in range(numTpts)]
RsupD          = [[] for l in range(numTpts)] # referred to as R^Z in the supplement
RsupH          = [[] for l in range(numTpts)]
# Dummy variables for levels of Z
Zsup1          = [[] for l in range(numTpts)]
Zsup2          = [[] for l in range(numTpts)]

# Indicators needed to construct nonadherence weights
IlagZeq0       = [[] for l in range(numTpts)]
IlagZeq1       = [[] for l in range(numTpts)]
IlagZeq2       = [[] for l in range(numTpts)]
IlLtl_23FV     = [[] for l in range(numTpts)]
IstateOfGrace  = [[] for l in range(numTpts)]
IfallFromGrace = [[] for l in range(numTpts)]

# component of at-risk indicator for outcome model
InonAdherentToA_v1 = [[] for l in range(numTpts)]

# IP weights and cumulative IP weights for treatment and dropout
IPT            = [[] for l in range(numTpts)]
IPC            = [[] for l in range(numTpts)]
Wjk            = [[] for l in range(numTpts)]

E              = [[] for j in range(numTrials)]

lpY_Aeq0       = [[] for j in range(numTrials)]
lpY_Aeq1       = [[] for j in range(numTrials)]

VE_jk          = [[] for j in range(numTrials)]
AUC_j          = [[] for j in range(numTrials)]

for j in range(numTrials):  # trial-specific 
  Wjk[j]       = [[] for k in range( ( numTpts-j ) )]
  
  lpY_Aeq0[j]  = [[] for k in range( ( numTpts-j ) )]
  lpY_Aeq1[j]  = [[] for k in range( ( numTpts-j ) )]
  
  VE_jk[j]     = [[] for k in range( ( numTpts-J ))]

def psi(theta):
  
  for j in range(numTrials):
    E[ j ]     = np.asarray( d[[ 'E_{}'.format( j ) ]])

  for l in range(numTpts):
    H[ l ]          = np.asarray( d[[ 'H_{}'.format( l ) ]])
    Y[ l ]          = np.asarray( d[[ 'Y_{}'.format( l ) ]])
    D[ l ]          = np.asarray( d[[ 'D_{}'.format( l ) ]])
    RsupD[ l ]      = np.asarray( d[[ 'RsupD_{}'.format( l )  ]]) # referred to as R^Z in the supplement
    Zsup1[ l ]      = np.asarray( d[[ 'Zsup1_{}'.format( l )  ]])
    Zsup2[ l ]      = np.asarray( d[[ 'Zsup2_{}'.format( l )  ]])
    IlagZeq0[l]     = np.asarray( d[[ 'IlagZeq0_{}'.format( l )  ]])
    IlagZeq1[l]     = np.asarray( d[[ 'IlagZeq1_{}'.format( l )  ]])
    IlagZeq2[l]     = np.asarray( d[[ 'IlagZeq2_{}'.format( l )  ]])
    
    IstateOfGrace[l]  = np.asarray( d[[ 'IstateOfGrace_{}'.format( l )  ]])
    IfallFromGrace[l] = np.asarray( d[[ 'IfallFromGrace_{}'.format( l ) ]])
    IlLtl_23FV[l]     = np.asarray( d[[ 'IlLtl_23.FV_{}'.format( l )    ]])
    
    RsupH[ l ]        = np.asarray( d[[ 'RsupH_{}'.format( l )          ]])
    InonAdherentToA_v1[ l ] = np.asarray( d[[ 'InonAdherentToA_v1_{}'.format( l )]]) 
    
    lpD[ l ] = np.asarray(d[['one',  
                    'lspl1_VUPM_{}'.format(l),
                    "lspl2_VUPM_{}".format(l),
                    "lspl3_VUPM_{}".format(l),
                    'X1_c', 'X2', 'X3', 
                    'IlagZeq1_{}'.format(l), 
                    'IlagZeq2_{}'.format(l),
                    'IpreGrace_Pfizer_{}'.format(l),    
                    'IpreGrace_Moderna_{}'.format(l),  
                    'IpreGrace_AstraZ_{}'.format(l),   
                    'Igrace_mRNA_{}'.format(l),        
                    'Igrace_AstraZ_{}'.format(l)]])
    
    lpH[ l ]   = np.asarray(d[['one', 
                    'lspl1_L2FU_{}'.format(l),
                    "lspl2_L2FU_{}".format(l), 
                    "lspl3_L2FU_{}".format(l), 
                    'Zsup1_{}'.format(l),
                    'Zsup2_{}'.format(l),
                    'X1_c', 'X2', 'X3', 'X1X3']])                    

  for j in range(numTrials):
    for k in range(1, (numTpts-j)):
      if estimator == "unstandardized":
        lpY_Aeq0[j][k] = np.asarray(d[['one',
          "oSpline_l_{}".format( ( j+k ) ),
          "oSpline_l'_{}".format( ( j+k ) ),
          "oSpline_l''_{}".format( ( j+k ) ),
          'zero',
          'zero',
          'zero',
          'zero',
          'zero',
          'zero',
          'zero']])
  
        lpY_Aeq1[j][k] = np.asarray(d[['one',
          "oSpline_l_{}".format( ( j+k ) ),
          "oSpline_l'_{}".format( ( j+k ) ),
          "oSpline_l''_{}".format( ( j+k ) ),
          'one',
          "oSpline_l_{}".format( ( j+k ) ),
          "oSpline_l'_{}".format( ( j+k ) ),
          "oSpline_l''_{}".format( ( j+k ) ),
          "oSpline_k_{}".format( k ),
          "oSpline_k'_{}".format( k ),
          "oSpline_k''_{}".format( k )]])
          
      elif estimator == "standardized":
        lpY_Aeq0[j][k] = np.asarray(d[['one',
          "oSpline_l_{}".format( ( j+k ) ),
          "oSpline_l'_{}".format( ( j+k ) ),
          "oSpline_l''_{}".format( ( j+k ) ),
          'X1_cspl1', 'X1_cspl2', 'X1_cspl3',
          'X2', 'X3', 'X1X3', 
          'zero',
          'zero',
          'zero',
          'zero',
          'zero',
          'zero',
          'zero']])
  
        lpY_Aeq1[j][k] = np.asarray(d[['one',
          "oSpline_l_{}".format( ( j+k ) ),
          "oSpline_l'_{}".format( ( j+k ) ),
          "oSpline_l''_{}".format( ( j+k ) ),
          'X1_cspl1', 'X1_cspl2', 'X1_cspl3',
          'X2', 'X3', 'X1X3', 
          'one',
          "oSpline_l_{}".format( ( j+k ) ),
          "oSpline_l'_{}".format( ( j+k ) ),
          "oSpline_l''_{}".format( ( j+k ) ),
          "oSpline_k_{}".format( k ),
          "oSpline_k'_{}".format( k ),
          "oSpline_k''_{}".format( k )]])
  ##----------------------------------------------------------------
  ##  keep track of theta indices
  ##----------------------------------------------------------------
  dParams       = theta[ 0:numDparams ]
  hi            = ( numDparams + numHparams )
  hParams       = theta[ numDparams:hi ]
  lo            = hi
  hi            = hi + numOparams
  oParams       = theta[ lo:hi ]
    
  if estimator == "standardized":
    lo            = hi
    hi            = hi + numPi0parms
    pi0params     = theta[lo:hi]
    lo            = hi
    hi            = hi + numPi1parms
    pi1params     = theta[lo:hi]

  # common to both estimators
  if testHypothesis:
    lo            = hi
    hi            = (hi + numBetaParams)
    betaParams    = theta[lo:hi]
  else:
    lo            = hi
    hi            = (hi + numlogRRparms)
    logRRparms    = theta[lo:hi]
  ##----------------------------------------------------------------
  ##  Vaccine uptake model
  ##----------------------------------------------------------------
  sumD   = np.zeros((numDparams, N))
  for l in range(numTpts):
    sumD = sumD  + ( RsupD[l] * ( D[l] - inverse_logit(np.dot(lpD[l], np.asarray(dParams).reshape(-1, 1)))) * lpD[l] ).T
  
  # predicted values for 1/(1-lambda^C)
  for l in range( numTpts ):
    p_l         = inverse_logit( np.dot( lpD[l], dParams ) ).reshape(-1,1)
    lagLambdaC  = ( RsupD[l] ) * ( 1-IlLtl_23FV[l] ) * (   
      ( IlagZeq0[l] * ( ( p_l * (1 - (Zsup1[l] + Zsup2[l] )) ) + 
      ( (1-p_l) * (Zsup1[l] + Zsup2[l] ) ) )) +
      ( IlagZeq1[l] * ( (IfallFromGrace[l] * (1-p_l)) +  
        ((1-IfallFromGrace[l]) * (1-IstateOfGrace[l]) * p_l )) ) +
      ( IlagZeq2[l] * p_l)  )                    
    IPT[ l ]    = 1 / ( 1-lagLambdaC )
  ##----------------------------------------------------------------
  ##  LTFU model
  ##----------------------------------------------------------------
  sumH   = np.zeros((numHparams, N))
  for l in range(numTpts):
    sumH = sumH + ( RsupH[l] * ( H[l] - inverse_logit(np.dot(lpH[l], np.asarray(hParams).reshape(-1, 1)))) * lpH[l] ).T
  
  for l in range(numTpts):
    hatLambdaH = inverse_logit( np.dot( lpH[l], hParams ) ).reshape(-1,1)
    IPC[ l ]   = 1 / ( 1-hatLambdaH )
  ##----------------------------------------------------------------
  ##  W_j(k)
  ##----------------------------------------------------------------
  for j in range(numTrials):
    curr             = np.ones(N).reshape(-1, 1)
    for k in range( 1, ( numTpts-j ) ):
      curr           = ( 
                       curr * 
                       E[ j ] * (
                         ( 
                           ( 1 - ( Zsup1[ j+1 ]   + Zsup2[ j+1 ] ) ) * # A_j = 0 case
                           ( 1 - ( Zsup1[ j+k ] + Zsup2[ j+k ] ) )     # adherent to A_j=0     
                         ) + ( 
                           ( Zsup1[ j+1 ] + Zsup2[ j+1 ] ) * ( 1 - ( Zsup1[ j ] + Zsup2[ j ] ) ) * # A_j = 1 case
                           ( 1-InonAdherentToA_v1[ (j+k) ] ) # adherent to A_j=1
                         )
                       ) * ( 1 - Y[ ( j+k-1 ) ] )      * # free of previous event
                       ( 1 - H[ ( j+k ) ] )            * # free of LTFU
                       IPT[ ( j+k ) ] * IPC[ ( j+k ) ] 
                       )
      Wjk[j][k]        = curr
  ##----------------------------------------------------------------
  ##  Outcome hazard model
  ##----------------------------------------------------------------
  sumO  = np.zeros((numOparams, N))
  for j in range(numTrials):
    for k in range(1, ( numTpts-j )):
      sumO = sumO + (
                    (
                    ( 1 - ( Zsup1[ j+1 ] + Zsup2[ j+1 ] ) ) * # A_j = 0 case
                    Wjk[ j ][ k ]    * # weight for treatment/LTFU history through j + k
                    ( Y[ ( j+k ) ] - inverse_logit(np.dot(lpY_Aeq0[j][k], np.asarray(oParams).reshape(-1, 1)))) * lpY_Aeq0[j][k] ).T +
                    ( 
                    ( Zsup1[ j+1 ] + Zsup2[ j+1 ] ) * ( 1 - ( Zsup1[ j ] + Zsup2[ j ] ) ) * # A_j = 1 case
                    Wjk[ j ][ k ]    * # weight for treatment/LTFU history through j + k
                    ( Y[ ( j+k ) ] - inverse_logit(np.dot(lpY_Aeq1[j][k], np.asarray(oParams).reshape(-1, 1)))) * lpY_Aeq1[j][k] ).T)
      
  if testHypothesis:
    ##----------------------------------------------------------------
    ##  beta
    ##----------------------------------------------------------------            
    beta          = np.zeros((numBetaParams, N))
    pi0out        = [[] for t in range(Jpls1K_j)] 
    pi1out        = [[] for t in range(Jpls1K_j)] 
    count         = 0 
    for j in range(numTrials):
      prod1       = np.ones(N)
      prod0       = np.ones(N)
      for k in range(1, ( numTpts-J )):  # goal is to compute VE_{j}(k) 
        prod1       = prod1 * (1-inverse_logit(np.dot(lpY_Aeq1[j][k], oParams)))
        prod0       = prod0 * (1-inverse_logit(np.dot(lpY_Aeq0[j][k], oParams)))
        if estimator=="standardized":
          pi0out[int(count)] = pi0params[int(count)] - (1-prod0)
          pi1out[int(count)] = pi1params[int(count)] - (1-prod1)
          VE_jk[j][k]   = 1 - ( pi1params[int(count)]/pi0params[int(count)] )
        else:
          VE_jk[j][k]   = 1 - ( (1 - prod1)/(1 - prod0) )
        count       = count + 1
    
    for i in range(numTrials):
      AUC_j[i] = 0 
      for k in range(1, (tau - J + 1) ): # for k in 1, 2, ..., K_J  
        AUC_j[i] = AUC_j[i] + VE_jk[i][k]
    
    sum_j_AUC_j = 0 # scalar
    for i in range(numTrials):
      sum_j_AUC_j = sum_j_AUC_j + (i * AUC_j[i])
    
    # scalar, algebraic simplification of num expression in supplement
    betaNum      = sum_j_AUC_j - ( numTrials * meanj * np.mean(AUC_j) )
    
    sum_jSquared = 0 # scalar
    for i in range(numTrials):
      sum_jSquared = sum_jSquared + (i**2)
    
    # scalar, algebraic simplification of denom expression in supplement
    betaDenom    = sum_jSquared - ( numTrials * (meanj**2) )
    
    # np.ones(N) times function of scalars
    beta[0]      = np.ones(N)*betaNum - betaDenom*betaParams[0]
    if estimator == "unstandardized":    
      return np.vstack((sumD,
                        sumH,
                        sumO,
                        beta))
    if estimator == "standardized":    
      return np.vstack((sumD,
                        sumH,
                        sumO,
                        np.vstack(pi0out),
                        np.vstack(pi1out),
                        beta))
  else:    
    ##----------------------------------------------------------------
    ##  logRR_jStar
    ##----------------------------------------------------------------
    logRRout = [[] for time in range(1, (numTpts-trialNo))]
    pi0out   = [[] for time in range(1, (numTpts-trialNo))]
    pi1out   = [[] for time in range(1, (numTpts-trialNo))]
    for time in range(1, (numTpts-trialNo)):
      prod1   = np.ones(N)
      prod0   = np.ones(N)
      for k in range(int(time)):
        prod1 = prod1 * (1-inverse_logit(np.dot(lpY_Aeq1[trialNo][k+1], oParams)))
        prod0 = prod0 * (1-inverse_logit(np.dot(lpY_Aeq0[trialNo][k+1], oParams)))
      
      if estimator == "standardized":    
        pi0out[(int(time)-1)]     = pi0params[(int(time)-1)]  - (1-prod0)
        pi1out[(int(time)-1)]     = pi1params[(int(time)-1)]  - (1-prod1)
        logRRout[(int(time)-1)]   = logRRparms[(int(time)-1)] - ( np.ones(N) * ( np.log(pi1params[(int(time)-1)]) - np.log(pi0params[(int(time)-1)]) ) )
      else:
        logRRout[(int(time)-1)]   = logRRparms[(int(time)-1)] - (np.log(1 - prod1) - np.log(1 - prod0))

    if estimator == "unstandardized":    
      return np.vstack((sumD,
                        sumH,
                        sumO,
                        np.vstack(logRRout)
                        ))
    if estimator == "standardized":    
      return np.vstack((sumD,
                        sumH,
                        sumO,
                        np.vstack(pi0out),
                        np.vstack(pi1out),
                        np.vstack(logRRout)
                        ))
                    
dInit         = r.initD
hInit         = r.initH
oInit         = r.initO
betaInit      = r.initBeta 
logRRinit     = r.initLogRR

numDparams    = len( dInit[0]  )
numHparams    = len( hInit[0]  )
numOparams    = len( oInit[0]  )
numBetaParams = len([betaInit])
numlogRRparms = len( logRRinit ) 

if estimator=="standardized":
  pi0init     = r.initPi0
  pi1init     = r.initPi1
  numPi0parms = len( pi0init )
  numPi1parms = len( pi1init )

# initial values
if estimator=="unstandardized":
  if testHypothesis:
    init = sum(dInit, []) + sum(hInit, []) + sum(oInit, []) + [betaInit]
  else:
    init = sum(dInit, []) + sum(hInit, []) + sum(oInit, []) + logRRinit

elif estimator=="standardized":
  if testHypothesis:
    init = sum(dInit, []) + sum(hInit, []) + sum(oInit, []) + pi0init + pi1init + [betaInit]
  else:
    init = sum(dInit, []) + sum(hInit, []) + sum(oInit, []) + pi0init + pi1init + logRRinit

checkIt = psi(init)
  
mest    = MEstimator(psi, init=init, subset=[(1)])
mest.estimate()

# output to R
deliEst = mest.theta
deliVar = mest.asymptotic_variance
deliCI  = mest.confidence_intervals()
bread   = mest.bread
meat    = mest.meat
