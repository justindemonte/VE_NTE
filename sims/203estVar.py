import numpy as np
from delicatessen import MEstimator
from delicatessen.utilities import inverse_logit
from delicatessen.utilities import logit

d           = r.deliData # wide data set, n rows

# Constants and option switches
numTpts     = int(r.numTimePts)
numTrials   = int(r.numTrials) 
J           = numTrials-1
tau         = numTpts-1
scenario    = r.scenario
simulation  = r.simulation
analysis    = r.analysis
spl         = r.spl
stdMethod   = simulation=="extension" and analysis=="proposed"

k_0         = int(r.k_0) 
j_0         = int(r.j_0) 

meanj       = J/2 # for setj={0, ..., J} mean(setj) = {J(J+1)/2}/(J+1) = J/2

# time points at which VE will be evaluated on cal.- and trial-time scales
seTimes_j     = r.seTimes_j
seTimes_k     = r.seTimes_k
numCalTimes   = len(seTimes_j)
numTrialTimes = len(seTimes_k)

# prepare data structures
d['one']    = 1
d['zero']   = 0

N           = d.shape[0]

# lpQ is linear predictor for the model with Q as outcome
Y           = [[] for l in range(numTpts)]
lagY        = [[] for l in range(numTpts)]
leadZ       = [[] for l in range(numTpts)]
Z           = [[] for l in range(numTpts)]
RsupLeadZ   = [[] for l in range(numTpts)]
lpZ         = [[] for l in range(numTpts)]

# IPT weight and cumulative IPT weight for unstandardized estimator
PrZ_x       = [[] for j in range(numTrials)] 
Wjk         = [[] for j in range(numTrials)] 

E           = [[] for j in range(numTrials)] 

logRR_jOut  = [[] for ct in range(numCalTimes)]
logRR_kOut  = [[] for k in range(numTrialTimes)]

pi0out      = [[] for ct in range(numCalTimes)]
pi1out      = [[] for ct in range(numCalTimes)]

lpY_Aeq0    = [[] for j in range(numTrials)]
lpY_Aeq1    = [[] for j in range(numTrials)]

VE_jk       = [[] for j in range(numTrials)]
AUC_j       = [[] for j in range(numTrials)]

for j in range( numTrials ):  # trial-specific quantities
  PrZ_x[j]  = [[] for k in range( ( numTpts-j-1 ))]
  Wjk[j]    = [[] for k in range( ( numTpts-j ) )]
  
  lpY_Aeq0[j] = [[] for k in range( ( numTpts-j ) )]
  lpY_Aeq1[j] = [[] for k in range( ( numTpts-j ) )]
  
  VE_jk[j]  = [[] for k in range( ( numTpts-j ))]

def psi(theta):

  for j in range( numTrials ):
    E[ j ]         = np.asarray( d[[ 'E_{}'.format( j ) ]])
    
  for l in range( numTpts ):
    Y[ l ]         = np.asarray( d[[ 'Y_{}'.format( l )     ]])
    lagY[ l ]      = np.asarray( d[[ 'lagY_{}'.format( l )  ]])
    leadZ[ l ]     = np.asarray( d[[ 'leadZ_{}'.format( l ) ]])
    Z[ l ]         = np.asarray( d[[ 'Z_{}'.format( l )     ]])
    RsupLeadZ[ l ] = np.asarray( d[[ 'RsupLeadZ_{}'.format( l )  ]])
      
    if spl:
      lpZ[ l ]  = np.asarray( d[[ 'one', 
        "Z_lspl1_{}".format(l), 
        "Z_lspl2_{}".format(l),
        "Z_lspl3_{}".format(l),
        'X1_c', 'X2', 'X3' ]])
    else:  
      lpZ[ l ]  = np.asarray( d[[ 'one', 
        'i_{}'.format(l), 'iSqr_{}'.format(l),
        'X1_c', 'X2', 'X3' ]])

  # Design matrix for outcome haz model  
  for j in range( numTrials ):
    for k in range( 1, (numTpts-j) ):
      if stdMethod:
        lpY_Aeq0[j][k] = np.asarray(d[['one',
          'i_{}'.format( ( j+k ) ),
          'iSqr_{}'.format( ( j+k ) ),
          'X1_c', 'X2', 'X3', 
          'zero',
          'zero',
          'zero',
          'zero',
          'zero',
          'zero']])
    
        lpY_Aeq1[j][k] = np.asarray(d[['one',
          'i_{}'.format( ( j+k ) ),
          'iSqr_{}'.format( ( j+k ) ),
          'X1_c', 'X2', 'X3', 
          'one',
          'i_{}'.format( ( j+k ) ),
          'iSqr_{}'.format( ( j+k ) ),
          'i_{}'.format( k ),
          'iSqr_{}'.format( k ),
          'X1_c']])
  
      elif simulation=="main" and analysis=="naive":
        lpY_Aeq0[j][k] = np.asarray(d[['one',
          'i_{}'.format( ( j+k ) ),
          'iSqr_{}'.format( ( j+k ) ),
          'zero',
          'zero',
          'zero']])

        lpY_Aeq1[j][k] = np.asarray(d[['one',
          'i_{}'.format( ( j+k ) ),
          'iSqr_{}'.format( ( j+k ) ),
          'one',
          'i_{}'.format( k ),
          'iSqr_{}'.format( k )]])
          
      elif spl:
        lpY_Aeq0[j][k] = np.asarray(d[['one',
          "oSpline_l_{}".format(   ( j+k ) ),
          "oSpline_l'_{}".format(  ( j+k ) ),
          "oSpline_l''_{}".format( ( j+k ) ),
          'zero',
          'zero',
          'zero',
          'zero',
          'zero',
          'zero',
          'zero']])
              
        lpY_Aeq1[j][k] = np.asarray(d[['one',
          "oSpline_l_{}".format(   ( j+k ) ),
          "oSpline_l'_{}".format(  ( j+k ) ),
          "oSpline_l''_{}".format( ( j+k ) ),
          'one',
          "oSpline_l_{}".format(   ( j+k ) ),
          "oSpline_l'_{}".format(  ( j+k ) ),
          "oSpline_l''_{}".format( ( j+k ) ),
          "oSpline_k_{}".format(   k ),
          "oSpline_k'_{}".format(  k ),
          "oSpline_k''_{}".format( k )]])  
      else: 
        lpY_Aeq0[j][k] = np.asarray(d[['one',
          'i_{}'.format( ( j+k ) ),
          'iSqr_{}'.format( ( j+k ) ),
          'zero',
          'zero',
          'zero',
          'zero',
          'zero']])
    
        lpY_Aeq1[j][k] = np.asarray(d[['one',
          'i_{}'.format( ( j+k ) ),
          'iSqr_{}'.format( ( j+k ) ),
          'one',
          'i_{}'.format( ( j+k ) ),
          'iSqr_{}'.format( ( j+k ) ),
          'i_{}'.format( k ),
          'iSqr_{}'.format( k )]])
  ##----------------------------------------------------------------
  ##  keep track of theta indices
  ##----------------------------------------------------------------
  zParams       = theta[ 0:numZparams ]
  hi            = ( numZparams + numOcParams )
  ocParams      = theta[ numZparams:hi ]
  
  if stdMethod:
    lo            = hi
    hi            = hi + numPi0parms
    pi0params     = theta[lo:hi]
    lo            = hi
    hi            = hi + numPi1parms
    pi1params     = theta[lo:hi]
  if not stdMethod and analysis=="proposed":
    lo            = hi
    hi            = hi + numBetaParams
    betaParams    = theta[lo:hi]
    
  lo            = hi
  hi            = hi + numlogRRparms
  logRRparms    = theta[lo:hi]
  ##----------------------------------------------------------------
  ##  Propensity score model
  ##----------------------------------------------------------------
  sumZ   = np.zeros((numZparams, N))
  for l in range( numTpts ):
    sumZ = sumZ + ( RsupLeadZ[l] * ( leadZ[l] - inverse_logit(np.dot(lpZ[l], np.asarray(zParams).reshape(-1, 1)))) * lpZ[l] ).T
  
  for j in range( numTrials ):
    for k in range( numTpts-j-1 ):
      propjk      = inverse_logit( np.dot( lpZ[ ( j+k ) ], zParams ) ).reshape(-1,1)
      PrZ_x[j][k] = ( Z[ ( j+k ) ] +           # weight is one if already vacced
                      ( ( 1 - Z[ ( j+k ) ] ) * # ow weight depends on leadZ:
                      ( ( leadZ[ (j+k) ] * ( 1 / propjk ) ) + ( ( 1 - leadZ[ (j+k) ] ) * ( 1 / ( 1 - propjk ) ) ) )
                      )
                    )
  for j in range( numTrials ):
    curr              = np.ones(N).reshape(-1,1)
    for k in range( ( numTpts-j-1 ) ):
      curr            = curr * PrZ_x[j][k]
      Wjk[j][ (k+1) ] = curr
  ##----------------------------------------------------------------
  ##  Outcome hazard model
  ##----------------------------------------------------------------
  outcome_reg  = np.zeros((numOcParams, N))
  for j in range( numTrials ):
    for k in range( 1, (numTpts-j) ):
      outcome_reg = outcome_reg + (
        ( E[ j ] * ( 1 - leadZ[ j ] ) * # A_j = 0 case
        ( 1 - Z[ ( j+k ) ] )        * # adherent to trial-j regimen (A_j=0)
        ( 1 - lagY[ ( j+k ) ] )     * # free of previous event
        Wjk[ j ][ k ]               * 
        ( Y[ ( j+k ) ] - inverse_logit(np.dot(lpY_Aeq0[j][k], np.asarray(ocParams).reshape(-1, 1)))) * lpY_Aeq0[j][k] ).T +
        ( E[j] * leadZ[j]           * # A_j = 1 case
        Z[ ( j+k ) ]                * # adherent to trial-j regimen (A_j=1)
        ( 1 - lagY[ ( j+k ) ] )     * # free of previous event
        Wjk[ j ][ k ]               * 
        ( Y[ ( j+k ) ] - inverse_logit(np.dot(lpY_Aeq1[j][k], np.asarray(ocParams).reshape(-1, 1)))) * lpY_Aeq1[j][k] ).T)
  ##----------------------------------------------------------------
  ##  logRR_j 
  ##----------------------------------------------------------------
  count_j      = 0
  for j in range( numTrials ):
    sTimej     = seTimes_j[count_j]
    prod1      = np.ones(N)
    prod0      = np.ones(N)
    count_k    = 0
    for k in range( 1, (numTpts-j) ):  # goal is to compute logRR_{j}(k)
      sTimek   = seTimes_k[count_k]
      prod1    = prod1 * (1-inverse_logit(np.dot(lpY_Aeq1[j][k], ocParams)))
      prod0    = prod0 * (1-inverse_logit(np.dot(lpY_Aeq0[j][k], ocParams)))
      if analysis=="proposed":
        VE_jk[j][k] = 1 - ((1 - prod1)/(1 - prod0))

      if k_0==k and sTimej==j:
        if stdMethod:
          pi0out[count_j]     = pi0params[count_j]  - (1-prod0)
          pi1out[count_j]     = pi1params[count_j]  - (1-prod1)
          logRR_jOut[count_j] = logRRparms[count_j] - ( np.ones(N) * ( np.log(pi1params[count_j]) - np.log(pi0params[count_j]) ) )
        else:
          logRR_jOut[count_j] = logRRparms[count_j] - ( np.log(1 - prod1) - np.log(1 - prod0) )
        if count_j < (numCalTimes-1):
          count_j  = count_j + 1
      
      if simulation=="main":
        if j_0==j and sTimek==k:
          logRR_kOut[count_k] = logRRparms[(numCalTimes+count_k)] - (np.log(1 - prod1) - np.log(1 - prod0))
          if count_k < (numTrialTimes-1):
            count_k  = count_k + 1
  
  if not stdMethod and analysis=="proposed": 
    ##----------------------------------------------------------------
    ##  beta
    ##----------------------------------------------------------------            
    beta          = np.zeros((numBetaParams, N))
    
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
    return np.vstack((sumZ,
                    outcome_reg,
                    beta,
                    np.vstack(logRR_jOut),
                    np.vstack(logRR_kOut)
                    ))
  elif stdMethod: # extension, proposed
    return np.vstack((sumZ,
                    outcome_reg, 
                    np.vstack(pi0out),
                    np.vstack(pi1out),
                    np.vstack(logRR_jOut)
                    ))
  elif simulation=="extension": # extension, naive   
    return np.vstack((sumZ,
                    outcome_reg,
                    np.vstack(logRR_jOut)
                    ))
  else: # main, naive   
    return np.vstack((sumZ,
                    outcome_reg,
                    np.vstack(logRR_jOut),
                    np.vstack(logRR_kOut)
                    ))

# initial values passed in from R                    
propInit        = r.initProp    # estimates of kappa params from R
oInit           = r.initOutcome # estimates of alpha params from R
logRRinit       = r.initLogRR   # estimates of log risk ratio params from R
if not stdMethod and analysis=="proposed": 
  betaInit      = r.initBeta 
if stdMethod:
  pi0init       = r.initPi0
  pi1init       = r.initPi1

numZparams      = len( propInit[0]  )
numOcParams     = len( oInit[0]     )
numlogRRparms   = len( logRRinit    )
if not stdMethod and analysis=="proposed": 
  numBetaParams = len([betaInit])
if stdMethod:
  numPi0parms   = len( pi0init      )
  numPi1parms   = len( pi1init      )

# initial values by simulation/analysis
if stdMethod: 
  init = sum(propInit, []) + sum(oInit, []) + pi0init + pi1init + logRRinit 
elif analysis=="proposed":
  init = sum(propInit, []) + sum(oInit, []) + [betaInit] + logRRinit
else:
  init = sum(propInit, []) + sum(oInit, []) + logRRinit

mest = MEstimator(psi, init=init, subset=[(1)])
mest.estimate()

# output returned to R program. 
deliEst = mest.theta
deliVar = mest.asymptotic_variance
deliCI  = mest.confidence_intervals()
bread   = mest.bread
meat    = mest.meat
