# Justin DeMonte
# 240305

# This copy is for obtaining multiple SE and CI estimates
import numpy as np
import pandas as pd
from delicatessen import MEstimator
from delicatessen.utilities import inverse_logit
from delicatessen.utilities import logit
from delicatessen.estimating_equations import ee_regression

d               = r.d # wide data set

# Constants and option switches
numTpts         = int(r.numTimePts)
numTrials       = int(r.numTrials) 
J               = numTrials-1
tau             = numTpts-1
scenario        = r.scenario

k_0             = int(r.k_0) 
j_0             = int(r.j_0) 

meanj           = J/2 # for setj={0, ..., J} mean(setj) = {J(J+1)/2}/(J+1) = J/2

# time points at which VE will be evaluated on cal.- and trial-time scales
seTimes_j       = r.seTimes_j
seTimes_k       = r.seTimes_k
numCalTimes     = len(seTimes_j) 
numTrialTimes   = len(seTimes_k)

# prepare data structures
d['one']        = 1
d['zero']       = 0

N               = d.shape[0]

G               = [[] for j in range(numTrials)]
Z               = [[] for j in range(numTrials)]
Xprop           = [[] for j in range(numTrials)]
pi              = [[] for j in range(numTrials)]

# censX    = [[] for j in range(numTrials)]
Gcens           = [[] for j in range(numTrials)]
Zcens           = [[] for j in range(numTrials)]
censQ           = [[] for j in range(numTrials)]
censY           = [[] for j in range(numTrials)]
atRisk          = [[] for j in range(numTrials)]
censPred        = [[] for j in range(numTrials)]
cumPr           = [[] for j in range(numTrials)]

ocX             = [[] for j in range(numTrials)]
ocY             = [[] for j in range(numTrials)]
surv0mat        = [[] for j in range(numTrials)]
surv1mat        = [[] for j in range(numTrials)]

VE_jk           = [[] for j in range(numTrials)]
logRR_jOut      = [[] for j in range(numCalTimes)]
logRR_kOut      = [[] for k in range(numTrialTimes)]
AUC_j           = [[] for j in range(numTrials)]
def psi(theta):
  
  for j in range(numTrials):
    G[j]        = np.asarray(d['G_{}'.format(j)])
    Z[j]        = np.asarray(d['Z_{}'.format(j)])
    
    Gcens[j]    = np.asarray(d[['G_{}'.format(j)]])
    Zcens[j]    = np.asarray(d[['Z_{}'.format(j)]])
    
    Xprop[j]    = np.asarray(d[['one', 
      'prop_j_{}'.format(j),
      "prop_j'_{}".format(j), 
      "prop_j''_{}".format(j), 
      'L1_c', 'L2', 'L3']])
  
  for j in range(numTrials):
    censQ[j]    = [[] for k in range(1, (numTpts-j))]
    censY[j]    = [[] for k in range(1, (numTpts-j))]
    atRisk[j]   = [[] for k in range(1, (numTpts-j))]
    ocX[j]      = [[] for k in range(1, (numTpts-j))]
    surv1mat[j] = [[] for k in range(1, (numTpts-j))]
    surv0mat[j] = [[] for k in range(1, (numTpts-j))]
    ocY[j]      = [[] for k in range(1, (numTpts-j))]
    VE_jk[j]    = [[] for k in range(1, (numTpts-j))]
    for k in range(1, (numTpts-j)):
      censQ[j][k-1]  = np.asarray(d[['one', # list of design matrices 
        "cens_kPlusJ_{}".format(j)   + '_{}'.format(k),
        "cens_kPlusJ'_{}".format(j)  + '_{}'.format(k),
        "cens_kPlusJ''_{}".format(j) + '_{}'.format(k),
        'L1_c', 'L2', 'L3']])
      censY[j][k-1]  = np.asarray(d[["oneMinusC_{}".format(j) + "_{}".format(k)]])
      atRisk[j][k-1] = np.asarray(d[["Y_{}".format(j) + "_{}".format(k)]])
      
      ### outcome hazard model design matrices
      ocX[j][k-1] = np.asarray(d[['one', 
        "oc_kPlusJ_{}".format(j)   + '_{}'.format(k),
        "oc_kPlusJ'_{}".format(j)  + '_{}'.format(k),
        "oc_kPlusJ''_{}".format(j) + '_{}'.format(k),
        'Z_{}'.format(j),                          
        "oc_Z_k_{}".format(j) + '_{}'.format(k),   
        "oc_Z_k'_{}".format(j) + '_{}'.format(k),  
        "oc_Z_k''_{}".format(j) + '_{}'.format(k), 
        "oc_Z_kPlusJ_{}".format(j) + '_{}'.format(k),                 
        "oc_Z_kPlusJ'_{}".format(j) + '_{}'.format(k),                
        "oc_Z_kPlusJ''_{}".format(j) + '_{}'.format(k)]])
        
      surv1mat[j][k-1] = np.asarray(d[['one', 
        "oc_kPlusJ_{}".format(j) + '_{}'.format(k),
        "oc_kPlusJ'_{}".format(j) + '_{}'.format(k),
        "oc_kPlusJ''_{}".format(j) + '_{}'.format(k),
        'one',  
        "oc_k_{}".format(k),
        "oc_k'_{}".format(k),
        "oc_k''_{}".format(k),
        "oc_kPlusJ_{}".format(j) + '_{}'.format(k),
        "oc_kPlusJ'_{}".format(j) + '_{}'.format(k),
        "oc_kPlusJ''_{}".format(j) + '_{}'.format(k)]])
        
      surv0mat[j][k-1] = np.asarray(d[['one', 
        "oc_kPlusJ_{}".format(j) + '_{}'.format(k),
        "oc_kPlusJ'_{}".format(j) + '_{}'.format(k),
        "oc_kPlusJ''_{}".format(j) + '_{}'.format(k),
        "zero", "zero", "zero", "zero", "zero", "zero", "zero"]])
    
      ocY[j][k-1] = np.asarray(d[["E_{}".format(j) + "_{}".format(k)]])
  ##----------------------------------------------------------------
  ##  keep track of theta indices
  ##----------------------------------------------------------------
  propParams    = theta[0:numPropParams]
  hi            = (numPropParams + numCensParams)
  censParams    = theta[numPropParams:hi]
  lo            = hi
  hi            = (hi + numHazParams)
  ocParams      = theta[lo:hi]
  lo            = hi
  hi            = hi + numBetaParams
  betaParams    = theta[lo:hi]
  lo            = hi
  hi            = hi + numlogRRparms
  logRRparms    = theta[lo:hi]
  ##----------------------------------------------------------------
  ##  Propensity score models
  ##----------------------------------------------------------------
  preds_reg     = np.zeros((numPropParams, N))
  for j in range(numTrials):
    preds_reg   = preds_reg + ( Gcens[j] * (Zcens[j] - inverse_logit(np.dot(Xprop[j], np.asarray(propParams).reshape(-1, 1)))) * Xprop[j] ).T
  
  for j in range(numTrials): # individual IPT weights
    pi[j]       = inverse_logit(np.dot(Xprop[j], propParams)).reshape(-1, 1)
  ##----------------------------------------------------------------
  ##  Dropout models
  ##----------------------------------------------------------------
  cens_reg      = np.zeros((numCensParams, N)) 
  for j in range(numTrials):
    for k in range(1, (numTpts-j)):
      cens_reg  = cens_reg + ( ( (atRisk[j][k-1]) * (1-Zcens[j]) * (Gcens[j]) ) * 
        (censY[j][k-1] - inverse_logit(np.dot(censQ[j][k-1], np.asarray(censParams).reshape(-1, 1)))) * 
        censQ[j][k-1] ).T
  
  for j in range(numTrials):  # predicted values from pooled censoring model  
    censPred[j]        = [[] for k in range(1, (numTpts-j))]
    for k in range(1, (numTpts-j)):
      censPred[j][k-1] = inverse_logit(np.dot(censQ[j][k-1], censParams))
    
  # Build IPC weights from predicted values
  for j in range(numTrials):
    cumPr[j]           = [[] for k in range(1, (numTpts-j))]
    curr               = np.ones(N)
    for k in range(1, (numTpts-j)):
      curr             = curr * censPred[j][k-1]
      cumPr[j][k-1]    = curr
  ##----------------------------------------------------------------
  ##  Hazard models
  ##----------------------------------------------------------------
  outcome_reg   = np.zeros((numHazParams, N))
  for j in range(numTrials):
    for k in range(1, (numTpts-j)):
      outcome_reg = outcome_reg + ( 
        ( Gcens[j] * atRisk[j][k-1] * ( 
        ( Zcens[j] / pi[j] ) +
        ( ( (1-Zcens[j]) * censY[j][k-1] ) / ( (1-pi[j]) * cumPr[j][k-1].reshape(-1, 1) ) ) 
        ) ) *
        ( ocY[j][k-1] - inverse_logit(np.dot(ocX[j][k-1], np.asarray(ocParams).reshape(-1, 1))) ) *
        ocX[j][k-1] 
        ).T
  ##----------------------------------------------------------------
  ##  logRR_j(k) and beta
  ##----------------------------------------------------------------
  beta         = np.zeros((numBetaParams, N))
  
  count_j      = 0
  for j in range(numTrials):
    sTimej     = seTimes_j[count_j]
    prod1      = np.ones(N)
    prod0      = np.ones(N)
    count_k    = 0
    for k in range(1, (numTpts-j)):  # goal is to compute logRR_{j}(k) 
      sTimek   = seTimes_k[count_k]
      prod1    = prod1 * (1-inverse_logit(np.dot(surv1mat[j][k-1], ocParams)))
      prod0    = prod0 * (1-inverse_logit(np.dot(surv0mat[j][k-1], ocParams)))
      
      VE_jk[j][k-1] = 1 - ((1 - prod1)/(1 - prod0))
      
      if k_0==k and sTimej==j:
        logRR_jOut[count_j] = logRRparms[count_j] - (np.log(1 - prod1) - np.log(1 - prod0))
        if count_j < numCalTimes-1:
          count_j  = count_j + 1
    
      if j_0==j and sTimek==k:
        logRR_kOut[count_k] = logRRparms[(numCalTimes+count_k)] - (np.log(1 - prod1) - np.log(1 - prod0))
        if count_k < numTrialTimes-1:
          count_k  = count_k + 1
          
  for i in range(numTrials):
    AUC_j[i] = 0 
    for k in range(1, (tau - J + 1) ): # for k in 1, 2, ..., K_J  
      AUC_j[i] = AUC_j[i] + VE_jk[i][k-1]
    
  sum_j_AUC_j = 0 # scalar
  for i in range(numTrials):
    sum_j_AUC_j = sum_j_AUC_j + (i * AUC_j[i])
  
  # scalar
  betaNum      = sum_j_AUC_j - ( numTrials * meanj * np.mean(AUC_j) )
  
  sum_jSquared = 0 # scalar
  for i in range(numTrials):
    sum_jSquared = sum_jSquared + (i**2)
  
  # scalar
  betaDenom    = sum_jSquared - ( numTrials * (meanj**2) )
  
  # np.ones(N) times function of scalars
  beta[0]      = np.ones(N)*betaNum - betaDenom*betaParams[0]
  return np.vstack((preds_reg,
                    cens_reg,
                    outcome_reg, 
                    beta,
                    np.vstack(logRR_jOut),
                    np.vstack(logRR_kOut)
                    ))
    
pInit           = r.initProp    # estimates for propensity params from R
cInit           = r.initCens    # estimates for censoring  params from R
oInit           = r.initOutcome # estimates for outcome    params from R
betaInit        = r.initBeta    # estimate  for beta       param  from R
logRRinit       = r.initLogRR   # estimates for log risk ratio params from R

numPropParams   = len(pInit[0])
numCensParams   = len(cInit[0])
numHazParams    = len(oInit[0])
numBetaParams   = len([betaInit])
numlogRRparms   = len(logRRinit)

init = sum(pInit, []) + sum(cInit, []) + sum(oInit, []) + [betaInit] + logRRinit
mest = MEstimator(psi, init=init, subset=[(1)])
mest.estimate()

deliEst = mest.theta
deliVar = mest.asymptotic_variance
deliCI  = mest.confidence_intervals()
