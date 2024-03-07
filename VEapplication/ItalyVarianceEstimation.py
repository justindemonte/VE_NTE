# Justin DeMonte  
# 240305

import numpy as np
import pandas as pd
from delicatessen import MEstimator
from delicatessen.utilities import inverse_logit
from delicatessen.utilities import logit
from delicatessen.estimating_equations import ee_regression

# Constants and option switches
d              = r.d
numTpts        = int(r.numTimePts)
numTrials      = int(r.numTrials) 
trialNo        = int(r.trialNo)
testHypothesis = r.testHypothesis

J              = numTrials-1
tau            = numTpts-1

meanj          = J/2 # for setj={0, ..., J} mean(setj) = {J(J+1)/2}/(J+1) = J/2

# prepare data structures
d['one']  = 1
d['zero'] = 0

censX    = [[] for j in range(numTrials)]
censY    = [[] for j in range(numTrials)]
atRisk   = [[] for j in range(numTrials)]
censPred = [[] for j in range(numTrials)]
cumPr    = [[] for j in range(numTrials)]

censQ    = [[] for j in range(numTrials)]

ocX      = [[] for j in range(numTrials)]
ocY      = [[] for j in range(numTrials)]

surv1mat = [[] for j in range(numTrials)]
surv0mat = [[] for j in range(numTrials)]
surv1    = [[] for j in range(numTrials)]
surv0    = [[] for j in range(numTrials)]

G        = [[] for j in range(numTrials)]
Z        = [[] for j in range(numTrials)]
Gcens    = [[] for j in range(numTrials)]
Zcens    = [[] for j in range(numTrials)]
Q        = [[] for j in range(numTrials)]
caltime  = [[] for j in range(numTrials)]
pi       = [[] for j in range(numTrials)]
# spline   = [[] for j in range(numTrials)]

if testHypothesis:
  VE_jk  = [[] for j in range(numTrials)]
  AUC_j  = [[] for j in range(numTrials)]
def psi(theta):
  # Ensuring correct typing
  W = np.asarray(d[['one', 'BLage_c', 'sex', 'comorbidity']])
  N = W.shape[0] # Calculate number of unique individuals
  
  for j in range(numTrials):
  
    G[j]        = np.asarray(d['G_{}'.format(j)])
    Z[j]        = np.asarray(d['Z_{}'.format(j)])

    Gcens[j]    = np.asarray(d[['G_{}'.format(j)]])
    Zcens[j]    = np.asarray(d[['Z_{}'.format(j)]])
    caltime[j]  = np.asarray(d[['one', 'trial_{}'.format(j)]])
    if testHypothesis:
      VE_jk[j]  = [[] for k in range(1, (numTpts-j+1))]

    Q[j] = np.asarray(d[['one', 
                      'prop_j_{}'.format(j),
                      "prop_j'_{}".format(j), 
                      "prop_j''_{}".format(j), 
                      'BLage_c', 'sex', 'comorbidity']])
    
    for k in range(1, (numTpts-j+1)):
      ocX[j].append(np.asarray(d[['one',
        "oc_kPlusJ_{}".format(j)   + '_{}'.format(k),
        "oc_kPlusJ'_{}".format(j)  + '_{}'.format(k),
        "oc_kPlusJ''_{}".format(j) + '_{}'.format(k),
        'Z_{}'.format(j),                          
        "oc_Z_k_{}".format(j) + '_{}'.format(k),   
        "oc_Z_k'_{}".format(j) + '_{}'.format(k),  
        "oc_Z_k''_{}".format(j) + '_{}'.format(k), 
        "oc_Z_kPlusJ_{}".format(j) + '_{}'.format(k),                    
        "oc_Z_kPlusJ'_{}".format(j) + '_{}'.format(k),                   
        "oc_Z_kPlusJ''_{}".format(j) + '_{}'.format(k)]]))            
      
      censQ[j].append(np.asarray(d[['one', # list of design matrices
        "cens_kPlusJ_{}".format(j)   + '_{}'.format(k),
        "cens_kPlusJ'_{}".format(j)  + '_{}'.format(k),
        "cens_kPlusJ''_{}".format(j) + '_{}'.format(k),
        "Z_{}".format(j),
        'BLage_c', 'sex', 'comorbidity']]))
      
      censY[j].append(np.asarray(d[["oneMinusC_{}".format(j) + "_{}".format(k)]]))
      atRisk[j].append(np.asarray(d[["Y_{}".format(j) + "_{}".format(k)]]))
          
      ocY[j].append(np.asarray(d[["E_{}".format(j) + "_{}".format(k)]]))

      surv1mat[j].append(np.asarray(d[['one',
        "oc_kPlusJ_{}".format(j) + '_{}'.format(k),
        "oc_kPlusJ'_{}".format(j) + '_{}'.format(k),
        "oc_kPlusJ''_{}".format(j) + '_{}'.format(k),
        'one',  
        "oc_k_{}".format(k),
        "oc_k'_{}".format(k),
        "oc_k''_{}".format(k),
        "oc_kPlusJ_{}".format(j) + '_{}'.format(k),
        "oc_kPlusJ'_{}".format(j) + '_{}'.format(k),
        "oc_kPlusJ''_{}".format(j) + '_{}'.format(k)]]))
      
      surv0mat[j].append(np.asarray(d[['one',
        "oc_kPlusJ_{}".format(j) + '_{}'.format(k),
        "oc_kPlusJ'_{}".format(j) + '_{}'.format(k),
        "oc_kPlusJ''_{}".format(j) + '_{}'.format(k),
        "zero", "zero", "zero", "zero", "zero", "zero", "zero"]]))
  ##----------------------------------------------------------------
  ##  keep track of theta indices
  ##----------------------------------------------------------------
  propParams      = theta[0:numPropParams]
  hi              = (numPropParams+numCensParams)
  censParams      = theta[numPropParams:hi]
  lo              = hi
  hi              = (hi + numHazParams)
  ocParams        = theta[lo:hi]
  if testHypothesis:
    lo            = hi
    hi            = (hi + numBetaParams)
    betaParams    = theta[lo:hi]
  else:
    lo            = hi
    hi            = (hi + numlogRRparms)
    logRRparms    = theta[lo:hi]
  ##----------------------------------------------------------------
  ##  Propensity score models
  ##----------------------------------------------------------------
  preds_reg = np.zeros((numPropParams, N)) 
  for j in range(numTrials):
    preds_reg = preds_reg + ( Gcens[j] * (Zcens[j] - inverse_logit(np.dot(Q[j], np.asarray(propParams).reshape(-1, 1)))) * Q[j] ).T
  
  for j in range(numTrials): # individual IPT weights
    pi[j] = inverse_logit(np.dot(Q[j], propParams)).reshape(-1, 1)
  ##----------------------------------------------------------------
  ##  Dropout models
  ##----------------------------------------------------------------
  cens_reg = np.zeros((numCensParams, N))
  for j in range(numTrials):
    for k in range(1, (numTpts-j+1)):
      cens_reg = cens_reg + ( ( (atRisk[j][k-1]) * (1-Zcens[j]) * (Gcens[j]) ) * 
        (censY[j][k-1] - inverse_logit(np.dot(censQ[j][k-1], np.asarray(censParams).reshape(-1, 1)))) * 
        censQ[j][k-1] ).T
  
  for j in range(numTrials):  # predicted values from pooled censoring model    
    for k in range(1, (numTpts-j+1)):
      censPred[j].append(inverse_logit(np.dot(censQ[j][k-1], censParams)))
  
  # Build IPC weights from predicted values
  for j in range(numTrials): 
    curr               = np.ones(N)
    for k in range(1, (numTpts - j + 1)):
      curr             = curr * censPred[j][k-1]
      cumPr[j].append(curr)
  ##----------------------------------------------------------------
  ##  Discrete hazard and time-varying VE models
  ##----------------------------------------------------------------
  outcome_reg = np.zeros((numHazParams, N)) 
  for j in range(numTrials):
    for k in range(1, (numTpts-j+1)):
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
  if testHypothesis:
    ##----------------------------------------------------------------
    ##  beta
    ##----------------------------------------------------------------            
    beta          = np.zeros((numBetaParams, N))
    
    for j in range(numTrials):
      prod1       = np.ones(N)
      prod0       = np.ones(N)
      for k in range(1, (numTpts-j+1)):  # goal is to compute VE_{j}(k) 
        prod1     = prod1 * (1-inverse_logit(np.dot(surv1mat[j][k-1], ocParams)))
        prod0     = prod0 * (1-inverse_logit(np.dot(surv0mat[j][k-1], ocParams)))
        VE_jk[j][k-1] = 1 - ((1 - prod1)/(1 - prod0))
          
    for i in range(numTrials):
      AUC_j[i] = 0 
      for k in range(1, (tau - J + 2) ): # for k in 1, 2, ..., K_J  
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
                      beta
                      ))
                      
  else:
    ##----------------------------------------------------------------
    ##  generate SE estimates and CI for logRR_{trialNo}(k) 
    ##----------------------------------------------------------------
    logRRout = [[] for seTime in range(1, (numTpts-trialNo+1))]
    for seTime in range(1, (numTpts-trialNo+1)):
      prod1   = np.ones(N)
      prod0   = np.ones(N)
      for k in range(1, (int(seTime) + 1)):
        prod1 = prod1 * (1-inverse_logit(np.dot(surv1mat[trialNo][k-1], ocParams)))
        prod0 = prod0 * (1-inverse_logit(np.dot(surv0mat[trialNo][k-1], ocParams)))

      logRRout[(int(seTime)-1)] = logRRparms[(int(seTime)-1)] - np.log(1 - prod1) - np.log(1 - prod0)

    return np.vstack((preds_reg,
                      cens_reg,
                      outcome_reg, 
                      np.vstack((logRRout))
                      ))

pInit           = r.pFit # estimates for propensity params from R
cInit           = r.cFit # estimates for censoring  params from R
oInit           = r.oFit # estimates for outcome    params from R
betaInit        = r.tFit # estimate  for beta       param  from R
logRRinit       = r.logRRinit # estimates for rho   params from R

numPropParams   = len(pInit[0])
numCensParams   = len(cInit[0])
numHazParams    = len(oInit[0])
numBetaParams   = len([betaInit])
numlogRRparms   = len(logRRinit)

if testHypothesis:
  init = sum(pInit, []) + sum(cInit, []) + sum(oInit, []) + [betaInit]
else:
  init = sum(pInit, []) + sum(cInit, []) + sum(oInit, []) + logRRinit
mest = MEstimator(psi, init=init, subset=[1])
mest.estimate()

deliEst = mest.theta
deliVar = mest.asymptotic_variance
deliCI  = mest.confidence_intervals()
