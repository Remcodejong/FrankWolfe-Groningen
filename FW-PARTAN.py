import numpy as np
from IPython.core.magics.osm import line_cell_magic
from inspect import isfunction
import matplotlib.pyplot as plt
import sympy as sym
import copy
import time
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
from graph import *
from otherFun import *
from paths import *


def findTauMin(lcur,lprev,tprev,tmin):
    bar = lambda X: 1-X
    if tprev<0:
        return 1+1/(bar(lprev)*bar(lcur)*(1-tprev/tmin)-1)
    if tprev >= 0:
        return 1+1/(bar(lprev)*bar(lcur)*bar(tprev)-1)
    
def computeDualGap(G,F):
    dualGap = 0
    for i in range(len(G.graph)):
        for j in range(len(G.graph)):
            dualGap += G    .costMatrix[i][j](F[i][j])*F[i][j]
    return dualGap

def FrankWolfePARTAN(graph,P,r,maxIt=10000,itTol=1e-2,LSTol=1e-3):
  errArr = np.array([])
  
  n = 0
  
  V_n = initialSetup(graph,P,r)
  V_prev = V_n
  err = 100000
  W_prev = V_n
  W_prev2 = V_n
  W_n = V_n
 
  l = 1000
  lprev = 1000
  
  tau = 0
  tauprev = 0
  taumin = 0
  
  while(abs(np.max(err)) and n < maxIt):
      n = n + 1
      newCost = updateCostMatrix(graph,W_n)
      
      E_n = AONUpdate(graph,P,newCost,r)
      
      lineSearchFn = lambda t: sum([sum([graph.costMatrix[i][j](W_n[i][j]+t*(E_n[i][j]-W_n[i][j]))* (E_n[i][j]-W_n[i][j])for j in range(len(graph.graph))])for i in range(len(graph.graph))])
      lprev = l
      l = bisection(lineSearchFn,0,1,LSTol)
      if l == None:
          if lineSearchFn(1)>0:
              l = 0
          else:
              l = 1
      V_prev = V_n
      V_n = W_n+l*(E_n-W_n)
      err = computeDualGap(graph,E_n-V_n)
      if n%10 ==0:
          print(abs(np.max(err)))
      if n<2:
          W_prev2 = W_prev
          W_prev = W_n
          W_n = V_n
          
    
      else:
          lineSearchFn = lambda t: sum([sum([graph.costMatrix[i][j]((1-t)*V_n[i][j]+t*W_prev2[i][j])* (W_prev2[i][j]-V_n[i][j])for j in range(len(graph.graph))])for i in range(len(graph.graph))])
          tauprev = tau

          taumin = findTauMin(l, lprev, tauprev, taumin)
          tau = bisection(lineSearchFn,taumin,0.999,LSTol)
          if tau == None:
              if lineSearchFn(1)>0:
                  tau = taumin
              else:
                  tau = 0.999
          W_prev2 = W_prev
          W_prev = W_n
          W_n = (1-tau)*V_n + tau*W_prev2
          
      if n==maxIt:
        print("maximum iterations")
    
  return V_n, errArr
