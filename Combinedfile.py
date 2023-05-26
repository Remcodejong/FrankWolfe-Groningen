# -*- coding: utf-8 -*-
"""
Created on Thu May 25 12:15:22 2023

@author: remco
"""
import numpy as np
from inspect import isfunction
import matplotlib.pyplot as plt
import sympy as sym
import copy
import time
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
from graph import *
from otherFun import *
from FW import *
from paths import *


#Set up the graph
G = Graph(29)
mu,sigma = 1,0.05
np.random.seed(5054572)

C1 =np.random.normal(1,0.05,size=(len(G.V),len(G.V)))
C2 =np.random.normal(1,0.05,size=(len(G.V),len(G.V)))
def goodLane(x,i,j):
  return 3*C1[i][j]*(1+0.15*(x/(7*C2[i][j]))**4)

def midLane(x,i,j):
  return 5*C1[i][j]*(1+0.15*(x/(6.5*C2[i][j]))**4)

def badLane(x,i,j):
  return 10*C1[i][j]*(1+0.15*(x/(6*C2[i][j]))**4)

G.addEdge(0,1,lambda x: goodLane(x,0,1))

G.addEdge(1,10,lambda x: midLane(x,1,10))
G.addEdge(1,2,lambda x: goodLane(x,1,2))

G.addEdge(2,3,lambda x: goodLane(x,2,3))
G.addEdge(2,7,lambda x: midLane(x,2,7))

G.addEdge(3,8,lambda x: badLane(x,3,8))
G.addEdge(3,4,lambda x: goodLane(x,3,4))

G.addEdge(4,5,lambda x: goodLane(x,4,5))
G.addEdge(4,9, lambda x: midLane(x,4,9))
G.addEdge(4,14, lambda x: midLane(x,4,14))

G.addEdge(5,6, lambda x: goodLane(x,5,6))
G.addEdge(5,22,lambda x: goodLane(x,5,22))

G.addEdge(7,11,lambda x: midLane(x,7,11))
G.addEdge(7,8,lambda x: badLane(x,7,8))
G.addEdge(7,12,lambda x: badLane(x,7,12))
G.addEdge(7,16,lambda x: midLane(x,7,16))


G.addEdge(8,9,lambda x: badLane(x,8,9))
G.addEdge(8,12, lambda x: badLane(x,8,12))
G.addEdge(8,13,lambda x: badLane(x,8,13))

G.addEdge(9,13, lambda x: badLane(x,9,13))

G.addEdge(10,11,lambda x: midLane(x,10,11))
G.addEdge(10,15,lambda x: midLane(x,10,15))

G.addEdge(11,15,lambda x: midLane(x,11,15))

G.addEdge(12,13,lambda x: badLane(x,12,13))
G.addEdge(12,16,lambda x: badLane(x,12,16))
G.addEdge(12,17,lambda x: badLane(x,12,17))

G.addEdge(13,14,lambda x: midLane(x,13,14))
G.addEdge(13,17,lambda x: badLane(x,13,17))
G.addEdge(13,21,lambda x: badLane(x,13,21))

G.addEdge(14,22,lambda x: goodLane(x,14,22))

G.addEdge(15,16,lambda x: midLane(x,15,16))
G.addEdge(15,18,lambda x: midLane(x,15,18))

G.addEdge(16,19,lambda x: goodLane(x,16,19))

G.addEdge(17,20,lambda x: badLane(x,17,20))

G.addEdge(18,19, lambda x: goodLane(x,18,19))

G.addEdge(19,20,lambda x: goodLane(x,19,20))
G.addEdge(19,23,lambda x: midLane(x,19,23))

G.addEdge(20,21,lambda x: badLane(x,20,21))
G.addEdge(20,24,lambda x: goodLane(x,20,24))
G.addEdge(20,25,lambda x: badLane(x,20,25))

G.addEdge(21,25,lambda x: badLane(x,21,25))
G.addEdge(21,22,lambda x: midLane(x,21,22))

G.addEdge(23,26,lambda x: midLane(x,23,26))
G.addEdge(23,24,lambda x: midLane(x,23,24))

G.addEdge(24,27,lambda x: goodLane(x,24,27))

G.addEdge(25,28,lambda x: midLane(x,25,28))

G.addEdge(27,28,lambda x: goodLane(x,27,28))

G = makeGraphSymmetric(G)

originNodes = [0,4,6,11,18,28]
destNodes = [2,8,12,13,20,22,26]
for u in originNodes:
  for v in destNodes:
    G.addODPair(u,v) 


# for pair in G.ODPairs:
#   G.pathfind(pair[0],pair[1],8)
#   print(f"\n")

TrafficRates = [12,16.17,14,6,40,30]

r = []
for i in range(0,len(TrafficRates)):
  for j in range(1,8):
    r.append(TrafficRates[i]/7)
    
start1 = time.time()
F1, _ = FrankWolfe(G,paths,r,itTol=1e-2,LSTol = 1e-2)
time1 = time.time()-start1

start2 = time.time()
F2, _ = FrankWolfe(G,paths,r,itTol=1e-3,LSTol = 1e-2)
time2 = time.time()-start2
print(f"The maximal difference is {np.max(np.abs(F1-F2))} and the time difference is {time2-time1}")