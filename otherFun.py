"""Other functions"""
import numpy as np
from IPython.core.magics.osm import line_cell_magic
from inspect import isfunction
import matplotlib.pyplot as plt
import sympy as sym
import copy
import time
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
from graph import *
from FW import *
from paths import *

#Function to print all the edge flows
#INPUT: - flow: The flow to print
#OUTPUT: - print the flow on all edges
def parseFlow(flow, marginal):
  for i in range(len(flow)):
    for j in range(len(flow)):
      if flow[i][j] != 0:
        if marginal:
          print(f"The marginal flow on edge ({i},{j}) is equal to {flow[i][j]}")
        else:
          print(f"The flow on edge ({i},{j}) is equal to {flow[i][j]}")


#Function to compute the cost of a flow, using the edge cost sum outlined in the
#paper.
#INPUTS: - flow: The flow to compute the cost of.
#         - graph: The graph structure used in computing the flow.
#
#OUTPUT: sum: The total cost of the flow.
def computeCost(flow,graph):
  sum = 0
  n = len(G.graph);
  for i in range(n):
    for j in range(n):
      sum += graph.costMatrix[i][j](flow[i][j])*flow[i][j]
  return sum


#Function to create a plot of the error
def errorPlot(arr, marginal):
  x = np.linspace(1,len(arr),len(arr))
  plt.figure(figsize=(10,7.5))
  p1  = plt.semilogy(x,1/x, label = "$f(x) = 1/k$",zorder=1)
  p2 = plt.semilogy(x,arr, label = "Approximate error",zorder=2)
  if marginal:
    plt.title("Error progression of marginal iteration")
  else:
    plt.title("Error progression of normal iteration")
  plt.ylabel("Error")
  plt.xlabel("Iteration number")
  plt.legend()
  plt.show()

def POAcalc(G,P,r):
  f,_ = FrankWolfe(G,P,r)
  MG = MarginalGraph(G)
  mf,_ = FrankWolfe(MG,P,r)
  POA = computeCost(f,G)/computeCost(mf,G)
  return POA


def fullRun(G,P,r,tol,lineTol = 1e-8):
  f,errF = FrankWolfe(G,P,r,LSTol = lineTol,itTol = tol)
  MG = MarginalGraph(G)
  mf, errMf = FrankWolfe(MG,P,r, LSTol = lineTol,itTol = tol)


  print(f"\n=====================\n")
  print(f"For the given network the following flow is a user equilibrium:")
  parseFlow(f,False)
  print(f"\n=====================\n")
  print(f"For the given network the following flow is a system optimum:")
  parseFlow(mf,True)
  print(f"\n=====================\n")
  print(f"The difference between these UE flow and optimal flow is as follows:")
  parseFlow((f-mf),False)
  print(f"\n=====================\n")
  errorPlot(errF,False)
  errorPlot(errMf,True)
  Cf = computeCost(f,G)
  Cmf = computeCost(mf,G)

  POA = Cf/Cmf
  print(f"The cost of the UE flow is {Cf}, while the cost of the system optimum flow is {Cmf}.\
  The price of anarchy is therefore {POA}")
  

