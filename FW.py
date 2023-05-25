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

#Function that implements the bisection method to find the root of a given function.
#INPUT: - f: Function to find the root of. It is assumed the function maps real numbers
#             to real numbers, and not vectors to real numbers.
#       - a,b: Beginning and end points of the interval to find a root over
#       - eps: Tolerance to end iterations, based on the criterion that |f(c)|<eps
#               as decribed in the article.
#
#OUTPUT: - mid: x coordinate of the root of the function f
def bisection(f,a,b, eps):
  funA = f(a)
  funB = f(b)
  if(funA*funB <= 0):
    mid = (a+b)/2
    funMid = f(mid)
    if abs(funMid) < eps:
      return mid
    else:
      if(funA*funMid < 0):
        return bisection(f,a,mid,eps)
      else:
        return bisection(f,mid,b,eps)
  else:
    print(f"cannot guarantee root, terminating bisection method.")
    return



#A function find the cheapest path between two vertices, based on no extra flow
#INPUT: - graph: The graph structure we are working with
#       - paths: List of all paths between the desired OD pair
#
#OUTPUT: - bestIndex: Index of the best possible path to take
#        - bestPath: The best path to take from the given list of paths
def AONinitial(graph,paths):
  #Set the best cost to infinity
  best = np.inf
  #If there is only one path, return this path
  if len(paths) == 1:
    return 0,paths[0]
  
  #For each path set the new cost to zero
  for i in range(0, len(paths)):
    newCost = 0
    #Add the cost at zero traffic of each edge to the newcost
    for s in range(0,len(paths[i])-1):
       newCost += graph.costMatrix[paths[i][s]][paths[i][s+1]](0)

    #If the new cost is lower than the previous best found, set the current to the 
    #best value
    if (newCost< best):
      bestPath = paths[i]
      bestIndex = i
      best = newCost
  return bestIndex,bestPath

#Function to find the initial flow based on the zero-flow all or nothing
# assignment.
#INPUT: - graph: The graph structure we are working with
#       - paths: List of all paths between all OD-pairs.
#       - rate: List that contains all the rates for all OD-pairs.
# Note that both the paths and rates should be in the same order as the OD pairs!
#
#OUTPUT: - x: Initial flow x(0) to be used as an initial feasible flow.
def initialSetup(graph,paths,rate):
  #Make a new matrix to store all the edge flows in
  N = len(graph.graph)
  x = np.zeros((N,N))

  #For all OD-pairs, perform the zero flow all or nothing assignment
  for j in range(0,len(graph.ODPairs)):
    indexOfPath, bestpath = AONinitial(graph,paths[j])
    for i in range(len(bestpath)-1):
      #For each edge in the found best path, add the rate of the current OD-pair
      # to the flow matrix.
      x[bestpath[i]][bestpath[i+1]] += rate[j]
  return x

#Function to update the cost matrix 
#INPUT: - graph: The graph structure we are working with
#       - x: Current flow vector
#
#OUTPUT: - newCostMatrix: matrix that stores the cost of each edge for a given
#         edge flow.
def updateCostMatrix(graph,x):
  #Make an array to store the costs in
  N = len(graph.graph)
  newCostMatrix = np.zeros((N,N))
  #For each edge, check if the entry in the cost matrix is a function, and 
  # if so, add the cost of the flow on that edge to the new cost matrix value.
  for i in range(N):
    for j in range(N):
      if isfunction(graph.costMatrix[i][j]):
        newCostMatrix[i][j] += graph.costMatrix[i][j](x[i][j])
  return newCostMatrix

#Function to perform the all or nothing step for the algorithm.
#INPUT: - graph: The graph structure we are working with
#       - paths: List of all paths between all OD-pairs.
#       - costMatrix: the current costs to base the all or nothing assignment on
#       - rate: List that contains all the rates for all OD-pairs.
#
#OUTPUT: -y: matrix that corresponds to the edge flow for the all or nothing
#           flow.
def AONUpdate(graph,paths,costMatrix,rate):
  #Make a matrix to store the flows in
  y = np.zeros_like(costMatrix)

  #For all OD pairs, store the number of paths, and make a vector
  #to store the cost of the flow on each path.
  for n in range(0,len(graph.ODPairs)):
    pathNr = len(paths[n])
    pathFlowCosts = np.zeros(pathNr)

    #For all paths, add the cost of the path, corresponding to the edge
    # flows provided
    for i in range(0,pathNr):
      for j in range(0,len(paths[n][i])-1):
        pathFlowCosts[i] += costMatrix[paths[n][i][j]][paths[n][i][j+1]]
    
    #Find the best possible path out of these path flows
      bestPathIndex = np.where(pathFlowCosts == min(pathFlowCosts))[0][0]
      bestPath = paths[n][bestPathIndex]

    #Add the full rate associated with the current OD pair to all
    # of the edges in this path
    for i in range(len(bestPath)-1):
        y[bestPath[i]][bestPath[i+1]] += rate[n]
    
  return y


#Function to perform the Frank-Wolfe algorithm for a given instance of a routing
# game. Algorithm finds a UE edge flow for the instance.
#INPUT:- graph: The graph structure we are working with
#       - P: List of all paths between all OD-pairs.
#       - r: List that contains all the rates for all OD-pairs.
#
#OPTIONAL INPUTS:
# - maxIt: Maximum number of iterations to perform before terminating
#                           standard set to 1e5 iterations
# - LSTol: The tolerance to use for the bisection method, standard 1e-6
# - itTol: The tolerance to use for the iterations of the FW-algorithm, tests
# the consecutive increase/decrease. That is, the iteration terminates if
# max(abs(x(n)-x(n+1)) < itTol. Standard set to 1e-3.

def FrankWolfe(graph,P,r,maxIt=int(1e5),LSTol=1e-8,itTol=1e-8):
  #Make an array to store the approximate error values in
  errorList = np.array([])
  err = 1e6
  print(f"Starting iterations with tolerance {LSTol} for the line search\
  and tolerance {itTol} for the error tolerance.\n")
  
  #Set the iteration counter to 1 and compute the initial flow matrix.
  n = 1
  x_n = initialSetup(graph,P,r)

  #While True, 
  while 0<1:
    #Perform the update for the cost matrix and based on this new cost,
    # find a feasible direction using all or nothing assignment.
    t_n = updateCostMatrix(graph,x_n)
    s_n = AONUpdate(graph,P,t_n,r)
    
    #Define the function to perform the line search on, as described in the paper.
    lineSearchFn = lambda a:  sum([sum([graph.costMatrix[i][j](x_n[i][j] + a*(s_n[i][j]-x_n[i][j]))*(s_n[i][j]-x_n[i][j]) for j in range(0,len(graph.V))]) for i in range(0,len(graph.V))])
    
    #Try to find the step size using the bisection method. If this fails, use a
    # static step size.
    alpha = bisection(lineSearchFn,0,1,LSTol)
    if alpha == None:
      if lineSearchFn(1)>0:
        alpha = 0
      else:
        alpha = 1
    
    #Compute the new flow as a convex combination of the current flow and the feasible
    # direction.
    x_new = x_n + alpha * (s_n - x_n)

    #Compute the "error" (that is the maximal difference between iterations)
    # and put this in the error array.
    err = (np.max(np.absolute(x_n-x_new)))
    errorList = np.append(errorList,err)
  
    #If this is less than the tolerance, we are done.
    if err< itTol:
      print(f"Sucessfully found an equilibrium flow at iteration: {n}: \nApproximate error: {err}, terminating.")
      return x_new,errorList

    #Else, print the error and iteration number and reset x_n, increase iteration count.
    if n%25 == 0:
      print(f"Iteration: {n}: \nApproximate error: {err} \n")
    x_n = x_new
    n+= 1

    #Check if we ran over the maximum number of iterations. If so, terminate.
    if n > maxIt:
      print(f"Exceeded maximum iteration count. Terminating with final error {err}.")
      return x_new,  errorList

