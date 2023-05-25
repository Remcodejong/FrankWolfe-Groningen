import numpy as np
from IPython.core.magics.osm import line_cell_magic
from inspect import isfunction
import matplotlib.pyplot as plt
import sympy as sym
import copy
import time
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)

#Define a class for the graph network that takes in a number of vertices
class Graph():
  def __init__(self,vertexnr):
    #Make arrays for the graph, the OD pairs associated with the network,
    # A matrix to store the cost function in and a list of vertices
    self.graph = []
    self.ODPairs = []
    self.costMatrix = []
    self.V = [i for i in range(0,vertexnr)]
    self.graph = [[0 for i in range(0,len(self.V))] for j in range(0,len(self.V))]
    self.costMatrix = [[lambda x: 0 for i in range(0,len(self.V))] for j in range(0,len(self.V))]
    
    #Create a function to add an edge.
    #INPUT: - i,j: Ordered pair to mark the vertices between which the edge runs.
    #       - costFN: A lambda function that provides a cost to the new edge.
    #OUTPUT: None

  def addEdge(self,i,j,costFn):
    #If the given vertices are part of the graph, put the cost function
    #in the correct spot in the cost matrix and put a 1 in the adjacency matrix.
    if i in self.V and j in self.V:
      self.costMatrix[i][j] = costFn
      self.graph[i][j] = 1
      return
  
  #Create a function to add an OD pair to the network
  #INPUT: - i,j: Ordered pair to give the starting and ending vertex
  #OUTPUT: None

  def addODPair(self,i,j):
    #Add the OD Pair to the OD array if the vertices are in the graph.
    if i in self.V and j in self.V:
      ODArr = [i,j]
      self.ODPairs.append(ODArr)

  #Create a function to find all paths between two vertices.
  #INPUT: - start: Vertex we want to start from
  #       - end: Vertex to end on
  #       - curpath: List to store the current path we are working on
  #       - visited: List to store the visited vertices
  #OUTPUT: Prints a list of all paths between start and end
  def pathfind(self,start,end, maxLen = 0, curpath=[],visited=[]):
    #Check if the starting and ending vertex are in the graph
    if start in self.V and end in self.V:
        #Add starting vertex to visited and current path
        visited.append(start)
        curpath.append(start)

        #If we reach the end, print the full path
        if start==end:
          if len(curpath) <= maxLen or  maxLen == 0:
            print(f"{curpath}, Vertex Pair: {(curpath[0],curpath[-1])}")

        else:
          #If not, loop over all vertices, and check if they are viable
          for i in range(0,len(self.V)):
            if(i not in visited and self.graph[start][i] == 1):
              #If so, use recursion on the same graph, with a new starting vertex
              self.pathfind(i,end,maxLen, curpath,visited)
        #If we did all these vertices (or reach the end), remove the last entry
        #of the current path and the visited array to "go up a layer"
        curpath.pop()
        visited.pop()
        return
    else:
      print(f"Error, vertex {start} or vertex {end} are not in the graph.")

#Create a function to make the graph that contains marginal costs instead of the
#normal cost functions.
#
#INPUT: G: Graph we wish to make a copy of
#OUTPUT: MG: Copy of G with marginal cost functions
def MarginalGraph(G):
  #Make a copy of the original graph
  MG = copy.deepcopy(G)
  #For all edges, differentiate the cost function and put it as the corresponding
  #marginal cost in the cost matrix of MG
  for i in range(len(MG.graph)):
    for j in range(len(MG.graph)):
      t = sym.symbols("t")
      eq = lambda x: G.costMatrix[i][j](t)
      MG.costMatrix[i][j] = sym.lambdify(t,sym.diff(t*eq(t),t))
  return MG

def makeGraphSymmetric(G):
  G2 = copy.deepcopy(G)
  for i in range(0,len(G2.V)):
    for j in range(0,len(G2.V)):
        G2.graph[j][i] = G2.graph[i][j]
        G2.costMatrix[j][i] = G2.costMatrix[i][j]
  return G2