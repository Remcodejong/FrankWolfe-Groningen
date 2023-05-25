#Braess's paradox
G1 = Graph(4)
G1.addODPair(0,3)
G1.addEdge(0,1,lambda x: 10*x)
G1.addEdge(1,3, lambda x: 50+x)
G1.addEdge(0,2,lambda x: 10*x)
G1.addEdge(2,3,lambda x: 50+x)

paths = [[[0,1,3],[0,2,3]]]
r = [6]
f, _ = FrankWolfe(G1,paths,r)


r = [[6]]
G2 = Graph(4)
G2.addODPair(0,3)
G2.addEdge(0,1,lambda x: 10*x)
G2.addEdge(0,2,lambda x: 50+x)
G2.addEdge(2,3,lambda x: 10*x)
G2.addEdge(1,3,lambda x: 50+x)
G2.addEdge(1,2,lambda x: 10+x)
paths2 = [[[0,1,2,3],[0,1,3],[0,2,3]]]
fNew, _ = FrankWolfe(G2,paths2,r)

print(f"\n=====================\n")
print(f"For the network without the extra edge the following flow is a user equilibrium:")
parseFlow(f,False)

print(f"\n=====================\n")
print(f"For the network with the extra edge the following flow is a user equilibrium:")
parseFlow(fNew,False)
print(f"\n=====================\n")


#Pigou example
G1 = Graph(4)
G1.addODPair(0,3)
G1.addEdge(0,1,lambda x: 10*x)
G1.addEdge(1,3, lambda x: 50+x)
G1.addEdge(0,2,lambda x: 10*x)
G1.addEdge(2,3,lambda x: 50+x)

paths = [[[0,1,3],[0,2,3]]]
r = [[1]]

pmax = 100
POAs = np.array([])

for pow in range(1,pmax):
  G = Graph(4)
  G.addODPair(0,3)
  G.addEdge(0,1,lambda x: 0.5)
  G.addEdge(1,3, lambda x: 0.5)
  G.addEdge(0,2,lambda x: 0.5*x**pow)
  G.addEdge(2,3,lambda x: 0.5*x**pow)
  POA = POAcalc(G,paths,r)
  POAs = np.append(POAs,POA)

    
x = np.arange(1,pmax,1)
plt.figure(figsize=(10,7.5))
p2 = plt.plot(x,POAs, label = "Price of Anarchy")
plt.title("Price of anarchy for different values of $p$ in the nonlinear Pigou example")
plt.ylabel("Price of anarchy")
plt.xlabel("$p$")
plt.legend()
plt.show()
plt.show()

plt.figure(figsize=(10,7.5))
p2 = plt.plot(x,np.absolute(POAs-1/(1-x*(x+1)**(-(x+1)/x))),label="Difference POA")
plt.title("Absolute difference between computed POA and expected POA")
plt.ylabel("Difference")
plt.xlabel("$p$")
plt.legend()
plt.show()
plt.show()


#Small example

#Define the rate value for each OD Pair
r = [1,2]

#Define the graph network
G = Graph(11)
G.addODPair(0,9)
G.addODPair(3,10)
G.addEdge(0,1,lambda x: 6*x+4)
G.addEdge(0,2,lambda x: 2*x+1)
G.addEdge(0,6,lambda x: x+1)
G.addEdge(1,3,lambda x: 5*x+2)
G.addEdge(1,7,lambda x: 3*x+1)
G.addEdge(2,3,lambda x: 2*x+6)
G.addEdge(2,5,lambda x: 2*x+3)
G.addEdge(3,7,lambda x: 3*x+1)
G.addEdge(3,6,lambda x: 4*x+1)
G.addEdge(4,2,lambda x: x+4)
G.addEdge(5,8,lambda x: 2*x+3)
G.addEdge(6,8,lambda x: x+2)
G.addEdge(7,8,lambda x: 2*x+2)
G.addEdge(7,10,lambda x: 8*x+6)
G.addEdge(8,4,lambda x: x+4)
G.addEdge(8,9,lambda x: 12*x+4)
G.addEdge(8,10,lambda x: 0.5*x**2+1)
G.addEdge(10,9,lambda x: 3*x+1)


paths = [[[0, 1, 3, 6, 8, 9], [0, 1, 3, 6, 8, 10, 9], [0, 1, 3, 7, 8, 9],
          [0, 1, 3, 7, 8, 10, 9], [0, 1, 3, 7, 10, 9], [0, 1, 7, 8, 9], 
           [0, 1, 7, 8, 10, 9], [0, 1, 7, 10, 9], [0, 2, 3, 6, 8, 9],[0, 2, 3, 6, 8, 10, 9],
           [0, 2, 3, 7, 8, 9], [0, 2, 3, 7, 8, 10, 9], [0, 2, 3, 7, 10, 9], [0, 2, 5, 8, 9],\
           [0, 2, 5, 8, 10, 9], [0, 6, 8, 4, 2, 3, 7, 10, 9], [0, 6, 8, 9], [0, 6, 8, 10, 9]],
[[3, 6, 8, 10], [3, 7, 8, 10], [3, 7, 10]]]

fullRun(G,paths,r,1e-2)
