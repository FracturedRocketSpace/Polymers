#==============================================================================
# Plot polymer
#==============================================================================


import matplotlib.pyplot as plt
import numpy as np
import config as c

def plotPolymers(polymers, endtoendDistances, averageEndtoend, averageEndtoendSq, minEp):
    plt.figure(1)
    plt.hold(True)
    plt.title('Polymer positions')
    plt.xlabel('x')
    plt.ylabel('y')
    for i in range(len(polymers)):
        r = polymers[i]
        plt.plot( r[:,0], r[:,1] )

    plt.figure(2)
    plt.hold(True)
    plt.title('End-to-end distance as a function of the number of beads')
    plt.xlabel('Number of beads')
    plt.ylabel('End-to-end distance')
    x=np.linspace(1,c.nBeads, c.nBeads)
    for i in range(len(endtoendDistances)):
        y=endtoendDistances[i][:,0]
        plt.plot ( x, y)

    plt.figure(3)
    plt.title('Average end-to-end distance as a function of the number of beads')
    plt.xlabel('Number of beads')
    plt.ylabel('End-to-end distance')
    plt.plot(x,averageEndtoend)
    
    plt.figure(4)
    plt.title('Average squared end-to-end distance vs number of beads')
    plt.xlabel('Number of beads')
    plt.ylabel('End-to-end distance squared')
    plt.plot(x,averageEndtoendSq)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([3,1000])

    plt.figure(5)
    plt.title('Minimal potential energy at each genetetic iteration')
    plt.xlabel('Iteration')
    plt.ylabel('Potential')
    plt.plot(minEp)