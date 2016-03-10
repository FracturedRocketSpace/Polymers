#==============================================================================
# Plot polymer
#==============================================================================


import matplotlib.pyplot as plt
import numpy as np
import config as c

def plotPolymers(polymers, endtoendDistances, averageEndtoend, errorEndtoend, averageEndtoendSq, errorEndtoendSq, minEp):
    plt.figure(1)

    plt.subplot(231)
    plt.title('Polymer positions')
    plt.xlabel('x')
    plt.ylabel('y')
    for i in range(len(polymers)):
        r = polymers[i]
        plt.plot( r[:,0], r[:,1] )

    plt.subplot(232)
    plt.title('End-to-end distance as a function of the number of beads')
    plt.xlabel('Number of beads')
    plt.ylabel('End-to-end distance')
    x=np.linspace(1,c.nBeads, c.nBeads)
    for i in range(len(endtoendDistances)):
        y=endtoendDistances[i][:,0]
        plt.plot ( x, y)

    plt.subplot(233)
    plt.title('Average end-to-end distance as a function of the number of beads')
    plt.xlabel('Number of beads')
    plt.ylabel('End-to-end distance')
    plt.errorbar(x,averageEndtoend,yerr=errorEndtoend)
    
    plt.subplot(234)
    plt.title('Average squared end-to-end distance vs number of beads')
    plt.xlabel('Number of beads')
    plt.ylabel('End-to-end distance squared')
    plt.errorbar(x,averageEndtoendSq,yerr=errorEndtoendSq)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([3,1000])

    plt.subplot(235)
    plt.title('Minimal potential energy at each genetetic iteration')
    plt.xlabel('Iteration')
    plt.ylabel('Potential')
    plt.plot(minEp)