#==============================================================================
# Plot polymer
#==============================================================================


import matplotlib.pyplot as plt
import numpy as np
import config as c

def plotPolymers(polymers, endtoendDistances, weightedEndtoend, weightedEndtoendStd, weightedEndtoendSq, weightedEndtoendSqStd, minEp, totalAttrition):
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
    plt.title('Weighted average end-to-end distance vs number of beads')
    plt.xlabel('Number of beads')
    plt.ylabel('End-to-end distance')
    plt.plot(x,weightedEndtoend)    
    #plt.errorbar(x,weightedEndtoend,yerr=weightedEndtoendStd)
    
    plt.subplot(234)
    plt.title('Weighted average squared end-to-end distance vs number of beads')
    plt.xlabel('Number of beads')
    plt.ylabel('End-to-end distance squared')
    plt.plot(x,weightedEndtoendSq)
    #plt.errorbar(x,weightedEndtoendSq,yerr=weightedEndtoendSqStd)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([3,250])
    plt.ylim([1,10000])

    plt.subplot(235)
    plt.title('Minimal potential energy at each genetetic iteration')
    plt.xlabel('Iteration')
    plt.ylabel('Potential')
    plt.plot(minEp)
    
    plt.subplot(236)
    plt.title('Attrition vs number of beads')
    plt.xlabel('Number of beads')
    plt.ylabel('Cumulative attrition number')
    plt.plot(x,totalAttrition)
    