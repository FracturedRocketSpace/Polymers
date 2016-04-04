#==============================================================================
# Plot polymer
#==============================================================================


import matplotlib.pyplot as plt
import numpy as np
import config as c

def plotPolymers(polymers, endtoendDistances, weightedEndtoendSq, weightedEndtoendSqStd, minEp, popSize, weightedGyradiusSq, weightedGyradiusSqStd, lp1, fittedWeightedEndtoendSq):
    # Set Font and do something else; Just when output is needed.
    #font = {'family' : 'normal',
    #        'weight' : 'normal',
    #        'size'   : 25}
    #import matplotlib
    #matplotlib.rc('font', **font)
    #from matplotlib import rcParams
    #rcParams.update({'figure.autolayout': True})

    # The plots
    plt.figure(1)

    plt.subplot(231)
    plt.title('Polymer positions')
    plt.xlabel('x')
    plt.ylabel('y')
    for i in range(min(len(polymers),c.plotMaxPolymers)):
        r = polymers[i]
        plt.plot( r[:,0], r[:,1] )

    plt.subplot(232)
    plt.title('End-to-end distance squared vs the number of beads')
    plt.xlabel('Number of beads')
    plt.ylabel('End-to-end distance squared')
    x=np.linspace(1,c.nBeads, c.nBeads)
    for i in range(min(len(polymers),c.plotMaxPolymers)):
        y=endtoendDistances[i][:,0]
        plt.plot ( x, y)

    plt.subplot(233)
    plt.title('Population size')
    plt.xlabel('Number of beads')
    plt.ylabel('Number of polymers')
    plt.plot(popSize)

    plt.subplot(234)
    plt.title('Weighted average end-to-end distance squared vs number of beads')
    plt.xlabel('Number of beads')
    plt.ylabel('End-to-end distance squared')
    #plt.plot(x,weightedEndtoendSq)
    plt.errorbar(x,weightedEndtoendSq,yerr=weightedEndtoendSqStd)
    plt.plot(x,fittedWeightedEndtoendSq,color = "r")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([3,1000])

    plt.subplot(235)
    plt.title('Weighted average gyradius squared vs number of beads')
    plt.xlabel('Number of beads')
    plt.ylabel('Gyradius squared')
    plt.errorbar(x,weightedGyradiusSq,yerr=weightedGyradiusSqStd)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([3,1000])

    if(c.minEp):
        plt.subplot(236)
        plt.title('Minimal potential energy at each genetic algorithm iteration')
        plt.xlabel('Iteration')
        plt.ylabel('Potential')
        plt.plot(minEp)

    plt.figure(2)
    plt.plot(lp1)
    plt.xlabel('Polymer number')
    plt.ylabel('Average local persistence length')

    # Population size plot
    plt.figure(3)
    plt.xlabel('Number of beads')
    plt.ylabel('Number of polymers')
    plt.locator_params(axis='x',nbins=4)
    plt.plot(popSize, linewidth=2)
    #plt.savefig('popSize.eps', bbox_inches='tight',  dpi=100)

