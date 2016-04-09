#==============================================================================
# Plot polymer
#==============================================================================


import matplotlib.pyplot as plt
import numpy as np
import config as c

def plotPolymerPositions(polymers,polWeights):
    # Polymers
    plt.figure(1)
    plt.xlabel('x')
    plt.ylabel('y')
    polPlotted = 0;
    #plot polymers starting from the center for a better representation when using PERM
    i=int(len(polymers)/2);
    while(polPlotted < min(len(polymers),c.plotMaxPolymers)):
        if(polWeights[i][-1] > 0):
            r = polymers[i];
            plt.plot( r[:,0], r[:,1] )
            polPlotted += 1;
        i += 1;

def plotPersistenceLength(lp1):

    # Persistence length
    plt.figure(2)
    plt.plot(lp1)
    plt.xlabel('Polymer number')
    plt.ylabel('Average local persistence length')

def plotPopulationSize(popSize):

    # Population size plot
    plt.figure(3)
    plt.xlabel('Number of beads')
    plt.locator_params(axis='x',nbins=4)
    plt.xlim([0,c.nBeads])
    if(c.PERM):
        plt.ylabel('Population')
        plt.plot(popSize, linewidth=2)
    else:
        plt.ylabel('Fraction alive')
        plt.plot(popSize/c.nPolymers, linewidth=2)
        plt.ylim(0,1)

def plotEnergyHistogram(sortedEnergy):

    # Energy distribution histogram
    lowerRange=sortedEnergy[0]
    upperRange=sortedEnergy[int(c.histFraction*len(sortedEnergy) )] # Disregard highest energy polymers, because they have much higher energy than the rest and ruin the histogram
    # Freedman-Diaconis method for determining optimal bin size
    q1=sortedEnergy[int(0.25*c.histFraction*len(sortedEnergy))]
    q3=sortedEnergy[int(0.75*c.histFraction*len(sortedEnergy))]
    IQR=q3-q1
    h=2*IQR*(c.histFraction*len(sortedEnergy))**(-1/3)
    b=(upperRange-lowerRange)/h

    plt.figure(4)
    n, bins, patches=plt.hist(sortedEnergy, int(b), range=(lowerRange,upperRange), facecolor='green')
    plt.ylim([0,1.5*np.max(n)])
    plt.xlabel('Potential')
    plt.ylabel('Number of polymers')
    plt.title('Total energy distribution')

def plotEndtoendSq(weightedEndtoendSq,weightedEndtoendSqStd,fittedWeightedEndtoendSq,popSize):

    x = np.linspace(1,c.nBeads, c.nBeads);

    # Plot the end-to-end distance squared
    plt.figure(5)
    plt.xlabel('Number of beads')
    plt.xlim(0,c.nBeads)
    plt.ylabel('End-to-end distance squared')
    plt.errorbar(x,weightedEndtoendSq,yerr=weightedEndtoendSqStd, label='Data')
    plt.plot(x,fittedWeightedEndtoendSq,color = "r", label='Fit')
    plt.plot(popSize, color = "g", label='Population')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([3,c.nBeads])
    plt.legend(loc='best')

def plotGyradiusSq(weightedGyradiusSq,weightedGyradiusSqStd,fittedGyradius,popSize):

    x = np.linspace(1,c.nBeads, c.nBeads);

    # Plot the gyradius
    plt.figure(6)
    plt.xlabel('Number of beads')
    plt.xlim(0,c.nBeads)
    plt.ylabel('Ensemble average $R_g^2$')
    plt.errorbar(x,weightedGyradiusSq,yerr=weightedGyradiusSqStd, label='Data')
    plt.plot(x,fittedGyradius,color = "r", label='Fit')
    plt.plot(popSize, color = "g", label='Population')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([3,c.nBeads])
    plt.legend(loc='best')

def plotMinEp(minEp):

    # Plot the minimal potential energy at every genetic algorithm iteration
    plt.figure(7)
    plt.title('Minimal potential energy at each genetic algorithm iteration')
    plt.xlabel('Iteration')
    plt.ylabel('Potential')
    plt.plot(minEp)

def plotEndtoendSqPolymers(polymers,endtoendDistances):

    x = np.linspace(1,c.nBeads, c.nBeads);    
    
    # Plot End-to-end distances of polymers
    plt.figure(8)
    plt.title('End-to-end distance squared vs the number of beads')
    plt.xlabel('Number of beads')
    plt.xlim(0,c.nBeads)
    plt.ylabel('End-to-end distance squared')
    for i in range(min(len(polymers),c.plotMaxPolymers)):
        y=endtoendDistances[i][:,0]
        plt.plot ( x, y)

def savePlots():
    # Save the plots to file
    # Set font and layout.
    font = {'family' : 'normal',
           'weight' : 'normal',
            'size'   : 25}
    import matplotlib
    matplotlib.rc('font', **font)
    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})
    # Save to file
    plt.figure(1)
    plt.savefig('polymers.eps', bbox_inches='tight',  dpi=200)
    plt.figure(2)
    plt.savefig('persistence.eps', bbox_inches='tight',  dpi=200)
    plt.figure(3)
    plt.savefig('popSize.eps', bbox_inches='tight',  dpi=200)
    plt.figure(4)
    plt.savefig('histogram.eps', bbox_inches='tight',  dpi=200)
    plt.figure(5)
    plt.savefig('e2e.eps', bbox_inches='tight',  dpi=200)
    plt.figure(6)
    plt.savefig('Gyradius.eps', bbox_inches='tight',  dpi=200)
    if(c.minEp):
        plt.figure(7)
        plt.savefig('EPmin.eps', bbox_inches='tight',  dpi=200)
    plt.figure(8)
    plt.savefig('e2ePolymers.eps', bbox_inches='tight',  dpi=200)