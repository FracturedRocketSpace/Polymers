#==============================================================================
# Plot polymer
#==============================================================================


import matplotlib.pyplot as plt
import numpy as np
import config as c

def plotPolymers(polymers, endtoendDistances, weightedEndtoendSq, weightedEndtoendSqStd, minEp, popSize, weightedGyradiusSq, weightedGyradiusSqStd, lp1, fittedWeightedEndtoendSq, fittedGyradius, totalEnergy, totalEnergyCount):
    # The plots

    sortEnergy=np.sort(totalEnergy)
    lowerRange=sortEnergy[0]
    upperRange=sortEnergy[int(c.histFraction*totalEnergyCount)] #Disregard highest energy polymers, because they have much higher energy than the rest and ruin the histogram
    q1=sortEnergy[int(0.25*c.histFraction*totalEnergyCount)]    #Freedman-Diaconis method for determining optimal bin size
    q3=sortEnergy[int(0.75*c.histFraction*totalEnergyCount)]
    IQR=q3-q1
    h=2*IQR*(c.histFraction*totalEnergyCount)**(-1/3)
    b=(upperRange-lowerRange)/h

    plt.figure(4)
    n, bins, patches=plt.hist(totalEnergy, int(b), range=(lowerRange,upperRange), facecolor='green')
    plt.ylim([0,1.5*np.max(n)])
    plt.xlabel('Potential')
    plt.ylabel('Number of polymers')

    plt.figure(1)
    plt.subplot(121)
    plt.title('Polymer positions')
    plt.xlabel('x')
    plt.ylabel('y')
    for i in range(min(len(polymers),c.plotMaxPolymers)):
        r = polymers[i]
        plt.plot( r[:,0], r[:,1] )

    plt.subplot(122)
    plt.title('End-to-end distance squared vs the number of beads')
    plt.xlabel('Number of beads')
    plt.xlim(0,c.nBeads)
    plt.ylabel('End-to-end distance squared')
    x=np.linspace(1,c.nBeads, c.nBeads)
    for i in range(min(len(polymers),c.plotMaxPolymers)):
        y=endtoendDistances[i][:,0]
        plt.plot ( x, y)

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

    plt.figure(6)
    plt.xlabel('Number of beads')
    plt.xlim(0,c.nBeads)
    plt.ylabel('Gyradius squared')
    plt.errorbar(x,weightedGyradiusSq,yerr=weightedGyradiusSqStd, label='Data')
    plt.plot(x,fittedGyradius,color = "r", label='Fit')
    plt.plot(popSize, color = "g", label='Population')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([3,c.nBeads])
    plt.legend(loc='best')

    if(c.minEp):
        plt.figure(7)
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
    plt.locator_params(axis='x',nbins=4)
    plt.xlim([0,c.nBeads])
    if(c.PERM):
        plt.ylabel('Population')
        plt.plot(popSize, linewidth=2)
    else:
        plt.ylabel('Fraction alive')
        plt.plot(popSize/c.nPolymers, linewidth=2)
        plt.ylim(0,1)

    # Save the plots to file
    if(c.savePlots):
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
        plt.savefig('random.eps', bbox_inches='tight',  dpi=200)
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
        plt.figure(7)
        plt.savefig('EPmin.eps', bbox_inches='tight',  dpi=200)