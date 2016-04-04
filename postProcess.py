#==============================================================================
# Post processing is done here
#==============================================================================

import config as c
import numpy as np
import math
from numba import jit
from scipy.optimize import curve_fit

def fitFunction(N,a):
    #the constant power is defined as 1.5 in the new version of the book
    return a * (N-1) ** 1.5

# Calculate weighted average end-to-end distance (squared) as a function of the number of beads
#TODO: Make it faster or make it compatibale with hyperdrive
#@jit( nopython=True )
def computeEndToEnd(endtoendDistances, polWeights):
    weightedEndtoendSq=np.zeros(c.nBeads)
    weightedEndtoendSqStd=np.zeros(c.nBeads)

    for z in range(c.nBeads):
        # Get data
        t1 = np.squeeze(np.asarray(endtoendDistances)[:,z])
        t2 = np.squeeze(np.asarray(polWeights)[:,z]);
        # Count non-Zero
        dataLength = len( np.flatnonzero(t2!=0) );

        weightedEndtoendSq[z]=np.average(t1, weights=t2)
        weightedEndtoendSqStd[z]=( (np.average((t1 - weightedEndtoendSq[z])**2, weights=t2))/(dataLength-1) )**(1/2)

    return weightedEndtoendSq, weightedEndtoendSqStd

 # Calculate gyradius and errors
#Interpretation: The larger the gyradius, the larger the mean squared difference with the average bead position, so the more stretched the polymer is.
#TODO: Make it faster or make it compatibale with hyperdrive
#@jit( nopython=True )
def computeGyradius(polymers, polWeights):
    # Calc gyradius
    gyradiusSq=np.zeros([len(polymers),c.nBeads])
    for w in range(len(polymers)):
        idxmax=np.max( np.argwhere(polymers[w][:,0]) )  # find highest nonzero index
        for v in range(1,idxmax+1):
            meanPosition=np.mean(polymers[w][0:v],0)             
            gyradiusSq[w,v]=1/(v)* np.trace((polymers[w][0:v]-meanPosition).dot((polymers[w][0:v]-meanPosition).T ))
    

    #Calc std        
    weightedGyradiusSq=np.zeros(c.nBeads)
    weightedGyradiusSqStd=np.zeros(c.nBeads)

    for z in range(c.nBeads):
        t3 = np.squeeze(np.asarray(polWeights)[:,z])
        dataLength = len( np.flatnonzero(t3!=0) )
        weightedGyradiusSq[z]=np.average(gyradiusSq[:,z], weights=t3)
        weightedGyradiusSqStd[z]=( (np.average((gyradiusSq[:,z] - weightedGyradiusSq[z])**2, weights=t3))/(dataLength-1) )**(1/2)

    return weightedGyradiusSq, weightedGyradiusSqStd

# Computes population at each bead
def computePopulation(polWeights):
    popSize = np.zeros(c.nBeads);
    for i in range(c.nBeads):
        temp = np.asarray(polWeights)[:,i];
        popSize[i] = len( np.flatnonzero(temp!=0) );

    return popSize


# Calculate pesistance length
#TODO: Cleanup
#TODO: Make it faster or make it compatibale with hyperdrive
#@jit( nopython=True )
def computePersistance(polymers, polWeights):
    lp1=np.zeros([len(polymers),1])
    totalWeight=np.zeros([len(polymers),1])
    for l in range(len(polymers)):
        idxmax=np.max( np.argwhere(polymers[l][:,0]) )  # find highest nonzero index
        lp1local=np.zeros([idxmax, 1 ])
        for k in range(idxmax):
            a=polymers[l][k+1,:]-polymers[l][k,:]
            lp1local[k] = np.dot(a, polymers[l][idxmax,:]-polymers[l][k,:]  ) / c.linkDistance
        if np.sum(polWeights[l][0:idxmax])>0:
            lp1[l]=np.dot( polWeights[l][0:idxmax].T , lp1local ) / np.sum(polWeights[l][0:idxmax])
            totalWeight[l]=np.sum(polWeights[l][0:idxmax])

    lp1Avg=np.dot( totalWeight.T , lp1 ) / np.sum(totalWeight);   # Take average over polymers
    # Calculate error persistence length
    nBlocks=10;
    lBlock=len(polymers)/nBlocks;
    lpBlocks=np.zeros([nBlocks,1])
    for j in range(nBlocks):
        lpBlocks[j]=np.dot( totalWeight[math.floor(j*lBlock) : math.floor( (j+1)*lBlock)].T , lp1[math.floor(j*lBlock):math.floor( (j+1)*lBlock)] ) / np.sum(totalWeight[math.floor(j*lBlock):math.floor( (j+1)*lBlock)])
    lpError=math.sqrt( (np.mean(lpBlocks*lpBlocks)- np.mean(lpBlocks)**2) / nBlocks );

    print ("Average persistence length method 1: ", lp1Avg, "Error: ", lpError)

    return lp1


def postProcess(polymers, polWeights, endtoendDistances):
    weightedEndtoendSq, weightedEndtoendSqStd = computeEndToEnd(endtoendDistances, polWeights);
    weightedGyradiusSq, weightedGyradiusSqStd = computeGyradius(polymers, polWeights);
    print("End to End done")
    print("Gyradius done")
    popt, pcov = curve_fit(fitFunction, np.arange(c.nBeads)+1 , weightedEndtoendSq);
    fittedWeightedEndtoendSq = fitFunction(np.arange(c.nBeads)+1,popt[0]);
    print("Fitting done")
    popSize = computePopulation(polWeights);
    print("Population calculated")
    lp1 = computePersistance(polymers, polWeights);
    print("Persistence length calculated")

    return weightedEndtoendSq, weightedEndtoendSqStd,  weightedGyradiusSq, weightedGyradiusSqStd, popSize, lp1, fittedWeightedEndtoendSq