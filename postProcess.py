#==============================================================================
# Post processing is done here
#==============================================================================

import config as c
import numpy as np
import math
from numba import jit
from scipy.optimize import curve_fit

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
        weightedEndtoendSqStd[z]= (np.average((t1 - weightedEndtoendSq[z])**2, weights=t2))**(1/2)

    return weightedEndtoendSq, weightedEndtoendSqStd

 # Calculate gyradius and errors
#Interpretation: The larger the gyradius, the larger the mean squared difference with the average bead position, so the more stretched the polymer is.
#TODO: Make it faster or make it compatibale with hyperdrive
@jit( nopython=True )
def computeGyradius(polymer, polWeight, gyradiusSq, deviation):
    idxmax= np.sum( polymer[:,1] != 0 ) +1 # find highest nonzero index
    for v in range(1,idxmax):
        meanPosition=np.array([0, 0])        
        meanPosition[0]= np.mean(polymer[0:v,0])
        meanPosition[1]= np.mean(polymer[0:v,1])
        
        deviation[:,0]= polymer[:,0] - meanPosition[0]
        deviation[:,1]= polymer[:,1] - meanPosition[1]
        
        gyradiusSq[v]=1/(v+1)  * np.sum(deviation[1:v]*deviation[1:v])
        
    return gyradiusSq.T

def computeGyradiusStd(gyradiusSq, polWeights):
#    #Calc Weighted gyradius and standard deviation
    weightedGyradiusSq=np.zeros(c.nBeads)
    weightedGyradiusSqStd=np.zeros(c.nBeads)
    
    for z in range(c.nBeads):
        t3 = np.squeeze(np.asarray(polWeights)[:,z])
        dataLength = len( np.flatnonzero(t3!=0) )
        weightedGyradiusSq[z]=np.average(gyradiusSq[:,z], weights=t3)
        weightedGyradiusSqStd[z]= (np.average((gyradiusSq[:,z] - weightedGyradiusSq[z])**2, weights=t3))**(1/2)

    return weightedGyradiusSq, weightedGyradiusSqStd

def fitEndToEnd(weightedEndtoendSq):
    # The fit function; The constant power is defined as 1.5 in the new version of the book
    def fitFunction(N,a):
        return a * (N-1) ** 1.5
    # Fit and return the fit
    popt, pcov = curve_fit(fitFunction, np.arange(c.nBeads)+1 , weightedEndtoendSq);
    print('Found coefficients for end-to-end: a =', popt[0],'; b fixed to 1.5')
    return fitFunction(np.arange(c.nBeads)+1,popt[0]);

def fitGyradius(weightedEndtoendSq):
    # The fit function; The constant power is defined as 1.5 in the new version of the book
    def fitFunction(N,a,b):
        return a * (N-1) ** b
    # Fit and return the fit
    popt, pcov = curve_fit(fitFunction, np.arange(c.nBeads)+1 , weightedEndtoendSq);
    print('Found coefficients for gyradius: a =', popt[0],'; b =', popt[1])
    return fitFunction(np.arange(c.nBeads)+1,popt[0],popt[1]);


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
    Weight=np.zeros([len(polymers),1])
    n=0
    for l in range(len(polymers)):
        idxmax=np.max( np.argwhere(polymers[l][:,0]) )  # find highest nonzero index
        if idxmax == (c.nBeads-1):                      # Only use polymers with maximum length
            lp1local=np.zeros([idxmax, 1 ])
            for k in range(idxmax):
                lref=polymers[l][k+1,:]-polymers[l][k,:]
                lp1local[k] = np.dot(lref, polymers[l][idxmax,:]-polymers[l][k,:] );
            lp1[n]=(np.mean(lp1local))
            Weight[n]=polWeights[l][idxmax]
            n+=1

    if n<10:
        print("Warning: small sample size, bad statistics")
    
    # Calculate average persistence length and standard deviation
    lp1Avg=np.average(lp1,weights=Weight);   # Take average over polymers
    lpStd=( (np.average((lp1 - lp1Avg)**2, weights=Weight)) )**(1/2)
    
    print("Persistence length: ", lp1Avg, ". Standard deviation: ", lpStd)

    return lp1
    
def crossingExists(p1,p2,p3,p4):
    #check if the lines are close to eachother    
    if(min(p1[0],p2[0]) > max(p3[0],p4[0])):
        return False;
        
    if(min(p3[0],p4[0]) > max(p1[0],p2[0])):
        return False;
        
    if(min(p1[1],p2[1]) > max(p3[1],p4[1])):
        return False;
        
    if(min(p3[1],p4[1]) > max(p1[1],p2[1])):
        return False;
    
    #regular intersection calculation    
    xDiff = [p1[0] - p2[0],p3[0] - p4[0]];
    yDiff = [p1[1] - p2[1],p3[1] - p4[1]];
    
    det = np.linalg.det([xDiff,yDiff]);
    
    #check if lines are parallel
    if(det == 0):
        return False;
        
    #determine intersection position
    d = (np.linalg.det([p1,p2]),np.linalg.det([p3,p4]));
    x = np.linalg.det([d, xDiff]) / det;
    y = np.linalg.det([d, yDiff]) / det;
    
    #check if the intersection is between the points
    if(x > min(p1[0],p2[0],p3[0],p4[0]) and x < max(p1[0],p2[0],p3[0],p4[0])):
        if(y > min(p1[1],p2[1],p3[1],p4[1]) and y < max(p1[1],p2[1],p3[1],p4[1])):
            return True;
                
    return False;

def computeAverageCrossings(polymers, polWeights):
    
    totalCrossings = 0;
    
    for p in range(len(polymers)):
        #only count completed polymers
        if(polWeights[p][-1] > 0):
            for l1 in range(1,c.nBeads-1):
                for l2 in range(l1-1):
                    if(crossingExists(polymers[p][l1],polymers[p][l1+1],polymers[p][l2],polymers[p][l2+1])):
                        totalCrossings += 1;
        print("Polymer ", p, " checked for crossings.")
    
    return totalCrossings/len(polymers)

def postProcess(polymers, polWeights, endtoendDistances):
    weightedEndtoendSq, weightedEndtoendSqStd = computeEndToEnd(endtoendDistances, polWeights);
    print("End to End done")
    
    gyradiusSq = np.zeros([len(polymers),c.nBeads])
    for polNum in range(len(polymers)):
        gyradiusSq[polNum, :] = computeGyradius(polymers[polNum], polWeights[polNum], np.zeros(c.nBeads) , np.zeros([c.nBeads,2]));
        if (polNum % 500 ==0):
            print('Gyradius polymer', polNum, 'done')
    weightedGyradiusSq, weightedGyradiusSqStd = computeGyradiusStd(gyradiusSq, polWeights)
    print("Gyradius done")
    
    fittedWeightedEndtoendSq = fitEndToEnd(weightedEndtoendSq);
    fittedGyradius = fitGyradius(weightedGyradiusSq);
    print("Fitting done")
    
    popSize = computePopulation(polWeights);
    print("Population calculated")
    
    lp1 = computePersistance(polymers, polWeights);
    print("Persistence length calculated")
    
#    averageCrossings = computeAverageCrossings(polymers, polWeights);
#    print("Average number of crossings:", averageCrossings)

    return weightedEndtoendSq, weightedEndtoendSqStd,  weightedGyradiusSq, weightedGyradiusSqStd, popSize, lp1, fittedWeightedEndtoendSq, fittedGyradius