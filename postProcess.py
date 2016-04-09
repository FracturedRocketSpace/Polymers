#==============================================================================
# Post processing is done here
#==============================================================================

import config as c
import numpy as np
from numba import jit
from scipy.optimize import curve_fit
from calculateEP import calculateEP2

# Calculate weighted average end-to-end distance (squared) as a function of the number of beads
def computeEndToEnd(endtoendDistances, polWeights):
    # Initiate variables
    weightedEndtoendSq=np.zeros(c.nBeads)
    weightedEndtoendSqStd=np.zeros(c.nBeads)
    
    # Loop over all possible polymer lengths
    for z in range(c.nBeads):
        # Get data
        t1 = np.squeeze(np.asarray(endtoendDistances)[:,z])
        t2 = np.squeeze(np.asarray(polWeights)[:,z]);
        # Count non-Zero
        dataLength = len( np.flatnonzero(t2!=0) );
        # Calculate mean and standard deveation of mean
        weightedEndtoendSq[z]=np.average(t1, weights=t2)
        weightedEndtoendSqStd[z]=( (np.average((t1 - weightedEndtoendSq[z])**2, weights=t2)) / (dataLength))**(1/2)

    return weightedEndtoendSq, weightedEndtoendSqStd

 # Calculate gyradius and errors
@jit( nopython=False )
def computeGyradius(polymer, polWeight, gyradiusSq, deviation):
    # find highest nonzero index
    idxmax= np.sum( polymer[:,1] != 0 ) +1
    
    # Determine gyradius for all polymer lengths
    for v in range(1,idxmax):
        # Get average coordinate in polymer
        meanPosition=np.array([0, 0])
        meanPosition[0]= np.mean(polymer[0:v,0])
        meanPosition[1]= np.mean(polymer[0:v,1])
        # Calculate deviation from mean position of all beads
        deviation[:,0]= polymer[:,0] - meanPosition[0]
        deviation[:,1]= polymer[:,1] - meanPosition[1]
        # Compute gyradius
        gyradiusSq[v]=1/(v+1)  * np.sum(deviation[1:v]*deviation[1:v])

    return gyradiusSq.T

def computeGyradiusStd(gyradiusSq, polWeights):  
    # Calc Weighted gyradius and standard deviation

    # Initiate variables
    weightedGyradiusSq=np.zeros(c.nBeads)
    weightedGyradiusSqStd=np.zeros(c.nBeads)
    
    # Loop over all possible polymer lengths
    for z in range(c.nBeads):
        # Get weight for bead number z in all polymers
        w = np.squeeze(np.asarray(polWeights)[:,z])
        # Count nonzero
        dataLength = len( np.flatnonzero(w!=0) )
        # Calculate mean and standard deveation of mean
        weightedGyradiusSq[z]=np.average(gyradiusSq[:,z], weights=w)
        weightedGyradiusSqStd[z]=( (np.average((gyradiusSq[:,z] - weightedGyradiusSq[z])**2, weights=w)) / (dataLength))**(1/2)

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
    # Initiate variables
    popSize = np.zeros(c.nBeads);
    
    # Sum nonzero weights to get population size
    for i in range(c.nBeads):
        temp = np.asarray(polWeights)[:,i];
        popSize[i] = len( np.flatnonzero(temp!=0) );

    return popSize

# Calculate pesistance length
def computePersistance(polymers, polWeights):
    # Initiate variables
    lp1=np.zeros([len(polymers),1])
    Weight=np.zeros([len(polymers),1])
    n=0
    
    # Loop over all polymers
    for l in range(len(polymers)):
        # find highest nonzero index
        idxmax=np.max( np.argwhere(polymers[l][:,0]) )
        # Only use polymers with maximum length
        if idxmax == (c.nBeads-1):
            # Initiate variable
            lp1local=np.zeros([idxmax, 1 ])
            # Determine local persistence length for all beads
            for k in range(idxmax):
                lref=polymers[l][k+1,:]-polymers[l][k,:]
                lp1local[k] = np.dot(lref, polymers[l][idxmax,:]-polymers[l][k,:] );
            
            # Calculate average of local persistence length and save weight
            lp1[n]=(np.mean(lp1local))
            Weight[n]=polWeights[l][idxmax]
            
            n+=1
            
    # Warn if small sample size
    if n<10:
        print("Warning: small sample size, bad statistics")

    # Calculate average persistence length and standard deviation
    lp1Avg=np.average(lp1,weights=Weight);
    lpStd=( (np.average((lp1 - lp1Avg)**2, weights=Weight)) / n )**(1/2)
    
    print("Persistence length: ", lp1Avg, ". Standard deviation: ", lpStd)

    return lp1

#numpy determinant does not work with numba for nonpython = True
@jit( nopython=True )
def det(a, b):
    return a[0] * b[1] - a[1] * b[0];
    
#determine if a crossing exists between two links
@jit( nopython=True )
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

    determinant = det(xDiff,yDiff);

    #check if lines are parallel
    if(determinant == 0):
        return False;

    #determine intersection position
    d = [det(p1,p2),det(p3,p4)];
    x = det(d, xDiff) / determinant;
    y = det(d, yDiff) / determinant;

    #check if the intersection is between the points
    if(x > min(p1[0],p2[0],p3[0],p4[0]) and x < max(p1[0],p2[0],p3[0],p4[0])):
        if(y > min(p1[1],p2[1],p3[1],p4[1]) and y < max(p1[1],p2[1],p3[1],p4[1])):
            return True;

    return False;

@jit( nopython=False )
def computeAverageCrossings(polymers, polWeights):
    # Initiate variables
    totalCrossings = 0;
    aliveCount = 0
    
    # Detect and note crossings for all polymers
    for p in range(len(polymers)):
        #only count completed polymers
        if(polWeights[p][-1] > 0):
            aliveCount += 1;
            for l1 in range(1,c.nBeads-1):
                for l2 in range(l1-1):
                    if(crossingExists(polymers[p][l1],polymers[p][l1+1],polymers[p][l2],polymers[p][l2+1])):
                        totalCrossings += 1;
        if (p % 500 == 0):
            print("Polymer ", p, " checked for crossings.");

    return totalCrossings/aliveCount;

def computeSortedEnergy(polymers):
    #Check energy distribution

    # Initiate variable
    totalEnergy = [];

    for a in range(len(polymers)):
        # find highest nonzero index
        idxmax=np.max( np.argwhere(polymers[a][:,0]) )
        # Only check full length polymers
        if(idxmax==c.nBeads-1):
            totalEnergy.append(calculateEP2(polymers[a]))

    # Sort energy and retunr
    return np.sort(totalEnergy)

def postProcess(polymers, polWeights, endtoendDistances):
    # Runs all individual functions  
    
    # Sorted energy
    sortedEnergy = computeSortedEnergy(polymers);
    
    # Ent to end distance
    weightedEndtoendSq, weightedEndtoendSqStd = computeEndToEnd(endtoendDistances, polWeights);
    print("End to End done")
    
    # Gyradius
    gyradiusSq = np.zeros([len(polymers),c.nBeads])
    # Loop over all polymers
    for polNum in range(len(polymers)):
        gyradiusSq[polNum, :] = computeGyradius(polymers[polNum], polWeights[polNum], np.zeros(c.nBeads) , np.zeros([c.nBeads,2]));
        # Print progress        
        if (polNum % 500 ==0):
            print('Gyradius polymer', polNum, 'done')
    weightedGyradiusSq, weightedGyradiusSqStd = computeGyradiusStd(gyradiusSq, polWeights)
    print("Gyradius done")
    
    # Fit end to end and gyradius
    fittedWeightedEndtoendSq = fitEndToEnd(weightedEndtoendSq);
    fittedGyradius = fitGyradius(weightedGyradiusSq);
    print("Fitting done")
    
    # Population size
    popSize = computePopulation(polWeights);
    print("Population calculated")
    
    # Persistence length
    lp1 = computePersistance(polymers, polWeights);
    print("Persistence length calculated")
    
    # Crossing detector
    averageCrossings = computeAverageCrossings(polymers, polWeights);
    print("Average number of crossings:", averageCrossings)

    return weightedEndtoendSq, weightedEndtoendSqStd,  weightedGyradiusSq, weightedGyradiusSqStd, popSize, lp1, fittedWeightedEndtoendSq, fittedGyradius, sortedEnergy