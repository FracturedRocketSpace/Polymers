#==============================================================================
# Main loop
#==============================================================================

import config as c
import numpy as np
import matplotlib.pyplot as plt
from addPolymers import addPolymers
from plotPolymers import plotPolymers
from minimizeEp import minimizeEp
from postProcess import postProcess
from calculateEP import calculateEP3

# Generate the polymers
polymers, polWeights, endtoendDistances = addPolymers();

#Check energy distribution
totalEnergy=np.zeros(len(polymers))

for a in range(len(polymers)):
    idxmax=np.max( np.argwhere(polymers[a][:,0]) )            # find highest nonzero index
    totalEnergy[a]=calculateEP3(polymers[a][0:idxmax],idxmax)

sortEnergy=np.sort(totalEnergy)
lowerRange=sortEnergy[0]
upperRange=sortEnergy[int(c.histFraction*len(polymers))] #Disregard highest energy polymers, because they have much higher energy than the rest and ruin the histogram
q1=sortEnergy[int(0.25*c.histFraction*len(polymers))]    #Freedman-Diaconis method for determining optimal bin size
q3=sortEnergy[int(0.75*c.histFraction*len(polymers))]
IQR=q3-q1
h=2*IQR*(c.histFraction*len(polymers))**(-1/3)
b=(upperRange-lowerRange)/h

plt.figure(4)
n, bins, patches=plt.hist(totalEnergy, int(b), range=(lowerRange,upperRange), facecolor='green')
plt.ylim([0,1.5*np.max(n)])
plt.xlabel('Total energy')
plt.ylabel('Number of polymers')
plt.title('Total energy distribution')
plt.show()

#Post processing
weightedEndtoendSq, weightedEndtoendSqStd,  weightedGyradiusSq, weightedGyradiusSqStd, popSize, lp1, fittedWeightedEndtoendSq = postProcess(polymers, polWeights, endtoendDistances);

# Ep minimalisation
minEp=0;
if(c.minEp):
    minEp = minimizeEp(polymers);

# Add plot of some/all polymers
plotPolymers(polymers, endtoendDistances, weightedEndtoendSq, weightedEndtoendSqStd, minEp, popSize, weightedGyradiusSq, weightedGyradiusSqStd, lp1, fittedWeightedEndtoendSq)
