#==============================================================================
# Main loop
#==============================================================================

import config as c
import numpy as np
from addPolymers import addPolymers
from plotPolymers import plotPolymers
from minimizeEp import minimizeEp
import matplotlib.pyplot as plt


# Generate the polymers
polymers, polWeights, endtoendDistances = addPolymers();

# Calculate average end-to-end distance as a function of the number of beads
#TODO: put this in another file
totalEndtoend=np.zeros(c.nBeads)
totalEndtoendSq=np.zeros(c.nBeads)
polymerEndtoend=np.zeros(len(polymers))
polymerEndtoendSq=np.zeros(len(polymers))
errorEndtoend=np.zeros(c.nBeads)
errorEndtoendSq=np.zeros(c.nBeads)

for a in range(c.nBeads):
    for z in range(len(polymers)):
        if a == c.nBeads-1:        
            totalEndtoend+=endtoendDistances[z][:,0]
            totalEndtoendSq+=endtoendDistances[z][:,1]
        polymerEndtoend[z]=endtoendDistances[z][a,0]
        polymerEndtoendSq[z]=endtoendDistances[z][a,1]
    errorEndtoend[a]=np.std(polymerEndtoend)
    errorEndtoendSq[a]=np.std(polymerEndtoendSq)

averageEndtoend=totalEndtoend/len(polymers)
averageEndtoendSq=totalEndtoendSq/len(polymers)
#Calculated errors are much larger than in the book!

# Calculate gyradius and errors
#TODO: put this in another file
gyradius=np.zeros(len(polymers))

for w in range(len(polymers)):
    meanX=np.mean(polymers[w][:,0])
    meanY=np.mean(polymers[w][:,1])
    totaldeviationSquared=0
    for v in range(c.nBeads):
        deviation=np.array([polymers[w][v,0]-meanX, polymers[w][v,1]-meanY])
        deviationSquared=np.inner(deviation,deviation)
        totaldeviationSquared+=deviationSquared
    gyradius[w]=totaldeviationSquared/c.nBeads
    averageGyradiusSquared=np.mean(gyradius)
    errorAverageGyradiusSquared=np.std(gyradius)

#Interpretation: The larger the gyradius, the larger the mean squared difference with the average bead position, so the more stretched the polymer is. 
#Of is het de bedoeling de gyradius na iedere add bead uit te rekenen, net zoals bij de end-to-end distance?

# Calculate attrition
#TODO: Put this in another file with the above 2
attrition = np.zeros(c.nBeads);
if(c.PERM is False):
    for i in range(c.nBeads):
        temp = np.asarray(polWeights)[:,i];
        attrition[i] = len( np.flatnonzero(temp!=0) )/len(polymers);

# Ep minimalisation
minEp=0;
if(c.minEp):
    minEp = minimizeEp(polymers);

# Add plot of some/all polymers
plotPolymers(polymers, endtoendDistances, averageEndtoend, errorEndtoend, averageEndtoendSq, errorEndtoendSq, minEp, attrition)
print ("Gyradius squared of each polymer:", gyradius)
print ("Average gyradius squared=", averageGyradiusSquared, "with standard deviation=", errorAverageGyradiusSquared)

# Calculate and plot pesistance length
lp1=np.zeros([c.nBeads, len(polymers) ])
lp2=np.zeros([c.nBeads, len(polymers) ])
for l in range(len(polymers)):
    for k in range(np.size(polymers[l],0)-1):
        a=polymers[l][k+1,:]-polymers[l][k,:]
        lp1[k,l] = np.dot(a, polymers[l][-1,:]-polymers[l][k,:]  ) / c.linkDistance # Definition from 2 papers
#        lp2[k,l] = np.dot(a, polymers[l][-1,:] ) / c.linkDistance # Definition from 1 paper

plt.figure(2)
plt.plot(np.mean(lp1,1))
#plt.plot(np.mean(lp2,1))
#plt.legend(["lp1", "lp2"])
print ("Average persistence length method 1: ",np.mean(lp1), "std: ", np.std(np.mean(lp1,1) ))
#print ("Average persistence length method 2: ",np.mean(lp2), "std: ", np.std(np.mean(lp2,1) ))