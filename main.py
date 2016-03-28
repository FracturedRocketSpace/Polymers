#==============================================================================
# Main loop
#==============================================================================

import config as c
import numpy as np
from addPolymer import addPolymer
from plotPolymers import plotPolymers
from minimizeEp import minimizeEp
import matplotlib.pyplot as plt


# Initiate lists
polymers = [];
endtoendDistances = [];
weights=[];
attrition=[];

# Initiate position vector
r = np.zeros((c.nBeads, 2));
r[1,1] = c.linkDistance;
polWeight=np.zeros((c.nBeads,1))

# Initialize end to end distance (and squared) vectors
endtoendDistance=np.zeros([c.nBeads,2])
endtoendDistance[1,:]=c.linkDistance

#Initialize polWeight
polWeight[0,0]=1
polWeight[1,0]=1

# Generate the polymers
for k in range(c.nPolymers):
    addPolymer(np.copy(r), 2, polymers, np.copy(endtoendDistance), endtoendDistances, np.copy(polWeight), weights, attrition)
    print("Polymer", k+1, "has been constructed")

# Calculate average end-to-end distance as a function of the number of beads
sigmaWeightedEndtoend=np.zeros(c.nBeads)
sigmaWeightedEndtoendSq=np.zeros(c.nBeads)
sigmaWeights=np.zeros(c.nBeads)
weightedEndtoendVarNumerator=np.zeros(c.nBeads)
weightedEndtoendSqVarNumerator=np.zeros(c.nBeads)

for z in range(len(polymers)):
    sigmaWeightedEndtoend+=weights[z][:,0]*endtoendDistances[z][:,0]
    sigmaWeightedEndtoendSq+=weights[z][:,0]*endtoendDistances[z][:,1]
    sigmaWeights+=weights[z][:,0]

weightedEndtoend=sigmaWeightedEndtoend/sigmaWeights
weightedEndtoendSq=sigmaWeightedEndtoendSq/sigmaWeights
#Fluctuations for large number of beads also visible in our simulation now (compare with Thijssen fig. 10.3)         

for z in range(len(polymers)):
    weightedEndtoendVarNumerator+=weights[z][:,0]*(endtoendDistances[z][:,0]-weightedEndtoend)**2
    weightedEndtoendSqVarNumerator+=weights[z][:,0]*(endtoendDistances[z][:,1]-weightedEndtoendSq)**2

weightedEndtoendVar=weightedEndtoendVarNumerator/sigmaWeights
weightedEndtoendSqVar=weightedEndtoendSqVarNumerator/sigmaWeights

weightedEndtoendStd=weightedEndtoendVar**0.5
weightedEndtoendSqStd=weightedEndtoendSq**0.5
#Not sure if this method for weighted variance is correct. Errors still do not correspond with Thijssen fig 10.3.

# Calculate attrition
totalAttrition=np.zeros(c.nBeads)
for w in range(len(attrition)):
    totalAttrition+=(attrition[w][:,0]==0)

# Calculate gyradius and errors -> Needs to be improved
gyradius=np.zeros(c.nPolymers)

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

# Ep minimalisation
minEp=0;
if(c.minEp):
    minEp = minimizeEp(polymers);

# Add plot of some/all polymers
plotPolymers(polymers, endtoendDistances, weightedEndtoend, weightedEndtoendStd, weightedEndtoendSq, weightedEndtoendSqStd, minEp, totalAttrition)
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