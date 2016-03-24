#==============================================================================
# Main loop
#==============================================================================

import config as c
import numpy as np
from addPolymers import addPolymers
from plotPolymers import plotPolymers
from minimizeEp import minimizeEp

# Generate the polymers
polymers, polWeights, endtoendDistances = addPolymers();

# Calculate average end-to-end distance as a function of the number of beads
#TODO: put this in another file
totalEndtoend=np.zeros(c.nBeads)
totalEndtoendSq=np.zeros(c.nBeads)
polymerEndtoend=np.zeros(c.nPolymers)
polymerEndtoendSq=np.zeros(c.nPolymers)
errorEndtoend=np.zeros(c.nBeads)
errorEndtoendSq=np.zeros(c.nBeads)

for a in range(c.nBeads):
    for z in range(c.nPolymers):
        if a == c.nBeads-1:        
            totalEndtoend+=endtoendDistances[z][:,0]
            totalEndtoendSq+=endtoendDistances[z][:,1]
        polymerEndtoend[z]=endtoendDistances[z][a,0]
        polymerEndtoendSq[z]=endtoendDistances[z][a,1]
    errorEndtoend[a]=np.std(polymerEndtoend)
    errorEndtoendSq[a]=np.std(polymerEndtoendSq)

averageEndtoend=totalEndtoend/c.nPolymers
averageEndtoendSq=totalEndtoendSq/c.nPolymers
#Calculated errors are much larger than in the book!

# Calculate gyradius and errors
#TODO: put this in another file
gyradius=np.zeros(c.nPolymers)

for w in range(c.nPolymers):
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
plotPolymers(polymers, endtoendDistances, averageEndtoend, errorEndtoend, averageEndtoendSq, errorEndtoendSq, minEp)
print ("Gyradius squared of each polymer:", gyradius)
print ("Average gyradius squared=", averageGyradiusSquared, "with standard deviation=", errorAverageGyradiusSquared)