#==============================================================================
# Main loop
#==============================================================================

import config as c
import numpy as np
from addPolymer import addPolymer
from plotPolymers import plotPolymers

# Initiate lists
polymers = [];
endtoendDistances = [];

# Initiate position vector
r = np.zeros((c.nBeads, 2));
r[1,1] = c.linkDistance;

# Initialize end to end distance vectors
endtoendDistance=np.zeros(c.nBeads)
endtoendDistance[1]=c.linkDistance

# Generate the polymers
for k in range(c.nPolymers):    
    addPolymer(np.copy(r), 2, 0, polymers, np.copy(endtoendDistance), endtoendDistances)
    
# Calculate average end-to-end distance as a function of the number of beads
totalEndtoend=np.zeros(c.nBeads)

for z in range(c.nPolymers):
    totalEndtoend+=endtoendDistances[z]
    
averageEndtoend=totalEndtoend/c.nPolymers #Errors not calculated yet

# Calculate gyradius
gyradius=np.zeros(c.nPolymers)

for w in range(c.nPolymers):
    meanX=np.mean(polymers[w][:,0])
    meanY=np.mean(polymers[w][:,1])
    for v in range(c.nBeads):
        totaldeviationSquared=0
        deviation=np.array([polymers[w][v,0]-meanX, polymers[w][v,1]-meanY])
        deviationSquared=np.inner(deviation,deviation)
        totaldeviationSquared+=deviationSquared
    gyradius[w]=totaldeviationSquared/c.nBeads

# Add plot of some/all polymers
plotPolymers(polymers, endtoendDistances, averageEndtoend)
print ("Gyradius of each polymer:", gyradius)
