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
    
averageEndtoend=totalEndtoend/c.nPolymers
    
# Add end-to-end calculation/gyradius/other things
    # end-to-end distance done (except for errors)  
    
# Add plot of some/all polymers
plotPolymers(polymers, endtoendDistances, averageEndtoend)
