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
    
    
# Add end-to-end calculation/gyradius/other things
    
    
# Add plot of some/all polymers
plotPolymers(polymers, endtoendDistances)
