#==============================================================================
# Main loop
#==============================================================================

import config as c
import numpy as np
from addPolymer import addPolymer

# Initiate list
polymers = [];

# Initiate postition vector
r = np.zeros((c.nBeads, 2));
r[1,1] = c.linkDistance;

# Generate the polymers
for k in range(c.nPolymers):    
    addPolymer(np.copy(r), 2, 0, polymers)
    
    
# Add end-to-end calculation/gyradius/other things
    
# Add plot of some/all polymers
