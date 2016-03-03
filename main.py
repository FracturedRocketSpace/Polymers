#==============================================================================
# Main loop
#==============================================================================

import config as c
import numpy as np

# Initiate postition vector
r = np.zeros((c.nBeads ,2, c.nPolymers));
r[1,1,:] = c.linkDistance;


for k in range(c.nPolymers):    
    print(k)
    
