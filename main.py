#==============================================================================
# Main loop
#==============================================================================

import config as c
import numpy as np

# Initiate postition vector
r = np.zeros((c.nBeads ,2, c.nPolymers));
r[1,1,:] = 1;

# Add more particles
for i in range(10):    
    print(i)
    
