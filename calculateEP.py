#==============================================================================
# Calculate eP
#==============================================================================

import config as c
import math as m

def calculateEP (r, angle, L):
    r[L, 0] = r[L-1, 0] + c.linkDistance*m.cos(angle)     # Position of the new bead for certain angle
    r[L, 1] = r[L-1, 1] + c.linkDistance*m.sin(angle)
    
    EP = 0;
    
    for a in range(L):                                      # Calculate interaction energy of new bead with each of the present beads
        dx= r[L,0] - r[a,0]
        dy= r[L,1] - r[a,1]
        
        r2=dx*dx + dy*dy
        r2i=1/r2
        r6i=r2i*r2i*r2i
        
        EP += 4*c.epsilon*((c.sigma**6*r6i)**2-c.sigma**6*r6i)  #Epsilon and sigma are not 1 in this problem. See book. Add up interaction energy with each present bead.
    
    return EP