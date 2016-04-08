#==============================================================================
# Calculate eP
#==============================================================================

import config as c
import math as m
from numba import jit

# Calculates EP between new bead and others
@jit
def calculateEP (r, angle, L):
    r[L, 0] = r[L-1, 0] + m.cos(angle)     # Position of the new bead for certain angle
    r[L, 1] = r[L-1, 1] + m.sin(angle)

    EP = 0;

    for a in range(L):                                      # Calculate interaction energy of new bead with each of the present beads
        dx= r[L,0] - r[a,0]
        dy= r[L,1] - r[a,1]

        r2=dx*dx + dy*dy
        r2i=1/r2
        r6i=r2i*r2i*r2i

        EP += 4*c.epsilon*((c.sigma**6*r6i)**2-c.sigma**6*r6i)  #Epsilon and sigma are not 1 in this problem. See book. Add up interaction energy with each present bead.

    return EP

# Calculates EP for whole polymer
@jit
def calculateEP2(r):
    EP = 0;

    for a in range(c.nBeads):                                      # Calculate interaction energy of new bead with each of the present beads
        for b in range(c.nBeads):
            if(a>b):
                dx= r[b,0] - r[a,0]
                dy= r[b,1] - r[a,1]

                r2=dx*dx + dy*dy
                r2i=1/r2
                r6i=r2i*r2i*r2i

                EP += 4*c.epsilon*((c.sigma**6*r6i)**2-c.sigma**6*r6i)  #Epsilon and sigma are not 1 in this problem. See book. Add up interaction energy with each present bead.

    return EP