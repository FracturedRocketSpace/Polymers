#==============================================================================
# Calculate eP
#==============================================================================

import config as c
import math as m

def calculateEP (r,angles,i,j,k):
    r[i,0,k]=r[i-1,0,k]+c.linkDistance*m.cos(angles[j])     # Position of the new bead for certain angle
    r[i,1,k]=r[i-1,1,k]+c.linkDistance*m.sin(angles[j])
    for a in range(i):                                      # Calculate interaction energy of new bead with each of the present beads
        dx=r[i,0,k]-r[a,0,k]
        dy=r[i,1,k]-r[a,1,k]
        r2=dx*dx+dy*dy
        r2i=1/r2
        r6i=r2i*r2i*r2i
        EP=+4*c.epsilon*((c.sigma**6*r6i)**2-c.sigma**6*r6i)  #Epsilon and sigma are not 1 in this problem. See book. Add up interaction energy with each present bead.
    return EP

    
 

#for p1 in range(config.nParticles):
#        for p2 in range(config.nParticles):
#            if p1 > p2:
#                # Calculate seperation and find nearest image
#                X = positions[p2,0] - positions[p1,0];
#                Y = positions[p2,1] - positions[p1,1];
#                Z = positions[p2,2] - positions[p1,2];
#                X -= np.rint(X/config.lCalc) * config.lCalc;
#                Y -= np.rint(Y/config.lCalc) * config.lCalc;
#                Z -= np.rint(Z/config.lCalc) * config.lCalc;
#                
#                # Calculate the total force, potential energy and virial
#                r2 = X*X + Y*Y + Z*Z;
#                r2i = 1 / r2;
#                r6i = r2i*r2i*r2i
#                
#                force = 24  * r6i * (2*r6i - 1) * r2i;
#                eP[i] += 4 * r6i * (r6i - 1);