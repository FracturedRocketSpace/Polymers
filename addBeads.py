#==============================================================================
# Add beads
#==============================================================================

import config as c
import random as rand
import numpy as np
import math as m
import calculateEP


for i in range(2, c.nBeads):
    angleOffset=rand.randrange(0, 2*np.pi/c.nAngles)    #Generate random angle offset
    
    angles=np.zeros(c.nAngles)                          #Initialize angle vector
    w=np.zeros(c.nAngles)                               #Initialize weight vector
    p=np.zeros(c.nAngles)                               #Initialize probability vector

    for j in range (c.nAngles):                         #Calculate evenly spaced angles
        angles[j]=angleOffset+j*2*np.pi/c.nAngles
        E=calculateEP(r, angles, i, j, k)               #Calculate potential energy for new configuration
        w[j]=m.exp(-E/(c.kB*c.T))                       #Calculate weights for each angle
        W=+w[j]                                         #Calculate sum of the weights
        
    p[j]=w[j]/W                                         #Calculate probabilities
    #Add roulette wheel algorithm etc.
                                    