#==============================================================================
# Add a new polymer
#==============================================================================

import config as c
import random as rand
import numpy as np
import math as m
from calculateEP import calculateEP

def calculateAngles(r, L):
    # Calculate angle weights
    angleOffset = rand.uniform(0, 2*np.pi/c.nAngles)   	 #Generate random angle offset
    angles = np.zeros(c.nAngles)                          #Initialize angle vector
    w = np.zeros(c.nAngles)                               #Initialize weight vector
    W=0;
    for j in range (c.nAngles):                         #Calculate evenly spaced angles
        angles[j] = angleOffset +j *2*np.pi/c.nAngles
        E = calculateEP(r, angles[j], L)               #Calculate potential energy for new configuration
        w[j] = m.exp(-E/(c.kB*c.T))                       #Calculate weights for each angle
        W += w[j]
    return angles, w, W

def chooseAngle(w, W, angles):
    p=w/W;
    number = rand.random();
    for i in range(c.nAngles):
        # Lower lim
        lower = 0;
        if(i > 0):
            lower = np.sum(p[0:i]);
        # Upper lim
        upper = lower + p[i]

        if( lower <= number < upper ):
            return angles[i];

    print('Problem with angle!');
    return 0;


def addBead(r, L, polWeight):
    # Calculate angles and weights
    angles, w, W = calculateAngles(r, L);

    # Choose Angle
    angle = chooseAngle(w, W, angles);

    # Add new bead!
    r[L,0] = r[L-1,0] + c.linkDistance*m.cos(angle)     # Position of the new bead for angle
    r[L,1] = r[L-1,1] + c.linkDistance*m.sin(angle)

    # Pruning
    if(c.pruning):
        print('Not implemented yet')

    # Do next recursive step
    if(L+1 < c.nBeads):
        polWeight *= W;
        addBead(r, L+1, polWeight)


def addPolymer(r, L, polWeight, polymers):
    # Build up new polymer using recursive addBead
    addBead(r, L, polWeight);
    # Save polymer
    polymers.append(r);