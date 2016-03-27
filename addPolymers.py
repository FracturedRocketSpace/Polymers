#==============================================================================
# Add a new polymer
#==============================================================================

import config as c
import random as rand
import numpy as np
import math as m
from calculateEP import calculateEP

def enrich(polPositions,polWeights,endtoendDistance,alive,L):
    enrichCounter = 0;
    # Sort polymers according to weight
    tempArray = np.array(polWeights)[:,L].flatten();
    order = np.argsort(tempArray);
    order = order[::-1];
    # Enrich
    while(alive.count(True) < c.nPolymers):
        # Set new weight
        polWeights[ order[enrichCounter] ][L] /= 2;
        # Multiply
        polPositions.append(np.copy(polPositions[ order[enrichCounter] ]));
        polWeights.append(np.copy(polWeights[ order[enrichCounter] ]));
        endtoendDistance.append(np.copy(endtoendDistance[ order[enrichCounter] ]));
        alive.append(True);
        #
        enrichCounter += 1;
    print("number of polymers duplicated: ", enrichCounter);

def prune(polPositions,polWeights,endtoendDistance,alive,L):
    pruneCounter = 0;
    # Sort polymers according to weight
    # Someone knows a very nice way to this without temporary arrays?
    tempArray = np.array(polWeights)[:,L].flatten();
    order = np.argsort(tempArray);
    # Do the pruning
    k = 0;
    while(pruneCounter < int(c.nPolymers*c.pruneFraction) ):
        if(alive[ order[k] ]):
            pruneCounter += 1;
            if(polWeights[ order[k] ][L] <= 0 or rand.uniform(0,1) < 0.5):
                alive[ order[k] ] = False;
                polWeights[ order[k] ][L] = 0;
            else:
                polWeights[ order[k] ][L] *= 2;
        k+=1;
    print("Number of prune steps: ", pruneCounter);


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
    if(W==0):
        print('Problem with angle! Chance for each angle is 0.') # High chance that polymer is crossing
        return angles[rand.randrange(0, c.nAngles)];
    else:
        p=w/W;
        number = rand.random();
        for i in range(c.nAngles):
            # Upper lim
            upper = np.sum(p[0:i]) + p[i]

            if( number <= upper ):
                return angles[i];
        # Could happen if w's nan
        print('Problem with angle! W is ', W);
        return angles[rand.randrange(0, c.nAngles)];


def addBead(r, polWeight, endtoendDistance, L):
    # Calculate angles and weights
    angles, w, W = calculateAngles(r, L);

    # Choose Angle
    angle = chooseAngle(w, W, angles);

    # Add new bead
    r[L,0] = r[L-1,0] + c.linkDistance*m.cos(angle)     # Position of the new bead for angle
    r[L,1] = r[L-1,1] + c.linkDistance*m.sin(angle)

    endtoendDistance[L,0]=m.sqrt(r[L,0]**2+r[L,1]**2)
    endtoendDistance[L,1]=r[L,0]**2+r[L,1]**2

    polWeight[L] = polWeight[L-1]* W;

def addPolymers():
    #initialize polymer list
    polPositions = [];
    polWeights = [];
    endtoendDistances = [];
    alive = [];

    #fill the lists with empty polymers and set initial beads
    for k in range(c.nPolymers):
        polPositions.append(np.zeros([c.nBeads,2]));
        polWeights.append(np.ones([c.nBeads,1]));
        endtoendDistances.append(np.zeros([c.nBeads,2]));
        alive.append(True);
        #
        polPositions[k][1,1] = c.linkDistance;
        endtoendDistances[k][1,1]=c.linkDistance

    #generate polymers and save the values in lists
    for L in range(2,c.nBeads):
        #add beads
        for k in range(len(polPositions)):
            if(alive[k]):
                addBead(polPositions[k],polWeights[k],endtoendDistances[k],L);

        # Enrich and prune
        if(c.PERM):
            prune(polPositions, polWeights, endtoendDistances, alive, L);
            enrich(polPositions, polWeights, endtoendDistances, alive, L);

        print("bead", L+1, " done!\tPolymers: ", len(polPositions), "\tAlive:", alive.count(True));
        print(" ");

    # remove all dead polymers
    alivePolPositions = [];
    alivePolWeights = [];
    aliveEndtoendDistances = [];

    for k in range(len(polPositions)):
        if(alive[k]):
            alivePolPositions.append(polPositions[k]);
            alivePolWeights.append(polWeights[k]);
            aliveEndtoendDistances.append(endtoendDistances[k]);

    return alivePolPositions, alivePolWeights, aliveEndtoendDistances
