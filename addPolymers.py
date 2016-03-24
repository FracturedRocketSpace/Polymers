#==============================================================================
# Add a new polymer
#==============================================================================

import config as c
import random as rand
import numpy as np
import math as m
from calculateEP import calculateEP

def calculateAvWeight(polWeight,alive,L):
    totalWeight = 0;
    for k in range(len(polWeight)):
        if(alive[k]):
            totalWeight += polWeight[k][L];
        
    return totalWeight/len(polWeight);

def enrich(polPositions,polWeights,endtoendDistance,alive,L, avWeight):
    for k in range(len(polPositions)):
        upLim = c.alphaUpLim * avWeight / polWeights[k][2];        
        
        if(polWeights[k][L] > upLim and alive[k]):
            newWeight = 0.5 * polWeights[k][L];
            polWeights[k][L] = newWeight;
            
            polPositions.append(np.copy(polPositions[k]));
            polWeights.append(np.copy(polPositions[k]));
            endtoendDistance.append(np.copy(polPositions[k]));
            alive.append(np.copy(alive[k]));

def prune(polPositions,polWeights,endtoendDistance,alive,L, avWeight):
    for k in range(len(polPositions)):
        downLim = c.alphaLowLim * avWeight / polWeights[k][2];        
        
        if(polWeights[k][L] < downLim and alive[k]):
            if(rand.uniform(0,1) < 0.5):
                newWeight = 2*polWeights[k][L];
                polWeights[k][L] = newWeight;
            else:
                alive[k] = False;
                

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
    
    polWeight[L] = polWeight[L-1] * W;

def addPolymers():
    #initialize polymer list
    polPositions = [];
    polWeights = [];
    endtoendDistances = [];
    alive = [];
    
    #fill the lists with empty polymers
    for k in range(c.nPolymers):
        polPositions.append(np.zeros([c.nBeads,2]));
        polWeights.append(np.ones([c.nBeads,1]));
        endtoendDistances.append(np.zeros([c.nBeads,2]));
        alive.append(True);
        
    #set initial beads
    for k in range(c.nPolymers):
        polPositions[k][1,1] = c.linkDistance;
        endtoendDistances[k][1,1]=c.linkDistance

    #generate polymers and save the values in lists
    for L in range(2,c.nBeads):
        #add beads
        for k in range(len(polPositions)):
            if(alive[k]):
                addBead(polPositions[k],polWeights[k],endtoendDistances[k],L);
            
        #enrich and prune
        avWeight = calculateAvWeight(polWeights,alive,L);
        prune(polPositions,polWeights,endtoendDistances,alive,L, avWeight);
        enrich(polPositions,polWeights,endtoendDistances,alive,L, avWeight);
    
    return polPositions, polWeights, endtoendDistances
