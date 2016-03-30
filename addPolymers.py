#==============================================================================
# Add a new polymer
#==============================================================================

import config as c
import random as rand
import numpy as np
import math as m
from calculateEP import calculateEP

def pruneANDenrich(polPositions,polWeights,endtoendDistance,alive,L,k, avWeight, avWeight3):
    lowLim = c.alphaLowLim * avWeight / avWeight3;
    upLim = c.alphaUpLim * avWeight / avWeight3;

    if(polWeights[k][L] < lowLim):
        if(rand.uniform(0,1) < 0.5):
            polWeights[k][L] *= 2;
            #print("Polymer not removed at step ", L);
        else:
            alive[k] = False;
            polWeights[k][L] = 0;
            #print("Polymer removed at step: ", L);

    elif(polWeights[k][L] > upLim and alive.count(True) < c.aliveLim):
        polWeights[k][L] /= 2;

        polPositions.append(np.copy(polPositions[k]));
        #TODO: In principle entry L should be zero as well for the new branch.
        polWeights.append(np.zeros([c.nBeads,1]));
        polWeights[len(polWeights)-1][L] = polWeights[k][L];
        endtoendDistance.append(np.copy(endtoendDistance[k]));

        alive.append(True);
        #print("Polymer doubled at step:", L);
    #else:
    #    print("Did nothing.")


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
    r[L,0] = r[L-1,0] + c.linkDistance*m.cos(angle)         # Position of the new bead for angle
    r[L,1] = r[L-1,1] + c.linkDistance*m.sin(angle)

    endtoendDistance[L]=(r[L,0]**2+r[L,1]**2)

    if (c.PERM==True):                                      # I think the correction factor is not necessary for normal Rosenbluth algorithm
        polWeight[L] = polWeight[L-1]* W/(0.75 * c.nAngles);
    else:
        polWeight[L] = polWeight[L-1]* W    
        
def addPolymers():
    #initialize polymer list
    polPositions = [];
    polWeights = [];
    endtoendDistances = [];
    alive = [];
    avWeight = np.zeros([c.nBeads,2])


    #generate polymers and save the values in lists
    k = 0;
    numCreated = 0;
    while numCreated < int(c.nPolymers):
        # Start new; Otherwise continue one already created.
        if(k == len(polPositions)):
            polPositions.append(np.zeros([c.nBeads,2]));
            polWeights.append(np.zeros([c.nBeads,1]));
            endtoendDistances.append(np.zeros([c.nBeads,1]));
            alive.append(True);

            polWeights[k][0] = 1;
            polWeights[k][1] = 1;
            polPositions[k][1,1] = c.linkDistance;
            endtoendDistances[k][1]=c.linkDistance**2

            numCreated += 1;

            start = 2;

        else:
            # Find starting position
            for L in range(2,c.nBeads):
                if(polWeights[k][L] > 0):
                    start = L + 1;
                    break;

        print('Iteration starts at: ', start)

        for L in range(start,c.nBeads):
            if(alive[k]):
                addBead(polPositions[k],polWeights[k],endtoendDistances[k],L);

                #enrich and prune
                if(c.PERM):
                    # Update avWeight
                    avWeight[L,0] += polWeights[k][L];
                    avWeight[L,1] += 1;


                    avWeight3 = avWeight[2,0];
                    avWeightL = avWeight[L,0];

                    pruneANDenrich(polPositions, polWeights, endtoendDistances, alive, L, k, avWeightL, avWeight3);
            else:
                print('Iteration stops at: ', L)
                break;

        print("polymer", k+1, " done!\tPolymers: ", len(polPositions), "\tAlive:", alive.count(True),"\tStarted:", numCreated);
        print(" ");
        k+=1;

    return polPositions, polWeights, endtoendDistances
