#==============================================================================
# Add a new polymer
#==============================================================================

import config as c
import random as rand
import numpy as np
import math as m
from calculateEP import calculateEP

# The pruning and enriching for PERM
def pruneANDenrich(polPositions,polWeights,endtoendDistance,alive,L,k, avWeight, avWeight3):
    lowLim = c.alphaLowLim * avWeight / avWeight3;
    upLim = c.alphaUpLim * avWeight / avWeight3;

    if(polWeights[k][L] < lowLim):
        if(rand.uniform(0,1) < 0.5):
            polWeights[k][L] *= 2;
        else:
            alive[k] = False;
            polWeights[k][L] = 0;
    elif(polWeights[k][L] > upLim and alive.count(True) < c.aliveLim):
        polWeights[k][L] /= 2;
        # Start new polymer at current bead
        polPositions.append(np.copy(polPositions[k]));
        polWeights.append(np.zeros([c.nBeads,1]));
        polWeights[len(polWeights)-1][L] = polWeights[k][L];
        endtoendDistance.append(np.copy(endtoendDistance[k]));
        alive.append(True);

# Calculates the discrete set of angles and their weights
def calculateAngles(r, L):
    angleOffset = rand.uniform(0, 2*np.pi/c.nAngles)   	 #Generate random angle offset
    angles = np.zeros(c.nAngles)                          #Initialize angle vector
    w = np.zeros(c.nAngles)                               #Initialize weight vector
    W=0;
    for j in range (c.nAngles):                         #Calculate evenly spaced angles
        angles[j] = angleOffset +j *2*np.pi/c.nAngles
        E = calculateEP(r, angles[j], L)               #Calculate potential energy for new configuration
        w[j] = m.exp(-E/(c.T))                       #Calculate weights for each angle
        W += w[j]
    return angles, w, W

# The roulette wheel algorithm for choosing an angle
def chooseAngle(w, W, angles):
    if(W==0):
        print('Problem with angle! Chance for each angle is 0.') # High chance that polymer is crossing
        return angles[rand.randrange(0, c.nAngles)];
    else:
        p=w/W;
        number = rand.random();
        for i in range(c.nAngles):
            # Upper lim
            upper = np.sum(p[0:i+1]);

            if( number <= upper ):
                return angles[i];
        # Just a failsafe
        print('Problem with angle! W is ', W);
        return angles[rand.randrange(0, c.nAngles)];

# Add a new bead to the polymer
def addBead(r, polWeight, endtoendDistance, L):
    # Calculate angles and weights
    angles, w, W = calculateAngles(r, L);

    # Choose Angle
    angle = chooseAngle(w, W, angles);

    # Add new bead to polymer
    r[L,0] = r[L-1,0] + m.cos(angle)
    r[L,1] = r[L-1,1] + m.sin(angle)

    # Calculate end-to-end distance
    endtoendDistance[L]=(r[L,0]**2+r[L,1]**2)

    # Update polWeight. Only include correction factor if doing PERM
    if (c.PERM==True):
        polWeight[L] = polWeight[L-1]* W/(0.75 * c.nAngles);
    else:
        polWeight[L] = polWeight[L-1]* W

# Generate all polymers
def addPolymers():
    #initialize polymer list
    #Lists of arrays are used so dynamically adding polymers is possible when enriching
    #In principle only lists could be used to save memory space, but arrays are more convenient to work with
    polPositions = [];
    polWeights = [];
    endtoendDistances = [];
    alive = [];
    avWeight = np.zeros([c.nBeads,2])

    # Generate polymers and save the values in lists
    k = 0;
    numCreated = 0;
    restarted = False;
    while numCreated < int(c.nPolymers) or restarted:
        restarted = False;
        # Start new polymer
        if(k == len(polPositions)):
            # Initialize storage arrays
            polPositions.append(np.zeros([c.nBeads,2]));
            polWeights.append(np.zeros([c.nBeads,1]));
            endtoendDistances.append(np.zeros([c.nBeads,1]));
            alive.append(True);
            # Set first values
            polWeights[k][0] = 1;
            polWeights[k][1] = 1;
            polPositions[k][1,1] = 1;
            endtoendDistances[k][1]= 1;

            numCreated += 1;
            start = 2;
        # Continue an already added polymer
        else:
            # Find starting position
            for L in range(1,c.nBeads):
                if(polWeights[k][L] > 0):
                    start = L + 1;
                    break;

        print('Iteration starts at: ', start)

        for L in range(start,c.nBeads):
            if(alive[k]):
                addBead(polPositions[k],polWeights[k],endtoendDistances[k],L);

                # Enrich and prune
                if(c.PERM):
                    # Update avWeight
                    avWeight[L,0] += polWeights[k][L];
                    avWeight[L,1] += 1;

                    pruneANDenrich(polPositions, polWeights, endtoendDistances, alive, L, k, avWeight[L,0], avWeight[2,0]);
                # Depending if we want to fix the population we restart the polymer or just stop building of polWeight has become zero
                elif(c.fixPop and polWeights[k][L] == 0):
                    restarted = True;
                    break;
                elif(polWeights[k][L] == 0):
                    alive[k]=False;
            else:
                print('Iteration stops at: ', L)
                break;

        if(restarted):
            print("Building of polymer ", k+1, "failed! Restarting build process.");
        else:
            print("Polymer", k+1, " done!\tPolymers: ", len(polPositions), "\tAlive:", alive.count(True),"\tStarted:", numCreated);
            print(" ");
            k+=1;
    # Return
    return polPositions, polWeights, endtoendDistances
