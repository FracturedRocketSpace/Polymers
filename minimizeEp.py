#==============================================================================
# Minimize Ep
#==============================================================================

import config as c
import numpy as np
import random as rand
from calculateEP import calculateTotalEP
import math as m

def pePolymers(polymers,Ep):
    order = np.argsort(Ep)
    polymers2 = [];
    # Make sure to keep population constant
    extra=0;
    if(len(polymers)>c.nPolymers):
        extra = len(polymers) - c.nPolymers;
    for j  in range(len(polymers)):
        if(j < c.enrichFrac ):
            polymers2.append( polymers[order[j]] )
        if(j < len(polymers) - c.pruneFrac - extra ):
            polymers2.append( polymers[order[j]] )
    return polymers2;

def matePolymers(polymers):
    # Randomly choose two polymers
    n1 = rand.randrange(0,len(polymers))
    n2 = rand.randrange(0,len(polymers))

    # Choose cut position
    cut = rand.randrange(0,c.nBeads)

    # Create first offspring
    r1 = polymers[n1][0:cut];
    r2 = polymers[n2][cut : :] - polymers[n2][cut-1] + polymers[n1][cut-1] ; # To make sure parts connect well
    rNew = np.zeros((c.nBeads,2))
    rNew[0:cut,:] = r1;
    rNew[cut : :,:] = r2;
    polymers.append(rNew);

    # Create second offspring
    r1 = polymers[n2][0:cut];
    r2 = polymers[n1][cut : :] - polymers[n1][cut-1] + polymers[n2][cut-1];
    rNew = np.zeros((c.nBeads,2))
    rNew[0:cut,:] = r1;
    rNew[cut : :,:] = r2;
    polymers.append(rNew);

def mutatePolymer(polymers):
    # Choose mutation
    mutation = rand.randrange(0,2)
    if(mutation==0):
        # Move one bead around
        pol = rand.randrange(0,len(polymers)) # Randomly choose a polymer
        bead = rand.randrange(1,c.nBeads-1) # Randomly choose a bead at which the angle is changed
        # Break polymer
        r1 = polymers[pol][0:bead];
        r2 = polymers[pol][bead+1 : :] - polymers[pol][bead];
        # Choose new angle
        angle = rand.uniform(0, 2*np.pi)
        # Build mutated polymer
        rNew = np.zeros((c.nBeads,2))
        rNew[0:bead,:] = r1;
        rNew[bead,0] = rNew[bead-1,0] + m.cos(angle);
        rNew[bead,1] = rNew[bead-1,1] + m.sin(angle);
        rNew[bead+1 : :,:] = r2 + rNew[bead];
        # Insert
        polymers[pol] = rNew;
    else:
        # Randomize a single angle
        pol = rand.randrange(0,len(polymers)) # Randomly choose a polymer
        bead = rand.randrange(1,c.nBeads-1) # Randomly choose a bead at which the angle is changed
        # Break polymer
        r1 = polymers[pol][0:bead];
        r2 = polymers[pol][bead : :] - polymers[pol][bead-1];
        r2new = polymers[pol][bead : :] - polymers[pol][bead-1];
        # Choose angle over wich to rotate second part
        angle = rand.uniform(0, 2*np.pi)
        # Rotate
        r2new[:,0] = r2[:,0]*m.cos(angle) - r2[:,1]*m.sin(angle);
        r2new[:,1] = r2[:,0]*m.sin(angle) + r2[:,1]*m.cos(angle);
        # Build mutated polymer
        rNew = np.zeros((c.nBeads,2))
        rNew[0:bead] = r1;
        rNew[bead : :] = r2new + rNew[bead-1];
        # Insert
        polymers[pol] = rNew;


# The genetic algorithm
def minimizeEp(polymers):
    minEp = np.zeros(c.minIter)

    for i in range(c.minIter):
        # Calculate EP for all
        Ep = []; # Reset
        for j in range(len(polymers)):
            Ep.append( calculateTotalEP(polymers[j]) )

        # Set minimum EP at start of this Iter
        minEp[i]=min(Ep)
        # Save best polymer
        if( min(minEp[0:i+1])==minEp[i] ):
            rBest=polymers[ Ep.index(min(Ep)) ];

        # Prune and enrich
        polymers = pePolymers(polymers,Ep)

        # Mate some pairs
        for j in range(c.numMates):
            matePolymers(polymers)

        # Mutate
        for j in range(c.numMutations):
            mutatePolymer(polymers);

        # Reinsert best polymer
        polymers.append(rBest);

        # Print
        print('Minimization step', i+1, 'completed')

    return minEp