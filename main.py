#==============================================================================
# Main loop
#==============================================================================

import config as c
from addPolymers import addPolymers
from plotPolymers import plotPolymers
from minimizeEp import minimizeEp
from postProcess import postProcess

# Generate the polymers
polymers, polWeights, endtoendDistances = addPolymers();

#Post processing
weightedEndtoendSq, weightedEndtoendSqStd,  weightedGyradiusSq, weightedGyradiusSqStd, popSize, lp1, fittedWeightedEndtoendSq, fittedGyradius, sortedEnergy = postProcess(polymers, polWeights, endtoendDistances);

# Ep minimalisation
minEp=0;
if(c.minEp):
    minEp = minimizeEp(polymers);

# Add plot of some/all polymers
plotPolymers(polymers, endtoendDistances, weightedEndtoendSq, weightedEndtoendSqStd, minEp, popSize, weightedGyradiusSq, weightedGyradiusSqStd, lp1, fittedWeightedEndtoendSq, fittedGyradius, sortedEnergy, polWeights)
