#==============================================================================
# Main loop
#==============================================================================

import config as c
import numpy as np
from addPolymers import addPolymers
from minimizeEp import minimizeEp
import postProcess as post
import plotPolymers as plot

# Generate the polymers
polymers, polWeights, endtoendDistances = addPolymers();

#Post processing
# Calculates and sort the energy for use in the histogram
sortedEnergy = post.computeSortedEnergy(polymers);

# Ent to end distance
weightedEndtoendSq, weightedEndtoendSqStd = post.computeEndToEnd(endtoendDistances, polWeights);
print("End to End done")

# Gyradius
# The implementation is done this way because numba does not like lists.
gyradiusSq = np.zeros([len(polymers),c.nBeads])
# Loop over all polymers
for polNum in range(len(polymers)):
    gyradiusSq[polNum, :] = post.computeGyradius(polymers[polNum], polWeights[polNum], np.zeros(c.nBeads) , np.zeros([c.nBeads,2]));
    # Print progress
    if (polNum % 500 ==0):
        print('Gyradius polymer', polNum, 'done')
weightedGyradiusSq, weightedGyradiusSqStd = post.computeGyradiusStd(gyradiusSq, polWeights)
print("Gyradius done")

# Fit end to end and gyradius
fittedWeightedEndtoendSq = post.fitEndToEnd(weightedEndtoendSq);
fittedGyradius = post.fitGyradius(weightedGyradiusSq);

# Population size
popSize = post.computePopulation(polWeights);
print("Population calculated")

# Persistence length
lp1 = post.computePersistance(polymers, polWeights);
print("Persistence length calculated")

# Crossing detector
averageCrossings = post.computeAverageCrossings(polymers, polWeights);
print("Average number of crossings:", averageCrossings)

# Ep minimalisation
minEp=0;
if(c.minEp):
    minEp = minimizeEp(polymers);

# Add plot of some/all polymers
plot.plotPolymerPositions(polymers,polWeights);
plot.plotPersistenceLength(lp1);
plot.plotPopulationSize(popSize);
plot.plotEnergyHistogram(sortedEnergy);
plot.plotEndtoendSq(weightedEndtoendSq,weightedEndtoendSqStd,fittedWeightedEndtoendSq,popSize);
plot.plotGyradiusSq(weightedGyradiusSq,weightedGyradiusSqStd,fittedGyradius,popSize);
if(c.minEp):
    plot.plotMinEp(minEp);
plot.plotEndtoendSqPolymers(polymers,endtoendDistances);
if(c.savePlots):
    plot.savePlots();