#==============================================================================
# Config file
#==============================================================================

linkDistance    = 1;                                # Distance between two successive particles in the polymer chain


nBeads          = 250;                              # Length of the polymer in beads
T               = 1;                               # Temperature
nAngles         = 6;                                # Number of possible angles
nPolymers       = 100;                              # Number of polymers
aliveLim        = nPolymers*5;                      # Limits maximum number of alive polymers
fixPop          = False;                             # Only if not PERM. Fixes population to nPolymers by restarting failed polymers.

# PERM
PERM            = False;                            # Enables or disables the pruning/enriching part
alphaUpLim      = 2.2;                                # Constant for enrichining
alphaLowLim     = 1.2;                              # Constant for pruning

# Energy minimization algorithm
minEp           = False;                            # Set minimization on or off
minIter         = 50;                               # Nunber of iterations
enrichFrac      = int(1/6 * nPolymers);             # Fraction that is multipled per iteration
pruneFrac       = int(1/4 * nPolymers);             # Fraction that is removed per iteration
numMates        = int((pruneFrac - enrichFrac)/2);  # Number of mating events per iteration; By making it equal to (pruneFrac-enrichFrac)/2 the number of polymers stay equal
numMutations    = int(nPolymers/2);                 # Number of mutations per iteration

# Lennard Jones parameters
epsilon         = 0.25;                             # Value proposed in the book. Depth of energy potential well
sigma           = 0.8;                              # Value proposed in the book. Distance for which potential is zero

# Plot things
plotMaxPolymers = 10;                               # Maximum number of polymers plotted.