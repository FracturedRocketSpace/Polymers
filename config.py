#==============================================================================
# Config file
#==============================================================================

linkDistance    = 1;                                # Distance between two successive particles in the polymer chain

nBeads          = 100;                              # Length of the polymer in beads
T               = 10;                               # Temperature
nAngles         = 6;                                # Number of possible angles
nPolymers       = 100;                              # Number of polymers

# PERM
PERM            = False;                            # Enables or disables the prunning/enriching part

# Energy minimization algorithm
minEp           = True;                             # Set minimization on or off
minIter         = 50;                               # Nunber of iterations
enrichFrac      = int(1/6 * nPolymers);             # Fraction that is multipled per iteration
pruneFrac       = int(1/4 * nPolymers);             # Fraction that is removed per iteration
numMates        = int((pruneFrac - enrichFrac)/2);  # Number of mating events per iteration; By making it equal to (pruneFrac-enrichFrac)/2 the number of polymers stay equal
numMutations    = int(nPolymers/2);                 # Number of mutations per iteration

# Constants
kB              = 1;                                # Boltzmann constant
mass            = 1;                                # Particle mass

# Lennard Jones parameters
epsilon         = 0.25;                             # Value proposed in the book. Depth of energy potential well
sigma           = 0.8;                              # Value proposed in the book. Distance for which potential is zero