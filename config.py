#==============================================================================
# Config file
#==============================================================================

linkDistance    = 1;    # Distance between two successive particles in the polymer chain
nBeads          = 100;  # Length of the polymer in beads
T               = 10;    # Temperature
nAngles         = 6;    # Number of possible angles
nPolymers       = 10;   # Number of polymers

# PERN
pruning         = False; #

# Constants
kB              = 1;    # Boltzmann constant
mass            = 1;    # Particle mass

# Lennard Jones parameters
epsilon =   0.25;         # Value proposed in the book. Depth of energy potential well
sigma   =   0.8;          # Value proposed in the book. Distance for which potential is zero