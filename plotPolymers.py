#==============================================================================
# Plot polymer
#==============================================================================


import matplotlib.pyplot as plt
import numpy as np
import config as c

def plotPolymers(polymers, endtoendDistances):
    plt.figure(1)
    plt.hold(True)
    plt.title('Polymer position')
    plt.xlabel('x')
    plt.ylabel('y')
    for i in range(len(polymers)):
        r = polymers[i]
        plt.plot( r[:,0], r[:,1] )

    plt.figure(2)
    plt.hold(True)
    plt.title('End-to-end distance as a function of the number of beads')
    plt.xlabel('Number of beads')
    plt.ylabel('End-to-end distance')
    x=np.linspace(1,c.nBeads, c.nBeads)
    for i in range(len(endtoendDistances)):
        y=endtoendDistances[i]
        plt.plot ( x, y)
        


