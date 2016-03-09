#==============================================================================
# Plot polymer
#==============================================================================

import numpy as np
import config as c
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plotPolymers(polymers):
    plt.figure(1)
    plt.hold(True)
    plt.title('Polymer position')
    plt.xlabel('x')
    plt.ylabel('y')
    for i in range(len(polymers)):
        r = polymers[i]
        plt.plot( r[:,0], r[:,1] )

