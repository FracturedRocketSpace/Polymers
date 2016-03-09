#==============================================================================
# Plot polymer
#==============================================================================


import matplotlib.pyplot as plt


def plotPolymers(polymers):
    plt.figure(1)
    plt.hold(True)
    plt.title('Polymer position')
    plt.xlabel('x')
    plt.ylabel('y')
    for i in range(len(polymers)):
        r = polymers[i]
        plt.plot( r[:,0], r[:,1] )

