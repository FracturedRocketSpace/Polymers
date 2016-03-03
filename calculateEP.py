#==============================================================================
# Calculate eP
#==============================================================================

for p1 in range(config.nParticles):
        for p2 in range(config.nParticles):
            if p1 > p2:
                # Calculate seperation and find nearest image
                X = positions[p2,0] - positions[p1,0];
                Y = positions[p2,1] - positions[p1,1];
                Z = positions[p2,2] - positions[p1,2];
                X -= np.rint(X/config.lCalc) * config.lCalc;
                Y -= np.rint(Y/config.lCalc) * config.lCalc;
                Z -= np.rint(Z/config.lCalc) * config.lCalc;
                
                # Calculate the total force, potential energy and virial
                r2 = X*X + Y*Y + Z*Z;
                r2i = 1 / r2;
                r6i = r2i*r2i*r2i
                
                force = 24  * r6i * (2*r6i - 1) * r2i;
                eP[i] += 4 * r6i * (r6i - 1);