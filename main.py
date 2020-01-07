import math, random
import numpy as np
# import matplotlib.pyplot as plt

from grid import grid3D
from spot import spot
from particle import particle



"""
1. Compute Charge Density
2. Compute Electric Potential
3. Compute Electric Field
4. Move Particles
5. Generate Particles - not needed
6. Output
7. Repeat
"""

"""
Setup: Create Particle distribution from dist. class in particle

sp_wt is ratio of macro to actual particles
"""

ActualParticles = 10000
MacroNumber = 100
sp_wt = ActualParticles/MacroNumber

T_e = 1.0e3
n0 = MacroNumber

xSize = 10
ySize = 10
zSize = 10

Grid = grid3D(xSize,ySize,zSize,parts=MacroNumber)

# Fill the grid with spots to help calculate charge density and fields
for i in range(len(Grid.g)):
    for j in range(len(Grid.g[0])):
        for k in range(len(Grid.g[0][0])):
            Grid.g[i][j][k] = YeeSpot(i,j,k)


# Create Particle dist.
rDist = randomDist( [Grid.xmin,Grid.xmax], [Grid.ymin,Grid.ymax], [Grid.zmin, Grid.zmax], vavg=2.0)

for i in range(MacroNumber):
        [r,v] = rDist.GetDist()
        part = particle(r,v)
        Grid.particles[i] = part


# Project the charges from macro particles onto YeeSpots
for i in range(MacroNumber):

    [x,y,z] = Grid.particles[i].getR()
    q = Grid.particles[i].getQ()

    i0 = math.floor(x); j0 = math.floor(y); k0 = math.floor(z)
    hx = x - i0; hy = y - j0; hz = z - k0;

    # Label specific, but if node moves in that direction then it loses the 1-x to just x
    w1 = (1-hx)*(1-hy)*(1-hz)*q
    w2 = (hx)*(1-hy)*(1-hz)*q
    w3 = (hx)*(hy)*(1-hz)*q
    w4 = (1-hx)*(hy)*(1-hz)*q
    w5 = (1-hx)*(hy)*(hz)*q
    w6 = (1-hx)*(1-hy)*(hz)*q
    w7 = (hx)*(1-hy)*(hz)*q
    w8 = (hx)*(hy)*(hz)*q

    Grid.g[i0][j0][k0].addQ(w1)
    Grid.g[i0+1][j0][k0].addQ(w2)
    Grid.g[i0+1][j0+1][k0].addQ(w3)
    Grid.g[i0][j0+1][k0].addQ(w4)
    Grid.g[i0][j0+1][k0+1].addQ(w5)
    Grid.g[i0][j0][k0+1].addQ(w6)
    Grid.g[i0+1][j0][k0+1].addQ(w7)
    Grid.g[i0+1][j0+1][k0+1].addQ(w8)


spaceing = 0.01
# Get Charges Dist
for i in range(len(Grid.g)):
    for j in range(len(Grid.g[0])):
        for k in range(len(Grid.g[0][0])):
            rho = Grid.g[i][j][k].q / spacing**3
            Grid.g[i][j][k].setRho(rho)


# Create YeeLattice from the YeeSpots made above ??? idk if this is needed
for i in range(len(Grid.cells)):
    for j in range(len(Grid.cells[0])):
        for k in range(len(Grid.cells[0][0])):
            spots = [Grid.g[i][j][k]), Grid.g[i+1][j][k], Grid.g[i][j+1][k], Grid.g[i+1][j+1][k], Grid.g[i][j][k+1], Grid.g[i+1][j][k+1], Grid.g[i][j+1][k+1], Grid.g[i+1][j+1][k+1]]
            Grid.cells[i][j][k] = YeeCube(spots,Grid.spacing)



N = xSize*ySize*zSize # flattened size i.e. i*j*K
A = [[0 for j in range(N)] for i in range(N)]
x = [0 for n in range(N)]
B = x

# For Making the A x = B to x = A' * B for solving potential

# Filling B

rhoScalar = 1/6 * spacing**3

for i in range(len(Grid.g)):
    for j in range(len(Grid.g[0])):
        for k in range(len(Grid.g[0][0])):
            index = (k*xSize*ySize) + (j*xSize) + i
            B[index] = Grid.g[i]][j][k].rho * rhoScalar

# Filling A
# for i in range(N):
#     for j in range(N):
for i in range(len(Grid.g)):
    for j in range(len(Grid.g[0])):
        for k in range(len(Grid.g[0][0])):

            n = (k*xSize*ySize) + (j*xSize) + i


            A[i][j] = 


# Solving for x
A = np.array(A)
B = np.array(B)
x = np.dot(A.transpose(), B)




"""
From here we should have the basic setup for a problem.

Next thing to do is start solving, using the discretized versions of our equations

we have the laplacian based on the values of each spot.

Gives:
    LAPLACIAN(phi) = phi_xx + phi_yy + phi_zz

    phi_xx[i][j][k] = ( phi[i+1][j][k] + phi[i-1][j][k] - 2*phi[i][j][k] ) / dx**2
    phi_yy[i][j][k] = ( phi[i][j+1][k] + phi[i][j-1][k] - 2*phi[i][j][k] ) / dy**2
    phi_zz[i][j][k] = ( phi[i][j][k+1] + phi[i][j][k-1] - 2*phi[i][j][k] ) / dz**2
"""

def laplacian(f):
    fpp = f
    for i in range(len(f)):
        for j in range(len(f[0])):
            for k in range(len(f[0][0])):

                fpp[i][j][k] =  ( phi[i+1][j][k] + phi[i-1][j][k] - 2*phi[i][j][k] ) / dx**2
                fpp[i][j][k] += ( phi[i][j+1][k] + phi[i][j-1][k] - 2*phi[i][j][k] ) / dy**2
                fpp[i][j][k] += ( phi[i][j][k+1] + phi[i][j][k-1] - 2*phi[i][j][k] ) / dz**2

    return fpp


"""
we know fpp =  - e/espilion * ( ni - n0*exp((phi - phi0)/k*Te) )

"""
































"""for each cube, you can do each laplacian of our phi equation """












"""
Charge Density:
    Using the filled in grid system with particles place everywhere, we can do sp_wt*q projection onto each spot surround the particle

    w4 = (hx)(hy)(hz)
    w3 = (1-hx)(hy)(hz)
    w2 = (hx)(1-hy)(hz)
    w1 = (1-hx)(1-hy)(hz)
    w8 = (hx)(hy)(1-hz)
    w7 = (1-hx)(hy)(1-hz)
    w6 = (hx)(1-hy)(1-hz)
    w5 = (1-hx)(1-hy)(1-hz)

    hx is fractional distance of the particel from the cell origin in x dir.

    Node volumes are cell volumes, dX*dY*dZ
"""







"""
Electric Potential:
    This is calculated by calling an elliptic equation solver. Used on the discretized Poisson equation

"""




"""
Electric Field:


"""






#eof
