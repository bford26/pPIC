import math, random
import numpy as np
# import matplotlib.pyplot as plt

from grid import grid2D
from spot import YeeSpot
from particle import particle, randomDist


# Plasma Freq
# w_p = math.sqrt(4*math.pi*n_e**2/m_e)
# Lam_d = V_therm / w_p
# Debye Length
# Skin Depth
# Lam_skin = c/w_p
# Larmor Freq
# w_c = q * B / (m_e*c)





ActualParticles = 10000
MacroNumber = 100
sp_wt = ActualParticles/MacroNumber

T_e = 1.0e3
n0 = MacroNumber

xSize = 10
ySize = 10

Grid = grid2D(xSize,ySize,parts=MacroNumber)

for i in range(len(Grid.g)):
    for j in range(len(Grid.g[0])):
            Grid.g[i][j] = YeeSpot(i,j)

# giving an average velocity we will get a unifrom distribution around that point
rDist = randomDist( [Grid.xmin,Grid.xmax], [Grid.ymin,Grid.ymax], vavg=2.0)


for i in range(MacroNumber):
    [r,v] = rDist.GetDist()
    part = particle(r,v)
    Grid.particles[i] = part


# This is where the simulation loop will begin


for i in range(MacroNumber):

    [x,y] = Grid.particles[i].getR()
    q = Grid.particles[i].getQ()

    i0 = math.floor(x); j0 = math.floor(y);
    hx = x - i0; hy = y - j0;

    w1 = (1-hx)*(1-hy)*q
    w2 = (hx)*(1-hy)*q
    w3 = (hx)*(hy)*q
    w4 = (1-hx)*(hy)*q

    Grid.g[i0][j0].addQ(w1)

    if i0 < xSize-1:
        Grid.g[i0+1][j0].addQ(w2)
        if j0 < ySize-1:
            Grid.g[i0+1][j0+1].addQ(w3)

    if j0 < ySize-1:
        Grid.g[i0][j0+1].addQ(w4)


spacing = Grid.spaceing
# Get Charges Dist
for i in range(len(Grid.g)):
    for j in range(len(Grid.g[0])):

        rho = Grid.g[i][j][k].q / spacing**3
        Grid.g[i][j].setRho(rho)

# because we are working in time steps the new potential will be phi
# as for finding phi[i][j] we use the grid.[][].phi because the time difference
phi = [[0 for j in range(len(Grid.g[0]))] for i in range(len(Grid.g))]

for i in range(xSize):
    for j in range(ySize):

    rho = Grid.g[i][j].getRho()

    if(i<xSize-1 and i>0 and j<ySize-1 and j>0):
        phi[i][j] = 0.25*(Grid.g[i-1][j].phi + Grid.g[i+1][j].phi + Grid.g[i][j-1].phi + Grid.g[i][j+1].phi + rho*spacing**2)

    elif(not i<xSize-1 and i>0 and j<ySize-1 and j>0):
        phi[i][j] = 0.33*(Grid.g[i-1][j].phi + Grid.g[i][j-1].phi + Grid.g[i][j+1].phi + rho*spacing**2)

    elif(i<xSize-1 and not i>0 and j<ySize-1 and j>0):
        phi[i][j] = 0.33*(Grid.g[i+1][j].phi + Grid.g[i][j-1].phi + Grid.g[i][j+1].phi + rho*spacing**2)

    elif(i<xSize-1 and i>0 and j<ySize-1 and not j>0):
        phi[i][j] = 0.33*(Grid.g[i-1][j].phi + Grid.g[i+1][j].phi + Grid.g[i][j+1].phi + rho*spacing**2)

    elif(i<xSize-1 and i>0 and not j<ySize-1 and j>0):
        phi[i][j] = 0.33*(Grid.g[i-1][j].phi + Grid.g[i+1][j].phi + Grid.g[i][j-1].phi + rho*spacing**2)

    #need to handel the boundaries

# solved phi push to grid spots this will be used later when nodal based
# objects are inserted to create elements
for i in range(len(Grid.g)):
    for j in range(len(Grid.g[0])):
        Grid.g[i][j].setPhi(phi[i][j])

"""
Now that the electric potential has been solved and store in the Grid object we
can solve for the E-field.

if bound:
    Ex,i = spacing*(phi[i+1][j] - phi[i][j])
else:
    Ex,i = 0.5*spacing*(phi[i-1][j] - phi[i+1][j])

if bound:
    Ey,j = spacing*(phi[i][j+1] - phi[i][j])
else:
    Ey,j = 0.5*spacing*(phi[i][j-1] - phi[i][j+1])

"""

Ex = [[0 for j in range(len(Grid.g[0]))] for i in range(len(Grid.g))]
Ey = [[0 for j in range(len(Grid.g[0]))] for i in range(len(Grid.g))]

for i in range(len(Grid.g)):
    for j in range(len(Grid.g[0])):

        if(i>0 and i<xSize-1):
            Ex[i][j] = 0.5*(phi[i-1][j] - phi[i+1][j])/spacing
        elif(i>0 and not i<xSize-1):
            Ex[i][j] = (phi[i-1][j] - phi[i][j])/spacing
        elif(i<xSize-1 and not i>0):
            Ex[i][j] = (phi[i][j] - phi[i+1][j])/spacing


        if(j>0 and j<ySize-1):
            Ey[i][j] = 0.5*(phi[i][j-1] - phi[i][j+1])/spacing
        elif(j>0 and not j<ySize-1):
            Ey[i][j] = (phi[i][j-1] - phi[i][j])/spacing
        elif(j<ySize-1 and not j>0):
            Ey[i][j] = (phi[i][j] - phi[i][j+1])/spacing



"""
With Efield moving particles using the Leapfrog method we see the equations

Note: this needs to be done for both coordinates:
    v_k+0.5 = v_k-0.5 + q/m * E * dt
    x_k+1 = x_k + v_k+0.5 * dt

this is the end of the loop and we can repeat the solving process
"""

for i in range(MacroNumber):

    [x,y] = Grid.particles[i].getR()
    [vx,vy] = Grid.particles[i].getV()
    q = Grid.particles[i].getQ()
    m = Grid.particles[i].getM()

    """
    getEx and getEy are place holding functions which will return the E field
    at posistion of the particle. maybe E[math.floor(x)][math.floor(y)]???
    """
    vx_new = vx + (q/m) * getEx(x) * dt
    vy_new = vy + (q/m) * getEy(y) * dt

    x_new = x + vx_new * dt
    y_new = y + vy_new * dt

    Grid.particles[i].vx,Grid.particles[i].vy = vx_new, vy_new
    Grid.particles[i].x,Grid.particles[i].y = x_new, y_new


















#eof
