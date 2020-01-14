import math, random, sys

import pygame
from pygame.locals import *

from grid import grid2D
from spot import YeeSpot
from particle import particle, randomDist



#===================#
pygame.init()
size = width, height = 1000, 1000
screen = pygame.display.set_mode(size)
run = True

black = 0, 0, 0
red = 255, 0, 0
green = 0, 255, 0
blue = 0, 0, 255
white = 255, 255, 255
grey = 50,50,50
#===================#



#========== Kin. Sim ===========#


def getE(x,y,phi,spacing):
    """Given a position and a potential, this function
     will return the components of the Efield"""
    i = math.floor(x/spacing)
    j = math.floor(y/spacing)

    if(i>0 and i<xSize-1):
        Ex = 0.5*(phi[i-1][j] - phi[i+1][j])/spacing
    elif(i>0 and not i<xSize-1):
        Ex = (phi[i-1][j] - phi[i][j])/spacing
    elif(i<xSize-1 and not i>0):
        Ex = (phi[i][j] - phi[i+1][j])/spacing

    if(j>0 and j<ySize-1):
        Ey = 0.5*(phi[i][j-1] - phi[i][j+1])/spacing
    elif(j>0 and not j<ySize-1):
        Ey = (phi[i][j-1] - phi[i][j])/spacing
    elif(j<ySize-1 and not j>0):
        Ey = (phi[i][j] - phi[i][j+1])/spacing

    return [Ex,Ey]


# Need to implement These factors for it to actually be sim for plasma
# Plasma Freq
# w_p = math.sqrt(4*math.pi*n_e**2/m_e)
# Lam_d = V_therm / w_p
# Debye Length
# Skin Depth
# Lam_skin = c/w_p
# Larmor Freq
# w_c = q * B / (m_e*c)

count = 0
run = True
ActualParticles = 10000
MacroNumber = 10
sp_wt = ActualParticles/MacroNumber

T_e = 1.0e3
n0 = MacroNumber

dt = 0.001

xSize = 50
ySize = 50

xMax = 10
xMin = 0
yMax = 10
yMin = 0


#======== Create Environment ========#

#Using only symmetric
cell_spacing = (xMax-xMin) / xSize

Grid = grid2D(xSize,ySize,spacing=cell_spacing,parts=MacroNumber)
Grid.set_xlim(xMin,xMax)
Grid.set_ylim(yMin,yMax)

for i in range(len(Grid.g)):
    for j in range(len(Grid.g[0])):
            Grid.g[i][j] = YeeSpot(i,j)


#======== Make distribution ========#

# giving an average velocity we will get a unifrom distribution around that point
rDist = randomDist( [Grid.xmin,Grid.xmax], [Grid.ymin,Grid.ymax], vavg=2.0)

for i in range(MacroNumber):
    [r,v] = rDist.GetDist()
    part = particle(r,v)
    Grid.particles[i] = part


#======== Simulation Loop ========#

phi = [[0 for j in range(len(Grid.g[0]))] for i in range(len(Grid.g))]

while run:

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            sys.exit()
        elif event.type == pygame.KEYDOWN and event.key == K_ESCAPE:
            run = False

    screen.fill(black)
    # Draw Parts
    for i in range(MacroNumber):
        [x,y] = Grid.particles[i].getR()

        x = width * x/xMax
        y = height * y/yMax

        pygame.draw.circle(screen, white, (math.floor(x),math.floor(y)),3)

    pygame.display.flip()

#======== Charge Projection & Density Cal. ========#

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


    # Get Charges Dist
    for i in range(len(Grid.g)):
        for j in range(len(Grid.g[0])):

            rho = Grid.g[i][j].q / cell_spacing**3
            Grid.g[i][j].setRho(rho)


#======== Solve for Electric Potential ========#

    for i in range(xSize):
        for j in range(ySize):

            rho = Grid.g[i][j].rho

            if(i<xSize-1 and i>0 and j<ySize-1 and j>0):
                phi[i][j] = 0.25*(Grid.g[i-1][j].phi + Grid.g[i+1][j].phi + Grid.g[i][j-1].phi + Grid.g[i][j+1].phi + rho*cell_spacing**2)

            elif(not i<xSize-1 and i>0 and j<ySize-1 and j>0):
                phi[i][j] = 0.33*(Grid.g[i-1][j].phi + Grid.g[i][j-1].phi + Grid.g[i][j+1].phi + rho*cell_spacing**2)

            elif(i<xSize-1 and not i>0 and j<ySize-1 and j>0):
                phi[i][j] = 0.33*(Grid.g[i+1][j].phi + Grid.g[i][j-1].phi + Grid.g[i][j+1].phi + rho*cell_spacing**2)

            elif(i<xSize-1 and i>0 and j<ySize-1 and not j>0):
                phi[i][j] = 0.33*(Grid.g[i-1][j].phi + Grid.g[i+1][j].phi + Grid.g[i][j+1].phi + rho*cell_spacing**2)

            elif(i<xSize-1 and i>0 and not j<ySize-1 and j>0):
                phi[i][j] = 0.33*(Grid.g[i-1][j].phi + Grid.g[i+1][j].phi + Grid.g[i][j-1].phi + rho*cell_spacing**2)


    # Commit Potential
    for i in range(len(Grid.g)):
        for j in range(len(Grid.g[0])):
            Grid.g[i][j].setPhi(phi[i][j])


#======== Half-Acceleration Move Particles ========#

    for i in range(MacroNumber):

        [x,y] = Grid.particles[i].getR()
        [vx,vy] = Grid.particles[i].getV()
        q = Grid.particles[i].getQ()
        m = Grid.particles[i].getM()

        [Ex, Ey] = getE(x,y,phi,cell_spacing)

        """
        getEx and getEy are place holding functions which will return the E field
        at posistion of the particle. maybe E[math.floor(x)][math.floor(y)]???
        """

        vx_new = vx + (q/m) * Ex * dt
        vy_new = vy + (q/m) * Ey * dt

        x_new = x + vx_new * dt
        y_new = y + vy_new * dt

        Grid.particles[i].vx,Grid.particles[i].vy = vx_new, vy_new
        Grid.particles[i].x,Grid.particles[i].y = x_new, y_new


    if count%25 == 0:
        print(count)

    count += 1

print("Done.")
pygame.quit()
