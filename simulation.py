import math, random, sys
import numpy as np

import pygame
from pygame.locals import *

from grid import grid2D
from spot import YeeSpot
from particle import particle, randomDist, ellipseDist



#===================#
pygame.init()
size = width, height = 750, 750
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

    if(i>0 and i<xSize-1 and j<ySize):
        Ex = 0.5*(phi[i-1][j] - phi[i+1][j])/spacing
    elif(i>0 and not i<xSize-1 and j<ySize):
        Ex = (phi[i-1][j] - phi[i][j])/spacing
    elif(i<xSize-1 and not i>0 and j<ySize):
        Ex = (phi[i][j] - phi[i+1][j])/spacing
    else:
        Ex = 0

    if(j>0 and j<ySize-1 and i<xSize):
        Ey = 0.5*(phi[i][j-1] - phi[i][j+1])/spacing
    elif(j>0 and not j<ySize-1 and i<xSize):
        Ey = (phi[i][j-1] - phi[i][j])/spacing
    elif(j<ySize-1 and not j>0 and i<xSize):
        Ey = (phi[i][j] - phi[i][j+1])/spacing
    else:
        Ey = 0

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

SingleParticle = False
count = 0
run = True
ActualParticles = 100
MacroNumber = 100
sp_wt = ActualParticles/MacroNumber

T_e = 1.0e3
n0 = MacroNumber

dt = 1.0e-10

xSize = 100
ySize = 100

xMax = 40
xMin = -40
yMax = 40
yMin = -40

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

if SingleParticle:
    r = [0.0, 0.0]
    v = [0.0, 1000.0]
    part = particle(r,v)
    Grid.particles[0] = part

else:
    rDist = randomDist( [-5,5], [-5,5], vavg=10000.0)
    # ellDist = ellipseDist([5,1],[1,1])

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
        # print(x,y)
        x = (8*x+width/2)
        y = (8*y+height/2)

        pygame.draw.circle(screen, white, (math.floor(x),math.floor(y)), 3)

    pygame.display.flip()

    # print average velcity
    avgvel = 0.0
    for i in range(MacroNumber):
        [vx,vy] = Grid.particles[i].getV()
        avgvel += math.sqrt(vx**2+vy**2)

    avgvel /= MacroNumber
    # print(avgvel)

#======== Charge Projection & Density Cal. ========#

    for i in range(MacroNumber):

        [x,y] = Grid.particles[i].getR()
        # print(x,y)
        if(x>xMin and y>yMin and x<xMax and y<yMax):

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
    # for i in range(len(Grid.g)):
    #     for j in range(len(Grid.g[0])):
    #         Grid.g[i][j].setPhi(phi[i][j])


#======== Half-Acceleration Move Particles ========#

    for i in range(MacroNumber):

        # print(x,y)

        [x,y] = Grid.particles[i].getR()
        [vx,vy] = Grid.particles[i].getV()
        q = Grid.particles[i].getQ()
        m = Grid.particles[i].getM()

        if(x>xMin and y>yMin and x<xMax and y<yMax):
            [Ex, Ey] = getE(x,y,phi,cell_spacing)


            Ez = 0.0
            Bx,By,Bz = 0.0, 0.0, 0.01 # Tesla


            #Boris Method

            """
            v0 = v_old + qdt/2m E
            v1 = v0 + 2* (v0 + v0 cross b0) / (1 + b0^2) ; b0 = qdt/2m B
            vnew = v1 + qdt/2m E
            """
            # -------------------

            E = np.array([Ex,Ey,Ez])
            B = np.array([Bx,By,Bz])
            v_old = np.array([vx,vy,0.0])

            v0 = np.add(v_old, np.multiply(E, q*dt/(2*m)))
            b0 = np.multiply(B,q*dt/(2*m))
            v1 = np.add( v0, np.cross( np.multiply( np.add(v0,np.cross(v0,b0)), 2.0/(1+np.linalg.norm(b0)**2) ), b0)  )
            v_new = np.add(v1, np.multiply(E, q*dt/(2*m)))

            vx_new = v_new[0]

            vy_new = v_new[1]

        else:
            vx_new = 0.0
            vy_new = 0.0

        x_new = x + vx_new * dt
        y_new = y + vy_new * dt

        Grid.particles[i].vx,Grid.particles[i].vy = vx_new, vy_new
        Grid.particles[i].x,Grid.particles[i].y = x_new, y_new


    if count%500 == 0:
        print(count)
    # strnum = "PIC%03d.png"%count
    # pygame.image.save(screen,strnum)


    count += 1


print("Done.")
pygame.quit()
