import random


class particle:

    def __init__(self, r, v):
        self.x = r[0]
        self.y = r[1]
        self.z = r[2]
        self.vx = v[0]
        self.vy = v[1]
        self.vz = v[2]

        self.q = -1.6e-19
        self.m = 1.0e-27

    def getR(self):
        return [self.x, self.y, self.z]

    def getV(self):
        return [self.vx, self.vy, self.vz]

    def getQ(self):
        return self.q

    def getM(self):
        return self.m



class randomDist:

    def __init__(self, xlims, ylims, zlims, vavg=1.0):
        self.xmin = xlims[0]
        self.xmax = xlims[1]

        self.ymin = ylims[0]
        self.ymax = ylims[1]

        self.zmin = zlims[0]
        self.zmax = zlims[1]

        self.v = vavg

    def GetDist(self):
        x = self.xmin + (random.random() * (self.xmax - self.xmin))
        y = self.ymin + (random.random() * (self.ymax - self.ymin))
        z = self.zmin + (random.random() * (self.zmax - self.zmin))

        vx = -0.5*self.v + random.random() * self.v
        vy = -0.5*self.v + random.random() * self.v
        vz = -0.5*self.v + random.random() * self.v

        r = [x,y,z]
        v = [vx,vy,vz]
        return [r, v]
