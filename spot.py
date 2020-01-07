
class YeeSpot:

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

        self.q = 0.0
        self.rho = 0.0

    def addQ(self, q):
        self.q = self.q + q
    def setRho(self, rho):
        self.rho = rho

class YeeCube:

    def __init__(self, spots, spacing):
        # center posistion
        self.r = self.getPos(spots)
        # Edge Values
        self.Ex = 0.0
        self.Ey = 0.0
        self.Ez = 0.0
        # Face Values
        self.Hx = 0.0
        self.Hy = 0.0
        self.Hz = 0.0

        self.spacing = spacing

        self.calValues(spots)

    def getPos(self,spots):

        sumx,sumy,sumz = 0.,0.,0.
        for s in spots:
            sumx = sumx + s.x
            sumy = sumy + s.y
            sumz = sumz + s.z

        sumx = sumx/8
        sumy = sumy/8
        sumz = sumz/8
        return [sumx,sumy,sumz]

    def setE(self,E):
        self.Ex = E[0]
        self.Ey = E[1]
        self.Ez = E[2]

    def setH(self,H):
        self.Hx = H[0]
        self.Hy = H[1]
        self.Hz = H[2]


# class YeeSquare:
#
#     def __init__(self):
#


#eof
