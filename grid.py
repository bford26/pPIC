

class grid3D:

    def __init__(self,Nx,Ny,Nz,spacing=None,parts=10):
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz

        self.spacing = 0.0001 # Update this value to match Debye length?

        self.xmin = -10
        self.xmax = 10
        self.ymin = -10
        self.ymax = 10
        self.zmin = -10
        self.zmax = 10

        self.particles = [None]*parts

        self.g = [[[ None for k in range(Nz)] for j in range(Ny)] for i in range(Nx)]
        self.cells = [[[ None for k in range(Nz-1)] for j in range(Ny-1)] for i in range(Nx-1)]
