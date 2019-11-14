"""
Creating the QuadTree

"""

import Gravity
#import SimpleGravity

# =============================================================================
# Classes
# =============================================================================


class Body():
    def __init__(self, x, y, z, vx, vy, vz, m):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.m  = m
        
        self.ax = None
        self.ay = None
        self.az = None
        self.ax_old = None
        self.ay_old = None
        self.az_old = None
        self.potE = None
        
        self.kinE = 0.5*self.m*(self.vx**2+self.vy**2+self.vz**2)
        self.AMom = self.m*(self.x*self.vy-self.y*self.vx) # remove this later
        
    
    def Update_Acc(self, QTree, theta, Pot = 0):
        self.ax, self.ay, self.az, self.potE = Gravity.TotalAcc(self, QTree.root, theta, Pot)
        #self.ax, self.ay, self.potE1 = SimpleGravity.TotalAccEasy(self, QTree, theta, Pot)
        ## direct summing algorithm
        
    def Update_Pos(self, dt):
        self.x = self.x + dt*self.vx + 0.5*self.ax*dt**2
        self.y = self.y + dt*self.vy + 0.5*self.ay*dt**2
        self.z = self.z + dt*self.vz + 0.5*self.az*dt**2
        return([self.x, self.y, self.z])
        
    def Update_Vel(self, dt):
        self.vx = self.vx + 0.5*dt*(self.ax+self.ax_old)
        self.vy = self.vy + 0.5*dt*(self.ay+self.ay_old)
        self.vz = self.vz + 0.5*dt*(self.az+self.az_old)
        
    def Update_Ekin(self):
        self.kinE = 0.5*self.m*(self.vx**2+self.vy**2+self.vz**2)
        
    def Update_AMom(self):  # remove this later
        self.AMom = self.m*(self.x*self.vy-self.y*self.vx)
        
        

class Node():
    def __init__(self, x0, y0, z0, w, h, d, points):
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.width = w
        self.height = h
        self.depth = d
        self.points = points
        self.children = []
    
        self.ComX, self.ComY, self.ComZ, self.TotM = self.Com()
        
    def Com(self):
        if not self.points:
            return(None, None, None, None)
        else:
            TotM = sum([self.points[n].m for n in range(len(self.points))])
            ComX = 1./TotM * sum([self.points[n].m*self.points[n].x for n in range(len(self.points))])
            ComY = 1./TotM * sum([self.points[n].m*self.points[n].y for n in range(len(self.points))])
            ComZ = 1./TotM * sum([self.points[n].m*self.points[n].z for n in range(len(self.points))])
            return(ComX, ComY, ComZ, TotM)

class OTree():
    """
    Creates a OctoTree.
    
    data   : position (3), velocity (3) and m, must be 5D array of shape [n,7]      
    k      : number of max points per node (integer)
    points
    """
    
    def __init__(self, data, k):        
        self.threshold = k
        self.bodies      = [Body(*[data[i,n] for n in range(data.shape[1])])                            
                             for i in range(data.shape[0])]
        
        self.__build_root()
        self.__subdivide()

    def __subdivide(self):
        recursive_subdivide(self.root, self.threshold)
        
    def __build_root(self):
        X_min   = min([self.bodies[n].x for n in range(len(self.bodies))])
        Y_min   = min([self.bodies[n].y for n in range(len(self.bodies))])
        Z_min   = min([self.bodies[n].z for n in range(len(self.bodies))])
        X_max   = max([self.bodies[n].x for n in range(len(self.bodies))])
        Y_max   = max([self.bodies[n].y for n in range(len(self.bodies))])            
        Z_max   = max([self.bodies[n].z for n in range(len(self.bodies))])
        
        X_dim     = X_max-X_min
        Y_dim     = Y_max-Y_min
        Z_dim     = Z_max-Z_min
        
        if X_dim == 0:
            X_dim = 0.1
        if Y_dim == 0:
            Y_dim = 0.1
        if Z_dim == 0:
            Z_dim = 0.1
        
        X_origin = X_min-0.1*X_dim
        Y_origin = Y_min-0.1*Y_dim
        Z_origin = Z_min-0.1*Z_dim
        
        X_dim    = (X_max-X_origin)*1.1
        Y_dim    = (Y_max-Y_origin)*1.1
        Z_dim    = (Z_max-Z_origin)*1.1
        
        self.root = Node(X_origin, Y_origin, Z_origin, X_dim, Y_dim, Z_dim, self.bodies)
    
    def update(self):
        self.__build_root()
        self.__subdivide()
    
    def calc_Ekin(self):
        Ekin = sum([self.bodies[n].kinE for n in range(len(self.bodies))])
        return(Ekin)
        
    def calc_Epot(self):
        Epot = 0.5*sum([self.bodies[n].potE for n in range(len(self.bodies))])
        return(Epot)
        
    def calc_AMom(self):
        AMom = sum([self.bodies[n].AMom for n in range(len(self.bodies))])
        return(AMom)

# =============================================================================
# Helper functions for QuadTree
# =============================================================================

import pdb

def recursive_subdivide(node, k):
    if len(node.points)<=k:
        return
    
    w_ = float(node.width/2)
    h_ = float(node.height/2)
    d_ = float(node.depth/2)
    
    points = node.points

        
    p = contains(node.x0, node.y0+h_, node.z0, w_, h_, d_, points)
    nwf = Node(node.x0, node.y0+h_, node.z0, w_, h_, d_, p)
    recursive_subdivide(nwf, k)
    points = [n for n in points if n not in p]

    p = contains(node.x0+w_, node.y0+h_, node.z0, w_, h_, d_, points)
    nef = Node(node.x0+w_, node.y0+h_, node.z0, w_, h_, d_, p)
    recursive_subdivide(nef, k)
    points = [n for n in points if n not in p]
    
    p = contains(node.x0, node.y0, node.z0, w_, h_, d_, points)
    swf = Node(node.x0, node.y0, node.z0, w_, h_, d_, p)
    recursive_subdivide(swf, k)
    points = [n for n in points if n not in p]
    
    p = contains(node.x0+w_, node.y0, node.z0, w_, h_, d_, points)
    sef = Node(node.x0+w_, node.y0, node.z0, w_, h_, d_, p)
    recursive_subdivide(sef, k)
    points = [n for n in points if n not in p]

    p = contains(node.x0, node.y0+h_, node.z0+d_, w_, h_, d_, points)
    nwb = Node(node.x0, node.y0+h_, node.z0+d_, w_, h_, d_, p)
    recursive_subdivide(nwb, k)
    points = [n for n in points if n not in p]

    p = contains(node.x0+w_, node.y0+h_, node.z0+d_, w_, h_, d_, points)
    neb = Node(node.x0+w_, node.y0+h_, node.z0+d_, w_, h_, d_, p)
    recursive_subdivide(neb, k)
    points = [n for n in points if n not in p]
    
    p = contains(node.x0, node.y0, node.z0+d_, w_, h_, d_, points)
    swb = Node(node.x0, node.y0, node.z0+d_, w_, h_, d_, p)
    recursive_subdivide(swb, k)
    points = [n for n in points if n not in p]
    
    p = contains(node.x0+w_, node.y0, node.z0+d_, w_, h_, d_, points)
    seb = Node(node.x0+w_, node.y0, node.z0+d_, w_, h_, d_, p)
    recursive_subdivide(seb, k)
    

    node.children = [nwf, nef, swf, sef, nwb, neb, swb, seb]
    
    
def contains(x, y, z, w, h, d, points):
    pts = []
    for point in points:
        if point.x >= x and point.x <= x+w and point.y>=y and point.y<=y+h and point.z >= z and point.z <= z+d:
            pts.append(point)
    return(pts)

def find_leaves(node):
    if not node.children:
        return [node]
    else:
        leaves = []
        for child in node.children:
            leaves += (find_leaves(child))
    return(leaves)