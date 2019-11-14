"""
Creating the QuadTree

"""
import matplotlib.pyplot as plt
from matplotlib import patches
import pdb

import Gravity
import SimpleGravity

# =============================================================================
# Classes
# =============================================================================


class Body():
    def __init__(self, x, y, vx, vy, m):
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy
        self.m  = m
        
        self.ax = None
        self.ay = None
        self.ax_old = None
        self.ay_old = None
        self.potE = None
        
        self.kinE = 0.5*self.m*(self.vx**2+self.vy**2)
        self.AMom = self.m*(self.x*self.vy-self.y*self.vx)
        
    
    def Update_Acc(self, QTree, theta, Pot = 0):
        self.ax, self.ay, self.potE = Gravity.TotalAcc(self, QTree.root, theta, Pot)
        #self.ax, self.ay, self.potE1 = SimpleGravity.TotalAccEasy(self, QTree, theta, Pot)
        ## direct summing
        
    def Update_Pos(self, dt):
        self.x = self.x + dt*self.vx + 0.5*self.ax*dt**2
        self.y = self.y + dt*self.vy + 0.5*self.ay*dt**2
        return([self.x,self.y])
        
    def Update_Vel(self, dt):
        self.vx = self.vx + 0.5*dt*(self.ax+self.ax_old)
        self.vy = self.vy + 0.5*dt*(self.ay+self.ay_old)
        
    def Update_Ekin(self):
        self.kinE = 0.5*self.m*(self.vx**2+self.vy**2)
        
    def Update_AMom(self):
        self.AMom = self.m*(self.x*self.vy-self.y*self.vx)
        
        

class Node():
    def __init__(self, x0, y0, w, h, points):
        self.x0 = x0
        self.y0 = y0
        self.width = w
        self.height = h
        self.points = points
        self.children = []
    
        self.ComX, self.ComY, self.TotM = self.Com()
        
    def Com(self):
        if not self.points:
            return(None,None,None)
        else:
            TotM = sum([self.points[n].m for n in range(len(self.points))])
            ComX = 1./TotM * sum([self.points[n].m*self.points[n].x for n in range(len(self.points))])
            ComY = 1./TotM * sum([self.points[n].m*self.points[n].y for n in range(len(self.points))])
            return(ComX,ComY,TotM)

class QTree():
    """
    Creates a QuadTree.
    
    data   : x and y coordinates, must be 2D array of shape [n,5]      
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
        X_max   = max([self.bodies[n].x for n in range(len(self.bodies))])
        Y_max   = max([self.bodies[n].y for n in range(len(self.bodies))])            

        X_dim     = X_max-X_min
        Y_dim     = Y_max-Y_min
        
        if X_dim == 0:
            X_dim = 0.1
        if Y_dim == 0:
            Y_dim = 0.1
        
        X_origin = X_min-0.1*X_dim
        Y_origin = Y_min-0.1*Y_dim
        
        X_dim    = (X_max-X_origin)*1.1
        Y_dim    = (Y_max-Y_origin)*1.1
        
        self.root = Node(X_origin, Y_origin, X_dim, Y_dim, self.bodies)
    
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
    
    def graph(self):
        fig, ax = plt.subplots(figsize=(12, 8))
        plt.title("Quadtree")
        leaves = find_leaves(self.root)
        for n in leaves:
            ax.add_patch(patches.Rectangle((n.x0, n.y0), n.width, n.height, fill=False))
        x = [point.x for point in self.bodies]
        y = [point.y for point in self.bodies]
        plt.plot(x, y, 'ro')
        plt.show()

# =============================================================================
# Helper functions for QuadTree
# =============================================================================

def recursive_subdivide(node, k):
    if len(node.points)<=k:
        return
    
    w_ = float(node.width/2)
    h_ = float(node.height/2)
    
    points = node.points

    p = contains(node.x0, node.y0, w_, h_, points)
    sw = Node(node.x0, node.y0, w_, h_, p)
    recursive_subdivide(sw, k)
    points = [n for n in points if n not in p] # remove already assigned points from list
    
    p = contains(node.x0, node.y0+h_, w_, h_, points)
    nw = Node(node.x0, node.y0+h_, w_, h_, p)
    recursive_subdivide(nw, k)
    points = [n for n in points if n not in p]

    p = contains(node.x0+w_, node.y0, w_, h_, points)
    se = Node(node.x0 + w_, node.y0, w_, h_, p)
    recursive_subdivide(se, k)
    points = [n for n in points if n not in p]

    p = contains(node.x0+w_, node.y0+h_, w_, h_, points)
    ne = Node(node.x0+w_, node.y0+h_, w_, h_, p)
    recursive_subdivide(ne, k)

    node.children = [nw, ne, se, sw]
    
    
def contains(x, y, w, h, points):
    pts = []
    for point in points:
        if point.x >= x and point.x <= x+w and point.y>=y and point.y<=y+h:
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