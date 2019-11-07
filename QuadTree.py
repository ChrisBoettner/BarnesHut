"""
Creating the QuadTree

"""
from numpy import amax, amin
import matplotlib.pyplot as plt
from matplotlib import patches

import Force

# =============================================================================
# Classes
# =============================================================================


class Body():
    def __init__(self, x, y, vx=None, vy=None, m=None):
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy
        self.m  = m
    
    def Calc_Acc(self, QTree, theta):
        Ax, Ay = Force(self.x,self.y, QTree, theta)
        return(Ax,Ay)
        

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
        
        self.X_dim    = (amax(data[:,0])-amin(data[:,0]))
        self.Y_dim    = (amax(data[:,1])-amin(data[:,1]))
        self.X_origin = amin(data[:,0])-0.1*self.X_dim
        self.Y_origin = amin(data[:,1])-0.1*self.Y_dim
        self.X_dim    = (amax(data[:,0])-self.X_origin)*1.1
        self.Y_dim    = (amax(data[:,1])-self.Y_origin)*1.1
        
        self.root = Node(self.X_origin, self.Y_origin, self.X_dim, self.Y_dim, self.bodies)
        
        self.__subdivide()

    def __subdivide(self):
        recursive_subdivide(self.root, self.threshold)
    
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
# Helper functions
# =============================================================================

def recursive_subdivide(node, k):
    if len(node.points)<=k:
        return
    
    w_ = float(node.width/2)
    h_ = float(node.height/2)

    p = contains(node.x0, node.y0, w_, h_, node.points)
    sw = Node(node.x0, node.y0, w_, h_, p)
    recursive_subdivide(sw, k)

    p = contains(node.x0, node.y0+h_, w_, h_, node.points)
    nw = Node(node.x0, node.y0+h_, w_, h_, p)
    recursive_subdivide(nw, k)

    p = contains(node.x0+w_, node.y0, w_, h_, node.points)
    se = Node(node.x0 + w_, node.y0, w_, h_, p)
    recursive_subdivide(se, k)

    p = contains(node.x0+w_, node.y0+h_, w_, h_, node.points)
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

# =============================================================================
# Main
# =============================================================================


def main():
    from numpy import random
    
    data = random.rand(6,5)
    Q = QTree(data, 1)
    return(data,Q)

if __name__ == "__main__":
    data, Q = main()