"""

Calculating the Force on a particle

"""

from math import sqrt

def Force(x, y, Tree, theta):
    TotFx = 0
    TotFy = 0
    
    
    return(Fx,Fy)
    
    
    # doesnt yet check if node is the same as the body itself
def CheckNode(x ,y, node, theta):
    if not node.children:
        if not node.points:
            return
        else:
            return(TwoBodyAcc(x, y, node.ComX, node.ComY, node.TotM))
    else:
        d = sqrt((x-node.ComX)**2+(y-node.ComY)**2)
        if node.width/d < theta:
            return(TwoBodyAcc(x, y, node.ComX, node.ComY, node.TotM))
        else:
            CheckNode(x, y, node.children[0], theta)
            CheckNode(x, y, node.children[1], theta)
            CheckNode(x, y, node.children[2], theta)
            CheckNode(x, y, node.children[3], theta)


def TwoBodyAcc(x1, y1, x2, y2, m2):
    """
    Acceleration of body 1
    """
    G = 1
    Fx = G*m2*1./sqrt((x2-x1)**2+(y2-y1)**2)**3*(x2-x1)
    Fy = G*m2*1./sqrt((x2-x1)**2+(y2-y1)**2)**3*(y2-y1)
    return(Fx,Fy)