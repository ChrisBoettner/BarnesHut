"""

Calculating the Force on a particle

"""

from math import sqrt
    
def TotalAcc(x ,y, node, theta):
    if not node.children:
        if not node.points:
            return([0,0])
        else:
            if (x == node.ComX and y == node.ComY):
                return([0,0])
            else:
                return(TwoBodyAcc(x, y, node.ComX, node.ComY, node.TotM))
    else:
        d = sqrt((x-node.ComX)**2+(y-node.ComY)**2)
        if node.width/d < theta:
            if (x == node.ComX and y == node.ComY):
                return([0,0])
            else:
                return(TwoBodyAcc(x, y, node.ComX, node.ComY, node.TotM))
        else:
            TotF = [0,0]
            for child in node.children:
                TotF = [sum(n) for n in zip(TotF, TotalAcc(x, y, child, theta))]
    return(TotF)
    

def TwoBodyAcc(x1, y1, x2, y2, m2):
    """
    Acceleration of body 1
    """
    G = 1
    Fx = G*m2*1./sqrt((x2-x1)**2+(y2-y1)**2)**3*(x2-x1)
    Fy = G*m2*1./sqrt((x2-x1)**2+(y2-y1)**2)**3*(y2-y1)
    return([Fx,Fy])