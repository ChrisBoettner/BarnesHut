"""

Calculating the Force on a particle

"""

from math import sqrt
    
def TotalAcc(body, node, theta):
    x = body.x
    y = body.y
    if not node.children:
        if not node.points:
            return([0,0])
        else:
            if body == node.points[0]:
                return([0,0])
            else:
                d = sqrt((x-node.ComX)**2+(y-node.ComY)**2)
                return(TwoBodyAcc(x, y, node.ComX, node.ComY, d, node.TotM))
    else:
        d = sqrt((x-node.ComX)**2+(y-node.ComY)**2)
        if node.width/d < theta:
            if (x == node.ComX and y == node.ComY):
                return([0,0])
            else:
                return(TwoBodyAcc(x, y, node.ComX, node.ComY, d, node.TotM))
        else:
            TotF = [0,0]
            for child in node.children:
                TotF = [sum(n) for n in zip(TotF, TotalAcc(body, child, theta))]
    return(TotF)
    

def TwoBodyAcc(x1, y1, x2 ,y2, d, m2):
    """
    Acceleration of body 1
    """
    G = 1
    Fx = G*m2*1./d**3*(x2-x1)
    Fy = G*m2*1./d**3*(y2-y1)
    return([Fx,Fy])