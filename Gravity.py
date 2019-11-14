"""

Calculating the Force on a particle

"""

from math import sqrt
import pdb
   
def TotalAcc(body, node, theta, Pot = 0):
    """
    Calculating the Force on a particle
    If Pot = 0: only forces are calulated and returned
    If Pot = 1: forces and potential energy is returned
    
    """
    x = body.x
    y = body.y
    if not node.children:
        if not node.points:
            if Pot == 0:
                return(0, 0)
            if Pot == 1:
                return(0, 0, 0)
        else:
            if body == node.points[0]:
                if Pot == 0:
                    return(0,0)
                if Pot == 1:
                    return(0 ,0 ,0)
            else:
                d = sqrt((x-node.ComX)**2+(y-node.ComY)**2)
                if Pot == 0:
                    return(TwoBodyAcc(x, y, node.ComX, node.ComY, d, node.TotM)
                           )
                if Pot == 1:
                    return(TwoBodyAcc(x, y, node.ComX, node.ComY, d, node.TotM)+
                           TwoBodyPotential(body.m, node.TotM, d)
                           )
    else:
        d = sqrt((x-node.ComX)**2+(y-node.ComY)**2)
        if d == 0:
            return(0, 0, 0)
        
        elif node.width/d < theta:            
            if Pot == 0:
                return(TwoBodyAcc(x, y, node.ComX, node.ComY, d, node.TotM)
                       )
            if Pot == 1:
                return(TwoBodyAcc(x, y, node.ComX, node.ComY, d, node.TotM)+
                       TwoBodyPotential(body.m, node.TotM, d)
                       )
        else:
            if Pot == 0:            
                TotF = [0,0]
                TotV = None
                for child in node.children:
                    TotF = [sum(n) for n in zip(TotF, TotalAcc(body, child, theta, Pot))]
            if Pot == 1:
                TotF = [0,0]
                TotV = 0
                for child in node.children:
                    Fx, Fy, V = TotalAcc(body, child, theta, Pot)
                    TotV = TotV + V
                    TotF[0] = TotF[0] + Fx 
                    TotF[1] = TotF[1] + Fy   
    
    TotFx = TotF[0]
    TotFy = TotF[1]
    return(TotFx, TotFy, TotV)

def TwoBodyPotential(m1, m2, d):
    """
    Acceleration of body 1
    """
    G = 1.184e-4 # G in Au^3/(yr^2* earth masses)
    
    V = -G*m1*m2*1./d
    return([V])    

def TwoBodyAcc(x1, y1, x2 ,y2, d, m2):
    """
    Acceleration of body 1
    """
    
    G = 1.184e-4 # G in Au^3/(yr^2* earth masses)
    
    Ax = G*m2*1./d**3*(x2-x1)
    Ay = G*m2*1./d**3*(y2-y1)
    return([Ax,Ay])
    
