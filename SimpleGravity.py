from math import sqrt
import pdb


def TwoBodyAccEasy(x1, y1, x2 ,y2, d, m2):
    """
    Acceleration of body 1
    """
    G = 1.184e-4
    Ax = G*m2*1./d**3*(x2-x1)
    Ay = G*m2*1./d**3*(y2-y1)
    return([Ax,Ay])
    
def TotalAccEasy(body, Tree, theta, Pot):
    
    """
    Calculating the Force on a particle
    If Pot = 0: only forces are calulated and returned
    If Pot = 1: forces and potential energy is returned
    
    """
    x = body.x
    y = body.y
    
    TotAx = 0
    TotAy = 0
  
    for n in range(len(Tree.bodies)):
        if not Tree.bodies[n] == body:
          d = sqrt((x-Tree.bodies[n].x)**2+(y-Tree.bodies[n].y)**2)
          Ax, Ay = TwoBodyAccEasy(x, y, Tree.bodies[n].x ,Tree.bodies[n].y, d, Tree.bodies[n].m)
          TotAx  = TotAx + Ax
          TotAy  = TotAy + Ay
        else:
          continue
    
    return( TotAx, TotAy, None)
