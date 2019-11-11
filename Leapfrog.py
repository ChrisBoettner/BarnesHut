"""
Numerical Integrator

"""
import numpy as np


def Integration(Tree, t_end, dt, theta):
    PosEvo = np.empty([len(Tree.bodies),2,len(np.arange(0,t_end+dt,dt))])
    Q = Tree
    
    n = 0
    for t in np.arange(0,t_end+dt,dt):
        Q, Positions = Step(Q, theta, dt)
        PosEvo[:,:,n] = Positions
        Q.update()
        n += 1
    return(PosEvo)    



def Step(Tree, theta, dt):
    Positions = np.empty([len(Tree.bodies), 2])   
    
    if not Tree.bodies[0].ax:
        for n in range(len(Tree.bodies)):
            Tree.bodies[n].Update_Acc(Tree, theta)
    
    for n in range(len(Tree.bodies)):
        pos = Tree.bodies[n].Update_Pos(dt)
        Positions[n,:] = pos
        
    for n in range(len(Tree.bodies)):
        Tree.bodies[n].Update_Acc(Tree, theta)
            
    for n in range(len(Tree.bodies)):
        Tree.bodies[n].Update_Vel(dt)
    
    return(Tree, Positions)
        