"""
Numerical Integrator

"""
import numpy as np
from QuadTree import QTree
import matplotlib. pyplot as plt
from matplotlib.animation import FuncAnimation

def Integration(Tree, t_end, dt, theta):
    PosEvo = np.empty([len(Tree.bodies),2,len(np.arange(0,t_end+dt,dt))])
    Q = Tree
    
    n = 0
    for t in np.arange(0,t_end+dt,dt):
        Q, Positions = Step(Q, theta, dt)
        PosEvo[:,:,n] = Positions
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
        
    Q.update()
        
    for n in range(len(Tree.bodies)):
        Tree.bodies[n].Update_Acc(Tree, theta)
            
    for n in range(len(Tree.bodies)):
        Tree.bodies[n].Update_Vel(dt)
    
    return(Tree, Positions)
    
# =============================================================================
# Main
# =============================================================================
     
def main():
    data = np.array([[0,0,0,0,1000],[1,0,0,27,0.1]])
    Q = QTree(data, 1)
    return(data,Q)

if __name__ == "__main__":
    data, Q = main()
    T = Integration(Q,10,0.001,0)
    
    #plt.figure()
    #for n in range(T.shape[2]):
    #    plt.scatter(T[:,0,n],T[:,1,n])
    
    
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.set(xlim=(-10, 10), ylim=(-10, 10))
    
    scat = ax.scatter(T[:,0,0], T[:,1,0])
    
    def animate(i):
        x_i = T[:,0,i]
        y_i = T[:,1,i]
        scat.set_offsets(np.c_[x_i, y_i])
        
    anim = FuncAnimation(fig, animate, interval=0.1, frames=T.shape[2]-1)
 
    plt.draw()
    plt.show()