"""
Numerical Integrator

"""
import numpy as np
from QuadTree import QTree
import matplotlib. pyplot as plt
from matplotlib.animation import FuncAnimation, writers

def Integration(Tree, t_end, dt, theta, Pot = 0):
    Pos_Evo     = np.empty([len(Tree.bodies),2,len(np.arange(0,t_end+dt,dt))])
    Cons_Evo    = np.empty([len(np.arange(0,t_end+dt,dt)),3])
    Q           = Tree
    
    n = 0
    for t in np.arange(0,t_end+dt,dt):
        Q, Positions, Kin_E, Pot_E, AMom = Step(Q, theta, dt, Pot)
        Pos_Evo[:,:,n] = Positions
        if Pot == 1:
            Cons_Evo[n,:] = Kin_E, Pot_E, AMom       
        n += 1
        
    if Pot == 0:
        return(Pos_Evo, None, None)
       
    if Pot == 1:      
        return(Pos_Evo, Cons_Evo)    


def Step(Tree, theta, dt, Pot = 0):
    Positions = np.empty([len(Tree.bodies), 2]) 
    
    if not Tree.bodies[0].ax:
        for n in range(len(Tree.bodies)):
            Tree.bodies[n].Update_Acc(Tree, theta, Pot)
    
    Ekin = None
    Epot = None
    AMom = None
            
    if Pot == 1:
        Ekin = Q.calc_Ekin()
        Epot = Q.calc_Epot()
        AMom = Q.calc_AMom()
    
    for n in range(len(Tree.bodies)):
        pos = Tree.bodies[n].Update_Pos(dt)
        Positions[n,:] = pos
       
    Q.update()
    
    for n in range(len(Tree.bodies)):
        Tree.bodies[n].ax_old = Tree.bodies[n].ax
        Tree.bodies[n].ay_old = Tree.bodies[n].ay
        Tree.bodies[n].Update_Acc(Tree, theta, Pot)     
        
    for n in range(len(Tree.bodies)):
        Tree.bodies[n].Update_Vel(dt)
        if Pot == 1:
          Tree.bodies[n].Update_Ekin()  
          Tree.bodies[n].Update_AMom()
    
    return(Tree, Positions, Ekin, Epot, AMom)
    
# =============================================================================
# Main
# =============================================================================
     
def main():
    #ata = np.array([[1,1,0,0,100000],[0,0,1.1,1,100000]])
    data = np.append(np.random.rand(5,4), np.random.rand(5,1)*1e+3,axis=1)
    #data = np.array([[0,1,-6.28,0,1],[0,0,0,0,332948.6]])
    Q = QTree(data, 1)
    return(data,Q)

if __name__ == "__main__":
     data, Q = main()
     
     t_end  = 1
     dt     = 0.00001
     Pos, Cons = Integration(Q ,t_end ,dt ,0 ,Pot = 1)
     


     fig = plt.figure()
     ax = fig.add_subplot(111)     
     #ax.set_aspect(aspect=1)
     ax.plot(Pos[0,0,:],Pos[0,1,:])
     ax.plot(Pos[1,0,:],Pos[1,1,:])
     ax.plot(Pos[2,0,:],Pos[2,1,:])
     ax.plot(Pos[3,0,:],Pos[3,1,:])
     ax.plot(Pos[4,0,:],Pos[4,1,:])
     #ax.plot(Pos[2,0,:],Pos[2,1,:])
     plt.show()

      


# =============================================================================
#      fig, ax = plt.subplots(figsize=(9, 7))
#      ax.set(xlim=(-1.1, 1.1), ylim=(-1.1, 1.1))
#       
#      scat = ax.scatter(Pos[:,0,0], Pos[:,1,0])
#        
#      def animate(i):
#          x_i = Pos[:,0,i]
#          y_i = Pos[:,1,i]
#          scat.set_offsets(np.c_[x_i, y_i])  
#          
#      anim = FuncAnimation(fig, animate, interval=0.1, frames=Pos.shape[2]-1) 
#  
#      plt.draw()
#      plt.show()
# =============================================================================
   
     #Writer = writers['ffmpeg']
     #writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
     
     #anim.save('tada.mp4', writer=writer)
     



