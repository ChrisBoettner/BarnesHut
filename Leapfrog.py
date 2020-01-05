"""
Numerical Integrator

"""
import numpy as np
from OctoTree import OTree
import matplotlib. pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation, writers

def Integration(Tree, t_end, dt, theta, Pot = 0):
    """
    Leapfrog Integator
    
    Tree        : OctoTree to be used for force calculation      
    t_end       : End time of simulation
    dt          : Time step
    theta       : Accuracy measure
    Pot         : potential energy to be calculted? y/n
    """
    Pos_Evo     = np.empty([len(Tree.bodies), 3, len(np.arange(0,t_end+dt,dt))])
    Cons_Evo    = np.empty([len(np.arange(0,t_end+dt,dt)),3])
    O           = Tree
    
    n = 0
    for t in np.arange(0,t_end+dt,dt):
        O, Positions, Kin_E, Pot_E, AMom = Step(O, theta, dt, Pot)
        Pos_Evo[:,:,n] = Positions
        if Pot == 1:
            Cons_Evo[n,:] = Kin_E, Pot_E, AMom       
        n += 1
        
    if Pot == 0:
        return(Pos_Evo, None)
       
    if Pot == 1:      
        return(Pos_Evo, Cons_Evo)    


def Step(Tree, theta, dt, Pot = 0):
    Positions = np.empty([len(Tree.bodies), 3]) 
    
    if not Tree.bodies[0].ax:
        for n in range(len(Tree.bodies)):
            Tree.bodies[n].Update_Acc(Tree, theta, Pot)
    
    Ekin = None
    Epot = None
    AMom = None
            
    if Pot == 1:
        Ekin = O.calc_Ekin()
        Epot = O.calc_Epot()
        AMom = O.calc_AMom()
    
    for n in range(len(Tree.bodies)):
        pos = Tree.bodies[n].Update_Pos(dt)
        Positions[n,:] = pos
       
    O.update()
    
    for n in range(len(Tree.bodies)):
        Tree.bodies[n].ax_old = Tree.bodies[n].ax
        Tree.bodies[n].ay_old = Tree.bodies[n].ay
        Tree.bodies[n].az_old = Tree.bodies[n].az
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
    # sun earth moon
    #data = np.array([[0,1,0,-6.28,0,0,1],[0,0,0,0,0,0,332948],[0.002672,1,0,-6.27,-0.215,0,0.0123]])
    
    # periodic orbit // t_end  = 4.1761292190
    #mas_1 = 1
    #mas_2 = 1 
    #mas_3 = 0.5
    #pos_1 = np.array([-1,0,0])
    #pos_2 = -pos_1
    #pos_3 = np.array([0,0,0])
    #vel_1 = np.array([0.2869236336,0.0791847624,0])
    #vel_2 = vel_1
    #vel_3 = np.array([-2*vel_1[0]/mas_3,-2*vel_1[1]/mas_3,0])
    #data = np.array([[*pos_1,*vel_1,mas_1],[*pos_2,*vel_2,mas_2],[*pos_3,*vel_3,mas_3]])
    
    
    # solar system (vel in AU/day)
    Sun     = np.array([-3.832007941881110E-03, 7.431880253097086E-03, 2.395097284417125E-05,
                        -8.341984977273799E-06, -2.031943154612960E-06, 2.300637176641244E-07,
                        332900])
    Mercury = np.array([2.253097551375561E-02, -4.521457921378866E-01,-3.994851614208735E-02,
                        2.243656860429175E-02, 3.053835900923555E-03, -1.809047012661591E-03,
                        0.0553]) 
    Venus   = np.array([7.087246949289332E-01, 1.398154611393923E-01, -3.927914638625085E-02,
                        -3.778077813099616E-03, 1.979299280312870E-02, 4.893868670287644E-04,
                        0.815])
    Earth   = np.array([-2.386515578514603E-01, 9.622241452984469E-01, -1.657804954590847E-05,
                        -1.699194014801483E-02, -4.180818110741882E-03, 1.677666715049852E-07,
                        1])
    Mars    = np.array([-1.289882875979829E+00, -9.194852356024853E-01, 1.215494427135789E-02,
                        8.698035477086610E-03, -1.015667150375433E-02,-4.261712499600761E-04,
                        0.1074])
    Jupiter = np.array([5.520028408515556E-01, -5.189038800890025E+00, 9.171440063252453E-03,
                        7.412315334933891E-03, 1.157279202541518E-03,-1.706126866566051E-04,
                        317.8])
    Saturn  = np.array([3.812855875904778E+00, -9.272254556343645E+00, 9.436420173289210E-03,
                        4.849843275292275E-03, 2.105610128704437E-03, -2.297391872576618E-04,
                        95.2])
    Uranus  = np.array([1.621253938472379E+01, 1.139861571482007E+01, -1.677007583555294E-1,
                        -2.291114242051097E-03, 3.034197858121154E-03, 4.094754388027251E-05,
                        14.54])
    Neptune = np.array([2.924159647719122E+01, -6.347405223822967E+00, -5.431887489844970E-01,
                        6.455164412594517E-04, 3.086536946477342E-03, -7.861175595808179E-05,
                        17.15])
    
    data = np.array([Sun,Mercury,Venus,Earth,Mars,Jupiter,Saturn,Uranus,Neptune])
    O = OTree(data, 1)
    return(data,O)

if __name__ == "__main__":
     data, O = main()
     
     t_end  = 2*365
     dt     = 0.01
     Pos, Cons = Integration(O ,t_end ,dt ,0)
     
     #saving position data
     np.save("Position_data",Pos)
    


# =============================================================================
# plotting 2d
     fig = plt.figure()
     ax = fig.add_subplot(111)
     ax.set_aspect(aspect=1)
     for i in range(Pos.shape[0]):
         ax.plot(Pos[i,0,:],Pos[i,1,:])
     plt.show()
# =============================================================================
# plotting 3d
     #fig = plt.figure()
     #ax = fig.add_subplot(111, projection='3d')
     #for i in range(Pos.shape[0]):
         #ax.plot(Pos[i,0,:],Pos[i,1,:],Pos[i,2,:])
     #plt.show()

     



