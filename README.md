# Barnes-Hut N-Body Simulation

This project implements a Barnes-Hut algorithm for simulating the gravitational interaction of N bodies in a 3D space.  The Barnes-Hut algorithm offers a significant performance improvement over direct summation methods, especially for large N, by treating distant clusters of particles as single entities. This reduces the computational complexity from O(N²) to O(N log N).

## Algorithm

The core of the Barnes-Hut algorithm lies in recursively dividing the 3D simulation space into octants (octree). Each node in the octree represents a region of space and contains either a single body or a center of mass and total mass representing a cluster of bodies within that region.

When calculating the force on a particular body, the algorithm traverses the octree. For each node encountered:

1. If the node is a leaf (containing a single body) and that body is not the one being considered, the direct gravitational force between the two bodies is calculated.
2. If the node represents a cluster of bodies, the algorithm checks if the distance between the body and the node's center of mass is sufficiently large compared to the node's size (using a threshold value `theta`). 
    - If the distance is large enough, the cluster is treated as a single point mass located at the center of mass, and the force between the body and this point mass is calculated.
    - If the distance is too small, the algorithm recursively descends into the node's children, repeating the process for each child.

This approach significantly reduces the number of force calculations required for distant bodies, thereby speeding up the simulation.

## Code Structure

* **`OctoTree.py`**:  Implements the octree data structure.
    * `Body` class: Represents a celestial body with properties like position, velocity, mass, and acceleration.
    * `Node` class: Represents a node in the octree, storing information about the bodies or clusters it contains.
    * `OTree` class: Manages the construction and update of the octree.
* **`Gravity.py`**: Contains the functions for calculating gravitational forces and potential energy.
    * `TotalAcc`: Calculates the total acceleration on a body due to all other bodies in the system using the Barnes-Hut algorithm.
    * `TwoBodyAcc`: Calculates the gravitational acceleration between two individual bodies.
    * `TwoBodyPotential`: Computes the gravitational potential energy between two bodies.
* **`Leapfrog.py`**: Implements the leapfrog numerical integration method for updating the positions and velocities of the bodies over time.
    * `Integration`: Performs the main simulation loop, calling the `Step` function to advance the simulation by one timestep.
    * `Step`: Updates the positions, velocities, and accelerations of all bodies using the leapfrog algorithm.
* **`SimpleGravity.py`**: Provides a direct summation gravity calculation (O(N²)) for comparison and testing.
* **`Leapfrog.py`**:  Implements the Leapfrog integrator to evolve the system forward in time.
* **`README.md`**:  This file.

## Usage

The `main` function in `Leapfrog.py` sets up the initial conditions of the simulation (currently configured for a solar system-like setup). You can modify the `data` array in `main` to simulate different systems. The simulation parameters, such as the end time (`t_end`), timestep (`dt`), and accuracy parameter (`theta`), can also be adjusted.  The results are saved to `Position_data.npy`.

Example usage in `Leapfrog.py`:

```python
t_end  = 2*365  # Simulation end time (in days)
dt     = 0.01   # Timestep (in days)
Pos, Cons = Integration(O, t_end, dt, 0)  # Run the simulation

np.save("Position_data", Pos) # Save the position data

# Plotting (example 2D plot)
fig = plt.figure()
# ... (plotting code)
plt.show()
```
