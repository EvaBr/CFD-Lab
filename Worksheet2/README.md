## Worksheet2: LBM Method
In this worksheet, we are going to simulate a 3D lid-driven cavity scenario by means of the Lattice Boltzmann Method (LBM).


Overview of the modules (Refer to Section 6 & 7 of the worksheet): 
* `helper.c`  contains functions for dynamical allocation and freeing of the memory for used variables (compare with the same file in Worksheet 1),
* `initLB.c` includes a function "readParameters" for input file reading (using helper functions from `helper.h`), and a function "initialiseFields" for initialization of particle distribution function fields collideField and streamField,
* `streaming.c` includes the function "doStreaming" for the streaming step,
* `collision.c` includes a function "doCollision" for the collision loop calculation, and a function "computePostCollisionDistributions" computing the post-collision distributions (f_i)* according to the BGK update rule using the equilibrium distribution values and the relaxation parameter,
* `computeCellValues.c` includes a function "computeDensity" computing the density within currentCell and storing the result in density, a function "computeVelocity" computing the velocity within currentCell and storing the result in velocity, and a function "computeFeq" computing the equilibrium distribution from the density and the velocity and storing the result in feq,
* `boundary.c` includes the function "treatBoundary", setting distribution functions inside boundary values,
* `visualLB.c` includes the function "writeVtkOutput", writing density and velocity from the collision field to a file. For this task, we can copy and adapt parts of the write_vtkFile function of Worksheet 1.
* `LBDefinitions.h` contains static const definitions for the lattice velocities, the lattice weights, etc., and
* `main.c` is the main program. 
