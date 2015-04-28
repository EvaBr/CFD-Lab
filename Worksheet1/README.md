## Worksheet 1: Navier-Stokes Equations
Here we will model incompressible viscous flow over a driven cavity 
with no-slip boundary conditions on three sides and one moving wall.

We will simulate the flow using parameters given in the input file
_cavity100.dat_.

Additionally we will examine the influence of different weights on 
SOR solver's behaviour and time step sizes on algorithm's stability.
By increasing the Reynolds number, different viscosity scenarios will be 
simulated and compared.

Overview of the modules: 
* `helper.c`  contains functions for dynamical allocation and freeing of the 
memory for used variables,
* `init.c` includes a function for input file reading, and a function for
initialization of pressure and velocity arrays,
* `boundary_val.c` sets boundary values for velocity arrays,
* `uvp.c` consists of functions for calculating appropriate time step, 
right hand side of the pressure eq., new velocity arrays and discretized 
differential expressions of the momentum equations for the current 
iteration(timestep),
* `visual.c` contains functions for proper extraction of visualization data,
* `sor.c` holds the SOR solver and
* `main.c` is the main program. 
  
######################
testing commit: Wei
######################
