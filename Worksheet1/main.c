#include "helper.h"
#include "visual.h"
#include "init.h"
#include <stdio.h>


/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */
int main(int argn, char** args){
	//read the parameters
	read_parameters("cavity100.dat", Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, imax, 
			jmax, alpha, omg, tau, itermax, eps, dt_value)
	
	//set initial values
	double t = 0;
	int n = 0;
	init_uvp(UI, VI, PI, imax, jmax, U, V, P);
	
	//going through all time steps
	while(t<t_end){
		//adaptive time stepping
		calculate_dt(Re, tau, dt, dx, dy, imax, jmax, U, V);
		
		//setting bound.values
		boundaryvalues(imax, jmax, U, V);
		
		//computing F, G and right hand side of pressue eq.
		calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G);
		calculate_rs(dt, dx, dy, imax, jmax, F, G, RS);
		
		//iteration counter
		it = 0;
		
		do{
			//perform SOR iteration, at same time set bound.values for P and new residual value
			sor(omg, dx, dy, imax, jmax, P, RS, res);
			it++;
		}while(it<itermax && res>eps)
		
		//calculate U and V of this time step
		calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P);
		
		//indent time and number of time steps
		n++;
		t += dt;
	}
	
	//output of U, V, P for visualization
	/* if pics forall time steps needed, put this in the main loop... */
	write_vtkFile("DrivenCavity", n, xlength, ylength, imax, jmax, dx, dy, U, V, P);  
	
	return -1;
}
