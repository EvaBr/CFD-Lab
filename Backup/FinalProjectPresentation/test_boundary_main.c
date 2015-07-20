#include "helper.h"
#include "visual.h"
#include "init.h"
#include "sor.h"
#include "boundary_val.h"
#include "uvp.h"
#include <stdio.h>


int main(int argn, char** args){
 	if (argn !=2 ) {
        	printf("When running the simulation, please give a valid scenario file name!\n");
        	return 1;
        }
	//set the scenario
	char *filename = NULL;
	filename = args[1];

	//initialize variables
	double t = 0; /*time start*/
	int it, n = 0; /*iteration and time step counter*/
	double res; /*residual for SOR*/
		/*arrays*/
	double **U, **V, **P;
	double **RS, **F, **G;
	int **Flag; //additional data structure for arbitrary geometry
		/*those to be read in from the input file*/
	double Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, alpha, omg, tau, eps, dt_value;
	int  imax, jmax, itermax;

	double presLeft, presRight, presDelta; //for pressure stuff
	int wl, wr, wt, wb;
	char problem[32];
	double vel; //in case of a given inflow or wall velocity

	//read the parameters, using problem.dat, including wl, wr, wt, wb
	read_parameters(filename, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax,
			&jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value, &wl, &wr, &wt, &wb, problem, &presLeft, &presRight, &presDelta, &vel);

	int pics = dt_value/dt; //just a helping variable for outputing vtk


	//allocate memory, including Flag
	U = matrix(0, imax+1, 0, jmax+1);
	V = matrix(0, imax+1, 0, jmax+1);
	P = matrix(0, imax+1, 0, jmax+1);
	RS = matrix(1, imax, 1, jmax);
	F = matrix(0, imax, 1, jmax);
	G = matrix(1, imax, 0, jmax);
	Flag = imatrix(0, imax+1, 0, jmax+1); // or Flag = imatrix(1, imax, 1, jmax);

	int kmax = 20; //test no slip boundary value function
	double ***U3d  = (double ***) malloc((size_t)((imax+1)*(jmax+1)*(kmax+1) * sizeof(double*)) ); //test no slip boundary value function
	double ***V3d  = (double ***) malloc((size_t)((imax+1)*(jmax+1)*(kmax+1) * sizeof(double*)) ); //test no slip boundary value function
	double ***W3d  = (double ***) malloc((size_t)((imax+1)*(jmax+1)*(kmax+1) * sizeof(double*)) ); //test no slip boundary value function
	int ***Flag3d  = (int ***) malloc((size_t)((imax+1)*(jmax+1)*(kmax+1) * sizeof(int*)) ); //test no slip boundary value function

	//initialisation, including **Flag
	init_flag(problem, imax, jmax, presDelta, Flag);
	init_uvp(UI, VI, PI, imax, jmax, U, V, P, problem);

	//going through all time steps
	while(t < t_end){
		//adaptive time stepping
		calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V);
		
		//setting bound.values
		boundaryvalues(imax, jmax, U, V, P, wl, wr, wt, wb, F, G, problem, Flag, vel); //including P, wl, wr, wt, wb, F, G, problem

                //test no slip boundary value function
                for(int i=1; i<=imax; i++){
                        for(int j=1; j<=jmax; j++){
                                for(int k=1; k<=kmax; k++){
                                        boundaryvalues_no_slip(i, j, k, U3d, V3d, W3d, Flag3d); //test no slip boundary value function
                                }
                        }
		}
		
		//computing F, G and right hand side of pressue eq.
		calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G, Flag);
		calculate_rs(dt, dx, dy, imax, jmax, F, G, RS);
		
		//iteration counter
		it = 0;
		
		do{
			//
			//perform SOR iteration, at same time set bound.values for P and new residual value
			sor(omg, dx, dy, imax, jmax, P, RS, &res, Flag, presLeft, presRight);

			it++;
		}while(it<itermax && res>eps);
/*		if (it == itermax) {
			printf("Warning: sor while loop exits because it reaches the itermax. res = %f, time =%f\n", res, t);
		}
*/		//calculate U and V of this time step
		calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P, Flag);
		
		//indent time and number of time steps
		n++;
		t += dt;
		
		//output of pics for animation
		if (n%pics==0 ){
			write_vtkFile(filename, n, xlength, ylength, imax, jmax, dx, dy, U, V, P);
		}
	}
	//output of U, V, P at the end for visualization
	//write_vtkFile("DrivenCavity", n, xlength, ylength, imax, jmax, dx, dy, U, V, P);

	//free memory
	free_matrix(U, 0, imax+1, 0, jmax+1);
	free_matrix(V, 0, imax+1, 0, jmax+1);
	free_matrix(P, 0, imax+1, 0, jmax+1);
	free_matrix(RS, 1, imax, 1, jmax);
	free_matrix(F, 0, imax, 1, jmax);
	free_matrix(G, 1, imax, 0, jmax);
	free_imatrix(Flag, 0, imax+1, 0, jmax+1);

	free(U3d);
	free(V3d);
	free(W3d);
	free(Flag3d);
	return -1;
}
