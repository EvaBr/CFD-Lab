#include "helper.h"
#include "visual.h"
#include "init.h"
#include "sor.h"
#include "boundary_val.h"
#include "uvp.h"
#include "surface.h"
#include <stdio.h>




int main(int argn, char** args){
	if (argn !=2 ) {
		printf("When running the simulation, please give a valid geometry file name!\n");
		return 1;
	}


	//set the geometry file
	char *filename = NULL;
	filename = args[1];

	//initialize variables
	double t = 0; /*time start*/
	int it, n = 0; /*iteration and time step counter*/
	double res; /*residual for SOR*/
	/*arrays*/
	double ***U, ***V, ***W, ***P;
	double ***RS, ***F, ***G, ***H;

	int ***Flag; //additional data structure for arbitrary geometry
	int ***S; //additional data structure for arbitrary geometry
	int matrix_output = 0;
	/*those to be read in from the input file*/
	double Re, UI, VI, WI, PI, GX, GY, GZ, t_end, xlength, ylength, zlength, dt, dx, dy, dz, alpha, omg, tau, eps, dt_value;
	int  imax, jmax, kmax, itermax;

	//double presLeft, presRight, presDelta; //for pressure stuff  ...TODO: not allowed for now
	int wl, wr, wt, wb, wf, wh;
	char problemGeometry[200];
	//in case of a given inflow or wall velocity  TODO: will we have this? needs to be a vector?
	double velIN;
	double velMW[3]; // the moving wall velocity is a vector

	struct particleline *Partlines = (struct particleline *)malloc((unsigned)(1 * sizeof(struct particleline)));
	//char szFileName[80];

	//read the parameters, using problem.dat, including wl, wr, wt, wb
	read_parameters(filename, &Re, &UI, &VI, &WI, &PI, &GX, &GY, &GZ, &t_end, &xlength, &ylength, &zlength, &dt, &dx, &dy, &dz, &imax,
			&jmax, &kmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value, &wl, &wr,  &wf, &wh, &wt, &wb, problemGeometry, &velIN, &velMW[0]); //&presLeft, &presRight, &presDelta, &vel);


	//printf("d: %f, %f, %f\n", dx,dy,dz);


	//int pics = dt_value/dt; //just a helping variable for outputing vtk
	double last_output_t = -dt_value;
	//double every = t_end/10; //helping variable. Can be used for displaying info every tenth of the progress during simulation

	//allocate memory, including Flag
	U = matrix2(0, imax+1, 0, jmax+1, 0, kmax+1);
	V = matrix2(0, imax+1, 0, jmax+1, 0, kmax+1);
	W = matrix2(0, imax+1, 0, jmax+1, 0, kmax+1);
	P = matrix2(0, imax+1, 0, jmax+1, 0, kmax+1);
	RS = matrix2(1, imax, 1, jmax, 1, kmax);

	F = matrix2(0, imax, 1, jmax, 1, kmax);
	G = matrix2(1, imax, 0, jmax, 1, kmax);
	H = matrix2(1, imax, 1, jmax, 0, kmax);


	Flag = imatrix2(0, imax+1, 0, jmax+1, 0, kmax+1); // or Flag = imatrix(1, imax, 1, jmax);

	S    = imatrix2(0, imax+1, 0, jmax+1, 0, kmax+1);
	//S = Flag -> adjust C_x Flags in helper.h

	//initialisation, including **Flag


	init_flag(problemGeometry, imax, jmax, kmax, Flag, wl, wr, wf, wh, wt, wb);  //presDelta, Flag);

	init_particles(Flag,dx,dy,dz,imax,jmax,kmax,3,Partlines);

	printf("!!!!!!!!%d\n",binMatch(B_NO,B_N));


	init_uvwp(UI, VI, WI, PI, Flag,imax, jmax, kmax, U, V, W, P, problemGeometry);

	write_particles("particles_init",0, 1, Partlines);

	write_imatrix2("Flag_start.txt",t,Flag, 0, imax, 1, jmax, 1, kmax);
	//write_flag_imatrix("Flag.txt",t,Flag, 0, imax+1, 0, jmax+1, 0, kmax+1);
	write_vtkFile(filename, -1, xlength, ylength, zlength, imax, jmax, kmax, dx, dy, dz, U, V, W, P,Flag);

	//write_vtkFile("init", n, xlength, ylength, zlength, imax, jmax, kmax, dx, dy, dz, U, V, W, P);
	//going through all time steps
	/*
	write_matrix2("P_start.txt",t,P, 0, imax+1, 0, jmax+1, 0, kmax+1);
	write_matrix2("U_start.txt",t,U, 0, imax+1, 0, jmax+1, 0, kmax+1);
	write_matrix2("V_start.txt",t,V, 0, imax+1, 0, jmax+1, 0, kmax+1);
	write_matrix2("W_start.txt",t,W, 0, imax+1, 0, jmax+1, 0, kmax+1);
	*/
	//setting bound.values
	boundaryvalues(imax, jmax, kmax, U, V, W, P, F, G, H, problemGeometry, Flag, velIN, velMW); //including P, wl, wr, wt, wb, F, G, problem
	//    printf("calc bc \n");


	while(t < t_end){
		/*if(t - every >= 0){
      	  	  printf("Calculating time %f ... \n", t);
      	  	  every += t;
    	}*/

		//adaptive time stepping
		calculate_dt(Re, tau, &dt, dx, dy, dz, imax, jmax, kmax, U, V, W);
		//    printf("calc dt \n");

		/*
		mark_cells(Flag,dx,dy,dz,imax,jmax,kmax,1,Partlines);


		set_uvwp_surface(U,V,W,P,Flag,dx,dy,dz,imax,jmax,kmax,GX,GY,GZ,dt,Re);
*/

		//computing F, G, H and right hand side of pressue eq.
		calculate_fgh(Re, GX, GY, GZ, alpha, dt, dx, dy, dz, imax, jmax, kmax, U, V, W, F, G, H, Flag);
		//    printf("calc fgh \n");

		calculate_rs(dt, dx, dy, dz, imax, jmax, kmax, F, G, H, RS,Flag);
		//    printf("calc rs \n");

		//iteration counter
		it = 0;

		//write_matrix2("RS.txt",t,RS, 1, imax, 1, jmax, 1, kmax);
		//write_matrix2("F.txt",t,F, 0, imax, 1, jmax, 1, kmax);
		//write_matrix2("G.txt",t,G, 1, imax, 0, jmax, 1, kmax);
		//write_matrix2("H.txt",t,H, 1, imax, 1, jmax, 0, kmax);
		do{
			// sprintf( szFileName, "P_%d.txt",it );
			//write_matrix2(szFileName,1000+t,P, 0, imax+1, 0, jmax+1, 0, kmax+1);
			//perform SOR iteration, at same time set bound.values for P and new residual value
			sor(omg, dx, dy, dz, imax, jmax, kmax, P, RS, &res, Flag); //, presLeft, presRight);

			it++;
		}while(it<itermax && res>eps);

		/*if (it == itermax) {
	printf("Warning: sor while loop exits because it reaches the itermax. res = %f, time =%f\n", res, t);
		} */

		//write_matrix2("P.txt",t,P, 0, imax+1, 0, jmax+1, 0, kmax+1);
		//calculate U, V and W of this time step
		calculate_uvw(dt, dx, dy, dz, imax, jmax, kmax, U, V, W, F, G, H, P, Flag);
		//write_matrix2("U.txt",t,U, 0, imax+1, 0, jmax+1, 0, kmax+1);
		//write_matrix2("V.txt",t,V, 0, imax+1, 0, jmax+1, 0, kmax+1);
		//write_matrix2("W.txt",t,W, 0, imax+1, 0, jmax+1, 0, kmax+1);

		//    printf("calc uvw \n");


			//sprintf( szFileName, "simulation/%s.%i_debug.vtk", szProblem, timeStepNumber );


		//setting bound.values
		boundaryvalues(imax, jmax, kmax, U, V, W, P, F, G, H, problemGeometry, Flag, velIN, velMW); //including P, wl, wr, wt, wb, F, G, problem

		/*
		set_uvwp_surface(U,V,W,P,Flag,dx,dy,dz,imax,jmax,kmax,GX,GY,GZ,dt,Re);

		printf("advance_particles!\n");
		advance_particles(dx,dy, dz,imax,jmax,kmax, dt,U,V,W,1,Partlines);
*/

		//indent time and number of time steps
		printf("timer\n");
		n++;
		t += dt;

		//output of pics for animation
		if ( t-last_output_t  >= dt_value ){  //n%pics==0 ){
			printf("output\n!");









			write_particles(filename,n, 1, Partlines);
			write_vtkFile(filename, n, xlength, ylength, zlength, imax, jmax, kmax, dx, dy, dz, U, V, W, P,Flag);

			printf("output vtk (%d)\n",n);
			last_output_t = t;
			matrix_output++;
		}
		printf("timestep: %f   - next output: %f   (dt: %f) \n",t,dt_value- (t-last_output_t),dt);
	}

	//output of U, V, P at the end for visualization
	//write_vtkFile("DrivenCavity", n, xlength, ylength, imax, jmax, dx, dy, U, V, P);
	//free memory
	free_matrix2(U, 0, imax+1, 0, jmax+1, 0, kmax+1);
	free_matrix2(V, 0, imax+1, 0, jmax+1, 0, kmax+1);
	free_matrix2(W, 0, imax+1, 0, jmax+1, 0, kmax+1);
	free_matrix2(P, 0, imax+1, 0, jmax+1, 0, kmax+1);
	free_matrix2(RS, 1, imax, 1, jmax, 1, kmax);


	free_matrix2(F, 0, imax, 1, jmax, 1, kmax);
	free_matrix2(G, 1, imax, 0, jmax, 1, kmax);
	free_matrix2(H, 1, imax, 1, jmax,  0, kmax);


	free_imatrix2(Flag, 0, imax+1, 0, jmax+1, 0, kmax+1);
	free_imatrix2(S, 0, imax+1, 0, jmax+1, 0, kmax+1);
	//printf("\n-\n");
	return -1;
}
