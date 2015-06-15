#ifndef _MAIN_C_
#define _MAIN_C_

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "LBDefinitions.h"
//#include <mpi.h> //do not need since it is already included in parallel.h?
#include "materials/helper_functions/parallel.h"


int main(int argc, char *argv[]){

	double *collideField = NULL;
	double *streamField = NULL;
	int *flagField = NULL;
	double *swap = NULL;

	int xlength[3];
	double tau;
	double velocityWall[3];
	int timesteps;
	int timestepsPerPlotting;

	int rank;
	int number_of_ranks;

	int proc[3];

	double *sendBuffer[6];
	double *readBuffer[6];

	//initialize MPI run
	//initializeMPI (&rank, &number_of_ranks, argc, argv );   ????
	MPI_Init ( &argc, &argv);
	MPI_Comm_size ( MPI_COMM_WORLD, &number_of_ranks);
	MPI_Comm_rank ( MPI_COMM_WORLD, &rank);

	//first processor reads parameters, exits if program call not valid
	if (rank == 0) {
		//announce the start of simulation
		printf("LBM Simulation of Lid Driven Cavity (by CFD Lab group 9)\n");
		printf("----------------------------------------------------------\n");
		printf("Reading the parameters... \n");

		int fail;
		fail = readParameters ( xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, proc, argc, argv );
		//check if parameters could be read
		if (fail) {
			printf("The parameter file could not be read. Aborting the program. \n");
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}
		//check if read parameters were adequately chosen
		if (proc[0]*proc[1]*proc[2] != number_of_ranks) {
			printf("Number of available processors differs from the one in parameter file!\n Aborting the program. \n");
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}
		//do the same for the case when xlen not divisible by num.of procesors 
		if (xlength[0]%proc[0]!=0 || xlength[1]%proc[1]!=0 || xlength[2]%proc[2]!=0) {
			printf("Size of domain is not divisible by number of processes. \n Aborting the program. \n");
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}


		//if we got to here, parameters are read and ok. first: distribute them to all other procesors
		printf("Done reading.\n Simulation started (with MPI, using %i ranks).\n", number_of_ranks);

		//distribution of parameters: TODO
		distributeParameters ( &xlength, &tau, &velocityWall, &timesteps, &timestepsPerPlotting, &proc );
	}


	//calculate dimensions of subdomain, dealt with by one process
	int subdomain[3] = { xlength[0] / proc[0], xlength[1] / proc[1], xlength[2] / proc[2] };
	//calculate subdomain volume
	int vol = (subdomain[0] + 2)*(subdomain[1] + 2)*(subdomain[2] + 2);

	// initialize space for pointers
	collideField = calloc ( Q*vol, sizeof(double) );
	streamField = calloc ( Q*vol, sizeof(double) );
	flagField = calloc ( vol, sizeof(int) );

	//initialise fields for this subdomain
	initialiseFields ( collideField, streamField, flagField, subdomain, rank, proc );

	//initialize buffers  TODO
	initialiseBuffers(sendBuffer, readBuffer, subdomain);

	int t;
	for (t = 0; t < timesteps; t++){
		//extraction, swap, injection for x
		if ( flagField[  index(0, subdomain[1]/2, subdomain[2]/2) ] ==PARALLEL_BOUNDARY) {
			extraction(...);
			swap(..); //copy our send buffer to neighbour's read buffer, and copy neighbour's send buffer to our read buffer.
			//USE Sendrecv!
			injection(..);
		}
		//extraction, swap, injection for y
		//TODO
		//extraction, swap, injection for z
		//TODO



		doStreaming ( collideField, streamField, flagField, subdomain );

		// swap the stream and collide arrays
		swap = collideField;
		collideField = streamField;
		streamField = swap;

		doCollision ( collideField, flagField, &tau, subdomain);

		treatBoundary ( collideField, flagField, velocityWall, subdomain);

		// write vtk data   TODO!
		/*if (t%timestepsPerPlotting==0){
			writeVtkOutput ( collideField, flagField, "DrivenCavity", t, xlength );
		}*/
	}


	//we finished with the simulation and output
	if (rank==0){
		printf("Simulation finished. Visualisation data written. \n Freeing allocated memory...\n");
	}

	// free the initialized space
	free ( collideField );
	free ( streamField );
	free ( flagField );
	//(dont forget the buffers)
	for (int i = 0; i < 6; ++i) {
		free ( sendBuffer[i] );
		free ( readBuffer[i] ); //a je treba tegale tko, če daš samo kazalce notr?
	}


	//program is done.
	if (rank==0){
		printf("...done.  :) \n\n");
	}


	//synchronization of processes and end of the MPI session
	//finalizeMPI(); 				 ????
	Programm_Stop("op_log.txt");

  	return 0;
}

#endif

