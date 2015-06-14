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
	//initializeMPI (&rank, &number_of_ranks, argc, argv );   ?????
	MPI_Init ( &argc, &argv);
	MPI_Comm_size ( MPI_COMM_WORLD, &number_of_ranks);
	MPI_Comm_rank ( MPI_COMM_WORLD, &rank);

	//read parameters, exit if program call not valid
	int fail;
	fail = readParameters ( xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, proc, argc, argv );
	if (fail) {
		return 1;
	}

	int len = (xlength[0] + 2*proc[0])*(xlength[1] + 2*proc[1])*(xlength[2] + 2*proc[2]);

	// initialize space for pointers
	collideField = calloc ( Q*len, sizeof(double) );
	streamField = calloc ( Q*len, sizeof(double) );
	flagField = calloc ( len, sizeof(int) );

	initialiseFields ( collideField, streamField, flagField, xlength, rank, number_of_ranks );

	//initialize buffers
	initialiseBuffers(sendBuffer, readBuffer, xlength, proc);

	int t;
	for (t = 0; t < timesteps; t++){
		//extraction, swap, injection for x
		//TODO
		//extraction, swap, injection for y
		//TODO
		//extraction, swap, injection for z
		//TODO



		doStreaming ( collideField, streamField, flagField, xlength );

		// swap the stream and collide arrays
		swap = collideField;
		collideField = streamField;
		streamField = swap;

		doCollision ( collideField, flagField, &tau, xlength);

		treatBoundary ( collideField, flagField, velocityWall, xlength);

		// write vtk data
	//	if (t%timestepsPerPlotting==0){
	//		writeVtkOutput ( collideField, flagField, "DrivenCavity", t, xlength );
	//	}
	}
	// free the initialized space
	free ( collideField );
	free ( streamField );
	free ( flagField );


	//synchronization of processes and end of the MPI session
	//finalizeMPI(); 				 ????
	Programm_Stop("op_log.txt");

  	return 0;
}

#endif

