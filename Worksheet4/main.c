#ifndef _MAIN_C_
#define _MAIN_C_

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "LBDefinitions.h"
#include "materials/helper_functions/parallel.h"
#include "helper.h"

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
	//initializeMPI (&rank, &number_of_ranks, argc, argv );
	//#############################
	MPI_Init ( &argc, &argv);
	MPI_Comm_size ( MPI_COMM_WORLD, &number_of_ranks);
	MPI_Comm_rank ( MPI_COMM_WORLD, &rank);
	//############################


	//first processor reads parameters, exits if program call not valid
	if (rank == 0) {

		//announce the start of simulation
		printf("\nLBM Simulation of Lid Driven Cavity (by CFD Lab group 9)\n");
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

	}

	//distribution of parameters to other ranks:
	distributeParameters ( xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, proc );

	//calculate dimensions of subdomain, dealt with by one process
	int subdomain[3] = { xlength[0] / proc[0], xlength[1] / proc[1], xlength[2] / proc[2] };
	//calculate subdomain volume
	int vol = (subdomain[0] + 2)*(subdomain[1] + 2)*(subdomain[2] + 2);

	// initialize space for pointers
	collideField = calloc ( Q*vol, sizeof(double) );
	streamField = calloc ( Q*vol, sizeof(double) );
	flagField = calloc ( vol, sizeof(int) );

	//initialise fields for this subdomain
	initialiseFields ( collideField, streamField, flagField, subdomain, rank, proc, sendBuffer, readBuffer);

	//initialize buffers - this is done in initialise fields
//	initialiseBuffers ( sendBuffer, readBuffer, subdomain, flagField );



// Program_Message("after initialisations");



	int t;
	for (t = 0; t < timesteps; t++){

	//########################################################################################
		//extraction, swap, injection for x
		// 1. LEFT; check, that rank doesn't have a no-slip on the left
		if (rank%proc[0]!=0) {
			extractionXleft ( sendBuffer, collideField, subdomain );
			swapXleft ( sendBuffer, readBuffer, subdomain, proc, rank); //copy our send buffer to neighbour's read buffer, and vice-versa
			injectionXleft ( readBuffer, collideField, subdomain );
		}
		// 2. RIGHT; check, that rank doesn't have a no-slip on the right
		if (rank%proc[0]!=proc[0]-1){
			extractionXright ( sendBuffer, collideField, subdomain );
			swapXright ( sendBuffer, readBuffer, subdomain, proc, rank);
			injectionXright ( readBuffer, collideField, subdomain );
		}


		//extraction, swap, injection for y
		// 3. FRONT; check, that rank doesn't have a no-slip at the front
		if (rank%(proc[0]*proc[1])>=proc[0]){
			extractionYfront ( sendBuffer, collideField, subdomain );
			swapYfront ( sendBuffer, readBuffer, subdomain, proc, rank );
			injectionYfront ( readBuffer, collideField, subdomain );
		}
		// 4. BACK; check, that rank doesn't have a no-slip  at the back
		if (rank%(proc[0]*proc[1])<proc[0]*(proc[1]-1)){
			extractionYback ( sendBuffer, collideField, subdomain );
			swapYback ( sendBuffer, readBuffer, subdomain, proc, rank );
			injectionYback ( readBuffer, collideField, subdomain );
		}


		//extraction, swap, injection for z
		// 5. TOP; check, that rank doesn't have a no-slip at the top
		if (rank<(proc[0]*proc[1]*(proc[2]-1))){
			extractionZtop ( sendBuffer, collideField, subdomain );
			swapZtop ( sendBuffer, readBuffer, subdomain, proc, rank );
			injectionZtop ( readBuffer, collideField, subdomain );
		}
		// 6. BOTTOM; check, that rank doesn't have a no-slip at the bottom
		if (rank >= proc[0]*proc[1]){
			extractionZbottom ( sendBuffer, collideField, subdomain );
			swapZbottom ( sendBuffer, readBuffer, subdomain, proc, rank );
			injectionZbottom ( readBuffer, collideField, subdomain );
		}
	//#########################################################################################


		// do the actual streaming step
		doStreaming ( collideField, streamField, flagField, subdomain );

		// swap the stream and collide arrays
		swap = collideField;
		collideField = streamField;
		streamField = swap;

		doCollision ( collideField, flagField, &tau, subdomain);

		treatBoundary ( collideField, flagField, velocityWall, subdomain);

		// write vtk (partial) data
		if (t%timestepsPerPlotting==0){
			writeVtkOutput ( collideField, flagField, "DrivenCavity", rank, t, subdomain, xlength, proc );
		}
	}



	//we finished with the simulation and output
	if (rank==0){
		printf("Simulation finished. Visualisation data written. \n Freeing allocated memory...\n");
		/*for (int i =0; i<6; i++){
			for (int j = 0; j<subdomain[0]*subdomain[1]; j++){
				printf("sendBuff [%d,%d] =  %f\n", i,j,sendBuffer[i][j]);
			}
		}*/
	}



	MPI_Barrier(MPI_COMM_WORLD);  // to make sure we don't free too quickly

	// free the initialized space
	free ( collideField );
	free ( streamField );
	free ( flagField );
	// (dont forget the buffers)
/*	for (int i = 0; i < 6; ++i) {
                sendBuffer[5] = calloc(kx*kz*5, sizeof(double));
                readBuffer[5] = calloc(kx*kz*5, sizeof(double));
	}*/
	int kx = subdomain[0]+2;
	int ky = subdomain[1]+2;
	int kz = subdomain[2]+2;
 	if (flagField[compute_index(0, ky/2, kz/2, subdomain)]==PARALLEL_BOUNDARY){
                //left buffers
		free ( sendBuffer[0] );
		free ( readBuffer[0] );
        }

        if (flagField[compute_index(kx-1, ky/2, kz/2, subdomain)]==PARALLEL_BOUNDARY){
                // right buffers
		free ( sendBuffer[1] );
		free ( readBuffer[1] );
        }

        if (flagField[compute_index(kx/2, ky/2, kz-1, subdomain)]==PARALLEL_BOUNDARY){
                //top  buffers
		free ( sendBuffer[2] );
		free ( readBuffer[2] );
        }

        if (flagField[compute_index(kx/2, ky/2, 0, subdomain)]==PARALLEL_BOUNDARY){
                //bottom buffers
		free ( sendBuffer[3] );
		free ( readBuffer[3] );
        }

        if (flagField[compute_index(kx/2, 0, kz/2, subdomain)]==PARALLEL_BOUNDARY){
		//front buffers
		free ( sendBuffer[4] );
		free ( readBuffer[4] );
        }

        if (flagField[compute_index(kx/2, ky-1, kz/2, subdomain)]==PARALLEL_BOUNDARY){
                //back buffers
		free ( sendBuffer[5] );
		free ( readBuffer[5] );
        }



	// program is done.
	if (rank==0){
		printf("...done.  :) \n\n");
	}


	// synchronization of processes and end of the MPI session
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	/*Programm_Stop("done");*/

  	return 0;
}

#endif

