#include "initLB.h"
#include "helper.h"
#include "LBDefinitions.h"
#include <mpi.h> // we need this right??

int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int *proc,  int argc, char *argv[]){
        // argument handling
        if (argc !=2 ) {
                printf("When running the simulation, please give a valid file name to read from!\n");
                return 1;
        }

	// get file name
	const char *FileName = NULL;
	FileName = argv[1];

	// read parameters, using read functions from helper
	read_int ( FileName, "xlengthx", &xlength[0] );
	read_int ( FileName, "xlengthy", &xlength[1] );
	read_int ( FileName, "xlengthz", &xlength[2] );

	read_double ( FileName, "tau", tau );

	read_double ( FileName, "velocityWallx", &velocityWall[0] );
	read_double ( FileName, "velocityWally", &velocityWall[1] );
	read_double ( FileName, "velocityWallz", &velocityWall[2] );

	read_int ( FileName, "timesteps", timesteps );
	read_int ( FileName, "timestepsPerPlotting", timestepsPerPlotting );

	read_int ( FileName, "iProc", &proc[0] );
	read_int ( FileName, "jProc", &proc[1] );
	read_int ( FileName, "kProc", &proc[2] );

	return 0;
}


void distributeParameters ( int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int *proc ){
	// broadcasting the read parameters to other ranks
	MPI_Bcast( xlength, 3, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( velocityWall, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
	MPI_Bcast( timesteps, 1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( timestepsPerPlotting, 1, MPI_INT, 0, MPI_COMM_WORLD );
	MPI_Bcast( proc, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );

	return 0;
}



void initialiseFields(double *collideField, double *streamField, int *flagField, int *subdomain, int rank, int *proc){
	// initialization of particle distribution func fields						TODO!!! change for xlength[3], add ghostlayer

	for (int a=0; a<subdomain[0]+2; a++){
		for (int b=0; b<subdomain[1]+2; b++){
			for (int c=0; c<subdomain[2]+2; c++){
				for (int i=0; i<Q; i++){
					// initialize streamField and collideField arrays
					streamField [ Q*compute_index(a, b, c, subdomain) + i ] = LATTICEWEIGHTS [i];
					collideField [ Q*compute_index(a, b, c, subdomain) + i ] = LATTICEWEIGHTS [i];
				}
				// set all as inner points
				flagField [ compute_index(a, b, c, subdomain) ] = FLUID;
			}
		}
	}


	// overwrite fluid to no-slip at other boundaries: this depends on the rank...
	//
	// 1. TOP BOUNDARY
	// if c==xlength+1, overwrite as moving wall:  now you have to check rank, not xlength, to know where you are - since
	// rank is growing with z coordinate (according to our numbering), c==xlength+1 arises at iproc*jproc*(kproc-1)<rank<=number_of_ranks
	if (rank >= proc[0]*proc[1]*(proc[2]-1)){ //plus, the ranks go from zero!
		for (int i=0; i<subdomain[0]+2; i++){
			for (int j=0; j<subdomain[1]+2; j++){
				flagField [ compute_index(i, j, subdomain[2]+1, subdomain)] = MOVING_WALL;
			}
		}
	} else { //in this case we are on the interior ghost boundary
		for (int i=1; i<subdomain[0]+1; i++){
			for (int j=1; j<subdomain[1]+1; j++){
				flagField [ compute_index(i, j, subdomain[2]+1, subdomain) ] = PARALLEL_BOUNDARY;
			}
		}
	}


	// 2. BOTTOM BOUNDARY
	if (rank < proc[0]*proc[1]){ //then our bottom boundary is outer; so no slip.
		for (int i=0; i<subdomain[0]+2; i++){
			for (int j=0; j<subdomain[1]+2; j++){
				flagField [ compute_index(i, j, 0, subdomain) ] = NO_SLIP;
			}
		}
	} else { //in this case we are on the interior ghost boundary
		for (int i=1; i<subdomain[0]+1; i++){
			for (int j=1; j<subdomain[1]+1; j++){
				flagField [ compute_index(i, j, 0, subdomain) ] = PARALLEL_BOUNDARY;
			}
		}
	}




	// 3. FRONT BOUNDARY
	if (rank%(proc[1]*proc[0]) < proc[0]) {
	//in this case, our front boundary is no slip:
		for (int i=0;  i<subdomain[0]+2; i++){
			for (int k=0; k<subdomain[2]+2; k++){
				flagField [ compute_index(i, 0, k, subdomain) ] = NO_SLIP; //y=0
			}
		}
	} else { //we are on the inside boundary
		for (int i=1; i<subdomain[0]+1; i++){
			for (int k=1; k<subdomain[2]+1; k++){
				flagField [ compute_index(i, 0, k, subdomain) ] = PARALLEL_BOUNDARY;
			}
		}
	}


	// 4. BACK BOUNDARY
	if (rank%(proc[1]*proc[0]) >= proc[0]*(proc[1]-1)) { //(y-1)*x<=rnk<y*x
	//in this case, our back boundary is no slip:
		for (int i=0;  i<subdomain[0]+2; i++){
			for (int k=0; k<subdomain[2]+2; k++){ 	//y=ymax
				flagField [compute_index(i, (subdomain[1]+1), k, subdomain) ] = NO_SLIP;
			}
		}
	} else { //we are on the inside boundary
		for (int i=1; i<subdomain[0]+1; i++){
			for (int k=1; k<subdomain[2]+1; k++){
				flagField [ compute_index(i, (subdomain[1]+1), k, subdomain) ] = PARALLEL_BOUNDARY;
			}
		}
	}




	// 5. LEFT BOUNDARY
	if (rank % proc[0] == 0){ //we're outer: no slip
		for (int k=0; k<subdomain[2]+2; k++){
			for (int j=0; j<subdomain[1]+2; j++){
				flagField [compute_index(0, j, k, subdomain) ] = NO_SLIP;
			}
		}
	} else { //in this case we are on the interior ghost boundary
		for (int k=1; k<subdomain[2]+1; k++){
			for (int j=1; j<subdomain[1]+1; j++){
				flagField [ compute_index(0, j, k, subdomain) ] = PARALLEL_BOUNDARY;
			}
		}
	}



	// 6. RIGHT BOUNDARY
	if (rank%proc[0] == proc[0]-1){ //we're outer: no slip
		for (int k=0; k<subdomain[2]+2; k++){
			for (int j=0; j<subdomain[1]+2; j++){
				flagField [ compute_index( subdomain[0]+1, j, k, subdomain) ] = NO_SLIP;
			}
		}
	} else { //in this case we are on the interior ghost boundary
		for (int k=1; k<subdomain[2]+1; k++){
			for (int j=1; j<subdomain[1]+1; j++){
				flagField [ compute_index( subdomain[0]+1, j, k, subdomain) ] = PARALLEL_BOUNDARY;
			}
		}
	}



	int ind;
	//Take care of the edges between  surfaces: the edges between two that are flagged as PARALLEL_BOUNDARY need yet to be set.
	//We will traverse through all these edges-points, just to avoid the if sentences... (see lower comment)
	for (int i=0; i<subdomain[0]+2; i++){  //front-bottom, back-bottom, front-top, back-top
				ind = compute_index(i, 0, 0, subdomain);
				flagField[ind] = max(PARALLEL_BOUNDARY, flagField[ind]); // we use max just to avoid if sentences... this is better, right?

				ind = compute_index(i, subdomain[1]+1, 0, subdomain);
				flagField[ind] = max(PARALLEL_BOUNDARY, flagField[ind]);

				ind = compute_index(i, 0, subdomain[2]+1, subdomain);
				flagField[ind] = max(PARALLEL_BOUNDARY, flagField[ind]);

				ind = compute_index(i, subdomain[1]+1, subdomain[2]+1, subdomain);
				flagField[ind] = max(PARALLEL_BOUNDARY, flagField[ind]);
	}
	for (int j=0; j<subdomain[1]+2; j++){ //left-bottom, right-bottom, left-top, right-top
				ind = compute_index(0, j, 0, subdomain);
				flagField[ind] = max(PARALLEL_BOUNDARY, flagField[ind]);

				ind = compute_index(subdomain[0]+1, j, 0, subdomain);
				flagField[ind] = max(PARALLEL_BOUNDARY, flagField[ind]);

				ind = compute_index(0, j, subdomain[2]+1, subdomain);
				flagField[ind] = max(PARALLEL_BOUNDARY, flagField[ind]);

				ind = compute_index(subdomain[0]+1, j, subdomain[2]+1, subdomain);
				flagField[ind] = max(PARALLEL_BOUNDARY, flagField[ind]);
	}
	for (int k=0; k<subdomain[2]+2; k++){ //front-left, front-right, back-left, back-right
				ind = compute_index(0, 0, k, subdomain);
				flagField[ind] = max(PARALLEL_BOUNDARY, flagField[ind]);

				ind = compute_index(subdomain[0]+1, 0, k, subdomain);
				flagField[ind] = max(PARALLEL_BOUNDARY, flagField[ind]);

				ind = compute_index(0, subdomain[1]+1, k, subdomain);
				flagField[ind] = max(PARALLEL_BOUNDARY, flagField[ind]);

				ind = compute_index(subdomain[0]+1, subdomain[1]+1, k, subdomain);
				flagField[ind] = max(PARALLEL_BOUNDARY, flagField[ind]);
	}

}
