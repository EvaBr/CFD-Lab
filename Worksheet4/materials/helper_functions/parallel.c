#include "parallel.h"
#include <mpi.h>

#define N_BORDER_VEL 5

void Program_Message(char *txt)
/* produces a stderr text output  */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
}


void Programm_Sync(char *txt)
/* produces a stderr textoutput and synchronize all processes */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */  
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
}


void Programm_Stop(char *txt)
/* all processes will produce a text output, be synchronized and finished */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                           /* synchronize output */
   fprintf(stderr,"-STOP- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(1);
}


void initialiseBuffers ( double **sendBuffer, double **readBuffer, int *subdomain) {
	int kx = subdomain[0] + 2; //size of buffer x
	int ky = subdomain[1] + 2; //size of buffer y
	int kz = subdomain[2] + 2; //size of buffer z


	//left and right buffers
	sendBuffer[0] = calloc(kz*ky*5, sizeof(double));
	sendBuffer[1] = calloc(kz*ky*5, sizeof(double));

	//top and bottom buffers
	sendBuffer[2] = calloc(kx*ky*5, sizeof(double));
	sendBuffer[3] = calloc(kx*ky*5, sizeof(double));

	//front and back buffers
	sendBuffer[4] = calloc(kx*kz*5, sizeof(double));
	sendBuffer[5] = calloc(kx*kz*5, sizeof(double));


	//read buffers - analogously
	readBuffer[0] = calloc(kz*ky*5, sizeof(double));
	readBuffer[1] = calloc(kz*ky*5, sizeof(double));
	sendBuffer[2] = calloc(kx*ky*5, sizeof(double));
	sendBuffer[3] = calloc(kx*ky*5, sizeof(double));
	sendBuffer[4] = calloc(kx*kz*5, sizeof(double));
	sendBuffer[5] = calloc(kx*kz*5, sizeof(double));
}

void extractionX ( double **sendBuffer, int rank, int *subdomain){
	//here we extract the five pdfs that would be streamed into our neighbour, x direction.  /just a comment - should extract left and right be their own functions?
	// 1. EXTRACT TO LEFT
	// pdfs, that need to be extracted are the ones with indices 1, 5, 8, 11, 15.
	int i = 0;
	if (rank%subdomain[0]!=0){ //check if left boundary is not a no-slip
		for (int z=0; z<subdomain[2]+2; z++){
			for (int y=0; y<subdomain[1]+2; y++){
				*(sendBuffer[0] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 1 ];
				*(sendBuffer[0] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 5 ];
				*(sendBuffer[0] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 8 ];
				*(sendBuffer[0] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 11 ];
				*(sendBuffer[0] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 15 ];
			}
		}
	}

	// 2. EXTRACT TO RIGHT
	// here, we extract pdfs with indeces 3, 7, 10, 13, 17.
	if (rank%subdomain[0]!=subdomain[0]-1){ //check if right boundary is not a no-slip
		i = 0;
		for (int z=0; z<subdomain[2]+2; z++){
			for (int z=0; z<subdomain[2]+2; z++){
				*(sendBuffer[0] + i++) = collideField [ Q*compute_index(subdomain[0], y, z, subdomain) + 3 ];
				*(sendBuffer[0] + i++) = collideField [ Q*compute_index(subdomain[0], y, z, subdomain) + 7 ];
				*(sendBuffer[0] + i++) = collideField [ Q*compute_index(subdomain[0], y, z, subdomain) + 10 ];
				*(sendBuffer[0] + i++) = collideField [ Q*compute_index(subdomain[0], y, z, subdomain) + 13 ];
				*(sendBuffer[0] + i++) = collideField [ Q*compute_index(subdomain[0], y, z, subdomain) + 17 ];
			}
		}
	}
}
void extractionY ( ){
	//here we extract the five pdfs that would be streamed into our neighbour, y direction.
	//TODO
}
void extractionZ ( ){
	//here we extract the five pdfs that would be streamed into our neighbour, z direction.
	//TODO
}

//Injecion
//here we inject the five pdfs that would be streamed into our neighbour, x direction.
void injectionXright ( double **readBuffer, int rank, int *subdomain){
	// 1. INJECT FROM RIGHT
	// pdfs, that need to be extracted are the ones with indices 1, 5, 8, 11, 15.
	int i = 0;
	if (rank%subdomain[0]!=subdomain[0]-1){ //check if right boundary is not a no-slip
		for (int z=0; z<subdomain[2]+2; z++){
			for (int y=0; y<subdomain[1]+2; y++){
				collideField [ Q*compute_index(1, y, z, subdomain) + 1 ] = *(readBuffer[1] + i++);
				collideField [ Q*compute_index(1, y, z, subdomain) + 5 ] = *(readBuffer[1] + i++);
				collideField [ Q*compute_index(1, y, z, subdomain) + 8 ] = *(readBuffer[1] + i++);
				collideField [ Q*compute_index(1, y, z, subdomain) + 11 ] = *(readBuffer[1] + i++);
				collideField [ Q*compute_index(1, y, z, subdomain) + 15 ] = *(readBuffer[1] + i++);
			}
		}
	}
}

void injectionXleft ( double **readBuffer, int rank, int *subdomain){
	// 2. INJECT FROM LEFT
	// here, we extract pdfs with indices 3, 7, 10, 13, 17.
	int i = 0;
	if (rank%subdomain[0]!=0){ //check if left boundary is not a no-slip
		for (int z=0; z<subdomain[2]+2; z++){
			for (int y=0; y<subdomain[1]+2; y++){
				collideField [ Q*compute_index(1, y, z, subdomain) + 3 ] = *(readBuffer[1] + i++);
				collideField [ Q*compute_index(1, y, z, subdomain) + 7 ] = *(readBuffer[1] + i++);
				collideField [ Q*compute_index(1, y, z, subdomain) + 10 ] = *(readBuffer[1] + i++);
				collideField [ Q*compute_index(1, y, z, subdomain) + 13 ] = *(readBuffer[1] + i++);
				collideField [ Q*compute_index(1, y, z, subdomain) + 17 ] = *(readBuffer[1] + i++);
			}
		}
	}
}
