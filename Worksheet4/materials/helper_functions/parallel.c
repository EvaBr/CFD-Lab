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



//	E X T R A C T I O N     F U N C T I O N S

// Here we extract the five pdfs that would be streamed into our neighbour, x direction.
void extractionX ( double **sendBuffer, double *collideField, int *subdomain){
	// 1. EXTRACT TO LEFT
	// pdfs, that need to be extracted are the ones with indices 1, 5, 8, 11, 15.
	int i = 0;
/* before calling it, we need to check that front boundary is not a no-slip: if (rank%proc[0]!=0){ */
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
void extractionXright (double **sendBuffer, double  *collideField, int *subdomain){
	// 2. EXTRACT TO RIGHT
	// here, we extract pdfs with indeces 3, 7, 10, 13, 17.
	int i = 0;
/* before calling it, we need to check that front boundary is not a no-slip: if (rank%proc[0]!=proc[0]-1){ */
	for (int z=0; z<subdomain[2]+2; z++){
		for (int y=0; y<subdomain[1]+2; y++){
			*(sendBuffer[1] + i++) = collideField [ Q*compute_index(subdomain[0], y, z, subdomain) + 3 ];
			*(sendBuffer[1] + i++) = collideField [ Q*compute_index(subdomain[0], y, z, subdomain) + 7 ];
			*(sendBuffer[1] + i++) = collideField [ Q*compute_index(subdomain[0], y, z, subdomain) + 10 ];
			*(sendBuffer[1] + i++) = collideField [ Q*compute_index(subdomain[0], y, z, subdomain) + 13 ];
			*(sendBuffer[1] + i++) = collideField [ Q*compute_index(subdomain[0], y, z, subdomain) + 17 ];
		}
	}
}




// Here we extract the five pdfs that would be streamed into our neighbour, y direction.
void extractionYfront (double **sendBuffer, double *collideField, int *subdomain){
	// 1. EXTRACT TO FRONT
	// pdfs, that need to be extracted are the ones with indices 0, 5, 6, 7, 14.
	int i = 0;
/* before calling it, we need to check that front boundary is not a no-slip: if (rank%proc[0]*proc[1]>=proc[0]){  */
	for (int z=0; z<subdomain[2]+2; z++){
		for (int x=0; x<subdomain[0]+2; x++){
			*(sendBuffer[4] + i++) = collideField [ Q*compute_index(1, y, z, subdomain)  ];
			*(sendBuffer[4] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 5 ];
			*(sendBuffer[4] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 6 ];
			*(sendBuffer[4] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 7 ];
			*(sendBuffer[4] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 14 ];
		}
	}
}
void extractionYback (double **sendBuffer, double *collideField, int *subdomain){
	// 2. EXTRACT TO BACK
	// we need to extract pdfs with indices 4, 11, 12, 13, 18.
	int i = 0;
/* before calling it, we need to check that front boundary is not a no-slip: if (rank%proc[0]*proc[1]<proc[0]*(proc[1]-1)){  */
	for (int z=0; z<subdomain[2]+2; z++){
		for (int x=0; x<subdomain[0]+2; x++){
			*(sendBuffer[5] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 4 ];
			*(sendBuffer[5] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 11 ];
			*(sendBuffer[5] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 12 ];
			*(sendBuffer[5] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 13 ];
			*(sendBuffer[5] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 18 ];
		}
	}
}



// Here we extract the five pdfs that would be streamed into our neighbour, z direction.
void extractionZtop (double **sendBuffer, double *collideField, int *subdomain){
	// 1. EXTRACT TO TOP
	// pdfs, that need to be extracted are the ones with indices 14, 15, 16, 17, 18.
	int i = 0;
/* before calling it, we need to check that front boundary is not a no-slip: if (rank>=proc[0]*proc[1]*(proc[2]-1)){ */
	for (int y=0; y<subdomain[1]+2; y++){
		for (int x=0; x<subdomain[0]+2; x++){
			*(sendBuffer[2] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 14 ];
			*(sendBuffer[2] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 15 ];
			*(sendBuffer[2] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 16 ];
			*(sendBuffer[2] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 17 ];
			*(sendBuffer[2] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 18 ];
		}
	}
}
void extractionZbottom (double **sendBuffer, double *collideField, int *subdomain){
	// 2. EXTRACT TO BOTTOM
	// pdfs, that need to be extracted are the ones with indices 0, 1, 2, 3, 4.
	int i = 0;
/* before calling it, we need to check that front boundary is not a no-slip: if (rank<proc[0]*proc[1]){ */
	for (int y=0; y<subdomain[1]+2; y++){
		for (int x=0; x<subdomain[0]+2; x++){
			*(sendBuffer[3] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) ];
			*(sendBuffer[3] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 1 ];
			*(sendBuffer[3] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 2 ];
			*(sendBuffer[3] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 3 ];
			*(sendBuffer[3] + i++) = collideField [ Q*compute_index(1, y, z, subdomain) + 4 ];
		}
	}
}
