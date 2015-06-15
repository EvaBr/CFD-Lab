#include "parallel.h"

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


void initialiseBuffers ( double **sendBuffer, double **readBuffer, int *xlength, int *proc) {
	int kx = xlength[0]/proc[0] + 2; //size of buffer x
	int ky = xlength[1]/proc[1] + 2; //size of buffer y
	int kz = xlength[2]/proc[2] + 2; //size of buffer z

	//top and bottom buffers
	double *bt = calloc(kx*ky*5, sizeof(double));
	double *bbo = calloc(kx*ky*5, sizeof(double));

	//left and right buffers
	double *bl = calloc(kz*ky*5, sizeof(double));
	double *br = calloc(kz*ky*5, sizeof(double));

	//front and back buffers
	double *bf = calloc(kx*kz*5, sizeof(double));
	double *bba = calloc(kx*kz*5, sizeof(double));

	//gather all into the send buffer
	sendBuffer[0] = bl;
	sendBuffer[1] = br;
	sendBuffer[2] = bt;
	sendBuffer[3] = bbo;
	sendBuffer[4] = bf;
	sendBuffer[5] = bba;


	//read buffer only contains pointers to right send buffers
	readBuffer[0] = sendBuffer[1]; //left read buffer
	readBuffer[1] = sendBuffer[0]; //right read buffer
	readBuffer[2] = sendBuffer[3]; //top read buffer (reads from [index+xProc*yProc])
	readBuffer[3] = sendBuffer[2]; //bottom read buffer (reads from [index-xProc*yProc])
	readBuffer[4] = sendBuffer[5]; //front read buffer (reads from [index-xProc])
	readBuffer[5] = sendBuffer[4]; //back read buffer (reads from [index+xProc])
	//be careful - in the code, where you have explicit reading or writing of a buffer, change indices accordingly...
}
