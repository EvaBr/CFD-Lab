#ifndef _MAIN_C_
#define _MAIN_C_

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "LBDefinitions.h"

int main(int argc, char *argv[]){

	double *collideField = NULL;
	double *streamField = NULL;
	int *flagField = NULL;
	double *swap = NULL;

	int xlength=100;
	double tau=1.4;
	double velocityWall[3]={1.0, 0.0, 0.0};
	int timesteps=14;
//	int timestepsPerPlotting;


	/*//read parameters, exit if program call not valid
	int fail;
	fail = readParameters ( &xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, argc, argv );
	if (fail) {
		return 1;
	}*/

	int len = (xlength + 2)*(xlength + 2)*(xlength + 2);

	// initialize space for pointers
	collideField = calloc ( Q*len, sizeof(double) );
	streamField = calloc ( Q*len, sizeof(double) );
	flagField = calloc ( len, sizeof(int) );

	initialiseFields ( collideField, streamField, flagField, xlength );

	int t;
	for (t = 0; t < timesteps; t++){
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

  	return 0;
}

#endif

