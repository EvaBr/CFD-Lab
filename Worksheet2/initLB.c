#include "initLB.h"
#include "helper.h"
#include "LBDefinitions.h"

int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){
        // argument handling
        if (argc !=2 ) {
                printf("When running the simulation, please give a valid file name to read from!\n");
                return 1;
        }

	const char *FileName = NULL;
	FileName = argv[1];

	read_int ( FileName, "xlength", xlength );

	read_double ( FileName, "tau", tau );

	read_double ( FileName, "velocityWallx", &velocityWall[0] );
	read_double ( FileName, "velocityWally", &velocityWall[1] );
	read_double ( FileName, "velocityWallz", &velocityWall[2] );

	read_int ( FileName, "timesteps", timesteps );
	read_int ( FileName, "timestepsPerPlotting", timestepsPerPlotting );

	return 0;
}


void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength){
	//initialization of particle distribution func fields
	int xlen = xlength + 2;
	int xlen2 = xlen*xlen;

	for (int a=0; a<=xlength+1; a++){
		for (int b=0; b<=xlength+1; b++){
			for (int c=0; c<=xlength+1; c++){
				for (int i=0; i<Q; i++){
					streamField [ Q*(c*xlen2 + b*xlen + a) + i ] = LATTICEWEIGHTS [i];
					collideField [ Q*(c*xlen2 + b*xlen + a) + i ] = LATTICEWEIGHTS [i];
				}
				//set all as inner points
				flagField [ c*xlen2 + b*xlen + a ] = FLUID;
			}
			//if c==xlength+1, moving wall:
			flagField [ xlen2*(xlength+1) + b*xlength + a ] = MOVING_WALL;
		}
	}
	//correct fluid to no slip boundaries:
	for (int k=0; k<xlen; k++){
		for (int j=0; j<xlength+1; j++){
			// if a||b||c==0 :
			flagField [ j*xlen2 + k*xlen ] = NO_SLIP;
			flagField [ j*xlen2 + k] = NO_SLIP;
			flagField [ k*xlen +j ] = NO_SLIP;
			//if a||b==xlength+1 :
			flagField [ j*xlen2 + k*xlen + xlen - 1 ] = NO_SLIP;
			flagField [ j*xlen2 + xlen2 - xlen + k ] = NO_SLIP;
		}
		// what remains for (if a||b==0) - case j==xlength+1:
		flagField [ k*xlen + xlen -1 ] = NO_SLIP;
	}
}
