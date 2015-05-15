#include "initLB.h"
#include "helper.h"
#include "LBDefinitions.h"

int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){
        // argument handling
        if (argc < 2 || argc > 2) {
                printf("When running the simulation, please give a valid file name to read from!\n");
                return 1;
        }

	char *FileName;
	FileName = argv[1];

	read_int ( FileName, "xlength", xlength );

	read_double ( FileName, "tau", tau );
	read_double ( FileName, "velocityWall", velocityWall );

	read_int ( FileName, "timesteps", timesteps );
	read_int ( FileName, "timestepsPerPlotting", timestepsPerPlotting );

	return 0;
}


void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength){
	//initialization of particle distribution func fields
	for (int a=0; a<=xlength+1; a++){
		for (int b=0; b<=xlength+1; b++){
			for (int c=0; c<=xlength+1; c++){
				for (int i=0; i<19; i++){
					streamField [ Q*(c*(xlength+2)*(xlength+2) + b*(xlength+2) + a) + i ] = LATTICEWEIGHTS [i];
					collideField [ Q*(c*(xlength+2)*(xlength+2) + b*(xlength+2) + a) + i ] = LATTICEWEIGHTS [i];
				}
				if (a==0 || b==0 || c==0 || a==xlength+1 || b==xlength+1) {
					//boundary no slip
					flagField [ c*(xlength+2)*(xlength+2) + b*(xlength+2) + a ] = NO_SLIP;
				} else {
					//set all inner points
					flagField [ c*(xlength+2)*(xlength+2) + b*(xlength+2) + a ] = FLUID;
				}
			}
			//if c==xlength, moving wall:
			flagField [ (xlength+2)*(xlength+2)*(xlength+1) + b*xlength + a ] = MOVING_WALL;
		}
	}
}
