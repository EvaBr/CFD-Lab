#include "initLB.h"
#include "helper.h"

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
  /* TODO */
}

