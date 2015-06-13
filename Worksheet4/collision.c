#include "collision.h"
#include "LBDefinitions.h"
#include "stdlib.h"

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
        for (int i=0; i<Q; i++)
                *(currentCell+i) = *(currentCell+i) - ( *(currentCell+i)-(*(feq+i)) ) / (*tau);
}

void doCollision(double *collideField, int *flagField,const double * const tau,int xlength){


	double density;
	double velocity[D];
	double feq[Q];
	double *currentCell = NULL; // currentCell points to the first distribution function within the respective cell
        for (int iz=1; iz<=xlength; iz++){
		for (int iy=1; iy<=xlength; iy++){
			for (int ix=1; ix<=xlength; ix++){
				// set pointer to current cell
				currentCell = collideField + Q*(iz*(xlength+2)*(xlength+2) + iy*(xlength+2) + ix);

				// compute density, velocity and equilibrium prob. distrib. for this cell
				computeDensity (currentCell, &density);
				computeVelocity (currentCell, &density, velocity);
				computeFeq (&density, velocity, feq);

				computePostCollisionDistributions (currentCell, tau, feq);
			}
		}
	}
}

