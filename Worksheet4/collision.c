#include "collision.h"
#include "LBDefinitions.h"
#include "stdlib.h"

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
        for (int i=0; i<Q; i++)
                *(currentCell+i) = *(currentCell+i) - ( *(currentCell+i)-(*(feq+i)) ) / (*tau);
}

void doCollision(double *collideField, int *flagField,const double * const tau, int *subdomain){


	double density;
	double velocity[D];
	double feq[Q];
	double *currentCell = NULL; // currentCell points to the first distribution function within the respective cell
        for (int iz=1; iz<=subdomain[2]; iz++){
		for (int iy=1; iy<=subdomain[1]; iy++){
			for (int ix=1; ix<=subdomain[0]; ix++){
				// set pointer to current cell
				currentCell = collideField + Q * ( iz * (subdomain[0]+2) * (subdomain[1]+2) + iy * (subdomain[0]+2) + ix); //should we have
													// special func. to calculate the right index?

				// compute density, velocity and equilibrium prob. distrib. for this cell
				computeDensity ( currentCell, &density );
				computeVelocity ( currentCell, &density, velocity );
				computeFeq ( &density, velocity, feq );

				computePostCollisionDistributions ( currentCell, tau, feq );
			}
		}
	}
}

