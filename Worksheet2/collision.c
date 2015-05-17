#include "collision.h"
#include "LBDefinitions.h"

#include "stdlib.h" //added only because currentCell is initialised to NULL. to remove if unnecessary.


void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
  /* TODO */

        for (int i=0; i<Q; i++)
                *(currentCell+i) = *(currentCell+i) - ( *(currentCell+i)-(*(feq+i)) ) / (*tau);
}

void doCollision(double *collideField, int *flagField,const double * const tau,int xlength){


	double density;
  /* TODO: replace hardcoded dimension of velocity, e.g., changing from 3 to D after defining D=3?  */
	double velocity[D];
	double feq[Q];
	double *currentCell=NULL; //currentCell points to the first distribution function within the respective cell
	
        for (int iz=1; iz<=xlength; iz++){
		for (int iy=1; iy<=xlength; iy++){
			for (int ix=1; ix<=xlength; ix++){
				currentCell = collideField + Q*(iz*xlength*xlength+iy*xlength+ix);
				computeDensity (currentCell, &density);
				computeVelocity(currentCell, &density, velocity);
				computeFeq(&density,velocity,feq);
				computePostCollisionDistributions(currentCell,tau,feq); //tau or &tau?
			}
		}
	}
  /* TODO */
}

