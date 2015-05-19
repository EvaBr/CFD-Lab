#include "collision.h"
#include "LBDefinitions.h"

#include "stdlib.h" //added only because currentCell is initialised to NULL. to remove if unnecessary.
#include <stdio.h> //added for printf

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
//				currentCell = collideField + Q*(iz*xlength*xlength+iy*xlength+ix);
//				currentCell = &collideField[Q*(iz*(xlength+2)*(xlength+2)+iy*(xlength+2)+ix)]; // change from xlength to xlength+2
				currentCell = collideField + Q*(iz*(xlength+2)*(xlength+2)+iy*(xlength+2)+ix); // change from xlength to xlength+2a
				computeDensity (currentCell, &density);
/*if (density>1.1 || density<0.9)
	printf("warning at loop %d\n", iz);
else printf("okay after computeDensity\n");*/
				computeVelocity(currentCell, &density, velocity);
/*if (density>1.1 || density<0.9)
        printf("warning at loop %d\n", iz);
else printf("okay after computeVelocity\n");*/
				computeFeq(&density,velocity,feq);
/*if (density>1.1 || density<0.9)
        printf("warning at loop %d\n", iz);
else printf("okay after computeFeq\n");*/
				computePostCollisionDistributions(currentCell,tau,feq); //tau or &tau?
/*double temp=0;
//if (iz==1 && iy==1 && ix<=2)

	for (int i=0; i<Q; i++)
	{
		temp += *(currentCell+i);
	}
       	printf("currentCell density summation is %f\n", temp);
       	printf("density is %f\n", density);
*/
			}
		}
	}
  /* TODO */
}

