#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

void treatBoundary(double *collideField, int* flagField,
                   const double * const wallVelocity, int *subdomain){
	int x, y, z, dx, dy, dz, i;
	//int len = xlength + 2;
	double cu = 0;
	double density;
	double *currentCell;
	int index;

	for (x = 0; x < subdomain[0]+2; x++){
        	for (y = 0; y < subdomain[1]+2; y++) {
			for (z = 0; z < subdomain[2]+2; z++) {
				index = z*(subdomain[0]+2)*(subdomain[1]+2) + y*(subdomain[0]+2) + x; //TODO
				currentCell = collideField + Q*index;

				for (i = 0; i < Q; ++i) {
					// save Lattice velocities for easier manipulation
					dx = LATTICEVELOCITIES[i][0];
					dy = LATTICEVELOCITIES[i][1];
					dz = LATTICEVELOCITIES[i][2];

					// check if index is valid (if particle is streamed to fluid)
					if ((x+dx > 0 && x+dx < subdomain[0]+1) && (y+dy > 0 && y+dy < subdomain[1]+1) && (z+dz > 0 && z+dz < subdomain[2]+1)){
						// deal with moving wall  boundary condition
						if (flagField[index] == MOVING_WALL){
							cu = wallVelocity[0]*LATTICEVELOCITIES[i][0] + wallVelocity[1]*LATTICEVELOCITIES[i][1] +
								wallVelocity[2]*LATTICEVELOCITIES[i][2];
							computeDensity(currentCell, &density);

							*(currentCell + i) = *(currentCell + Q*(dz*(subdomain[0]+2)*(subdomain[1]+2) + dy*(subdomain[0]+2) + dx) + Q-1-i)
										+ 2*LATTICEWEIGHTS[i]*density*cu/(C_S*C_S);
						}
						// deal with no-slip boundary condition
						else if (flagField[index] == NO_SLIP){
							*(currentCell + i) = *(currentCell + Q*(dz*(subdomain[0]+2)*(subdomain[1]+2) + dy*(subdomain[0]+2) + dx) + Q-1-i);
						}
					}
				}
			}
		}
	}
}
