#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

void treatBoundary(double *collideField, int* flagField,
                   const double * const wallVelocity, int xlength){
	int x, y, z, dx, dy, dz, i;
	int len = xlength + 2;
	double cu = 0;
	double density;
	double *currentCell;
	int index;

	for (x = 0; x < len; x++){
        	for (y = 0; y < len; y++) {
			for (z = 0; z < len; z++) {
				index = z*len*len + y*len + x;
				currentCell = collideField + Q*index;

				for (i = 0; i < Q; ++i) {
					// save Lattice velocities for easier manipulation
					dx = LATTICEVELOCITIES[i][0];
					dy = LATTICEVELOCITIES[i][1];
					dz = LATTICEVELOCITIES[i][2];

					// check if index is valid (if particle is streamed to fluid)
					if ((x+dx > 0 && x+dx<len-1) && (y+dy > 0 && y+dy< len-1) && (z+dz > 0 && z+dz< len-1)){
						// deal with moving wall  boundary condition
						if (flagField[index] == MOVING_WALL){
							cu = wallVelocity[0]*LATTICEVELOCITIES[i][0] + wallVelocity[1]*LATTICEVELOCITIES[i][1] +
								wallVelocity[2]*LATTICEVELOCITIES[i][2];
							computeDensity(currentCell, &density);

							*(currentCell + i) = *(currentCell + Q*(dz*len*len + dy*len + dx) + Q-1-i) +
											2*LATTICEWEIGHTS[i]*density*cu/(C_S*C_S);
						}
						// deal with no-slip boundary condition
						else if (flagField[index] == NO_SLIP){
							*(currentCell + i) = *(currentCell + Q*(dz*len*len + dy*len + dx) + Q-1-i );
						}
					}
				}
			}
		}
	}
}
