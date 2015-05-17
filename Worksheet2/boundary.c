#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"
void treatBoundary(double *collideField, int* flagField,
                   const double * const wallVelocity, int xlength){
    int x,y,z,dx,dy,dz,i;
    int N = xlength;
    double finv, cu = 0;
    double density;
    
    /* Two parallel walls */
    
    for (int k =0; k < 2; ++k) {
        x = k * (N + 1);
        for (y = 0; y < N + 2; ++y) {
            for (z = 0; z < N +1; ++z) {
                for (i = 0; i < Q; ++i) {
                    dx = LATTICEVELOCITIES[i][0];
                    dy = LATTICEVELOCITIES[i][1];
                    dz = LATTICEVELOCITIES[i][2];
                    
                    /* y-z plane*/
                    if (x+dx > 0 && x+dx< N+1 &&
                        y+dy > 0 && y+dy< N+1 &&
                        z+dz > 0 && z+dz< N+1) {
                       
                        finv = collideField[Q *((z+dz)*(N+2)*(N+2) + (y+dy)*(N+2)+x+dx)+ Q-i-1];
                        collideField[Q *( z*(N+2)*(N+2) + y*(N+2)+ x)+ i] =finv;
                    }
                    /*x-z plane swap x and y*/
                    if (y+dx > 0 && y+dx < N + 1 &&
                        x+dy > 0 && x+dy < N + 1 &&
                        z+dz > 0 && z+dz < N + 1) {
                            
                        finv = collideField[Q * ((z+dz)*(N+2)*(N+2) + (x+dy)*(N+2) + y+dx) + Q-i-1];
                        collideField[Q * (z*(N+2)*(N+2) + x*(N+2) + y) + i] = finv;
                    }
                    /*x-y plane swap x and z*/
                    if (z+dx > 0 && z+dx < N + 1 &&
                        y+dy > 0 && y+dy < N + 1 &&
                        x+dz > 0 && x+dz < N + 1) {
                        
                        finv = collideField[Q * ((x+dz)*(N+2)*(N+2) + (y+dy)*(N+2) + z+dx) + Q-i-1];
                        
                        
                        if(k == 1) {
                            currentCell = collideField + Q * ((x+dz)*(N+2)*(N+2) + (y+dy)*(N+2) + z+dx);
                            computeDensity(currentCell, &density);
                            cu = 0; /*dot product of c_i and velocity of wall*/
                            for (int d = 0; d < D; ++d) {
                                cu += LATTICEVELOCITIES[i][d] * wallVelocity[d];
                            }
                            finv += 2*LATTICEWEIGHTS[i]*density*cu/(C_S*C_S);
                        }
                        collideField[Q * (x*(N+2)*(N+2) + y*(N+2) + z) + i] = finv;
                   }
            }
        }
    }
}