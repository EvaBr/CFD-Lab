#include "streaming.h"
#include "LBDefinitions.h"

void doStreaming(double *collideField, double *streamField, int *flagField, int *subdomain){
    int dx, dy, dz;
    double fi;
    // setting distribution function for each moving direction/lattice velocity of every particle
    // TODO: subdomain[2] for z instead??
    for (int z = 1; z <= subdomain[2]; ++z) {
        for (int y = 1; y <= subdomain[1]; ++y) {
            for (int x = 1; x <= subdomain[0]; ++x) {
                for (int i = 0; i < Q; ++i) {

                    // dx = c_i_x*dt, where dt = 1, etc.
                    dx = LATTICEVELOCITIES[i][0];
                    dy = LATTICEVELOCITIES[i][1];
                    dz = LATTICEVELOCITIES[i][2];

                    /*New value for our distribution function (DF) of the index 'i'

                     (We set it to DF(i) of the next particle, whose i-th lattice velocity
                     points towards considered particle (x,y,z))

                     Position of that next particle is given by (x-dx, y-dy, z-dz)*/

                    fi = collideField [ Q * ((z-dz)*(subdomain[0]+2)*(subdomain[1]+2) + (y-dy)*(subdomain[0]+2) + x - dx) + i ];
                    streamField [Q * (z*(subdomain[0]+2)*(subdomain[1]+2) + y*(subdomain[0]+2) + x) + i ] = fi;
                }
            }
        }
    }
}
