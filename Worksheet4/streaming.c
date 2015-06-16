#include "streaming.h"
#include "LBDefinitions.h"
#include "helper.h"

void doStreaming(double *collideField, double *streamField, int *flagField, int *subdomain){
    int dx, dy, dz;
    double fi;
    // setting distribution function for each moving direction/lattice velocity of every particle
    for (int z = 1; z <= subdomain[2]; z++) {
        for (int y = 1; y <= subdomain[1]; y++) {
            for (int x = 1; x <= subdomain[0]; x++) {
                for (int i = 0; i < Q; i++) {

                    // dx = c_i_x*dt, where dt = 1, etc.
                    dx = LATTICEVELOCITIES[i][0];
                    dy = LATTICEVELOCITIES[i][1];
                    dz = LATTICEVELOCITIES[i][2];

                    /*New value for our distribution function (DF) of the index 'i'

                     (We set it to DF(i) of the next particle, whose i-th lattice velocity
                     points towards considered particle (x,y,z))

                     Position of that next particle is (x-dx, y-dy, z-dz)*/

                    fi = collideField [ Q * compute_index(x-dx, y-dy, z-dz, subdomain) + i ];
                    streamField [Q * compute_index(x, y, z, subdomain) + i ] = fi;
                }
            }
        }
    }
}
