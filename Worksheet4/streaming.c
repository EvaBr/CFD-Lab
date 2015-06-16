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


void  swap( double ** sendBuffer, double ** readBuffer, int *subdomain, int *flagField , int *rank, int *proc ){
    
    int x= subdomain[0];
    int y= subdomain[1];
    int z= subdomain[2];
    
    // Left
    if (flagField[compute_index(0, y/2, z/2, subdomain)] == PARALLEL_BOUNDARY) {
        /*MPI_Send(&sendBuffer[0][0], 5 * (subdomain[1] + 2) * (subdomain[2]+ 2), MPI_DOUBLE, rank - 1, 0,
                 MPI_COMM_WORLD);
        
        MPI_Recv(&readBuffer[0][0], 5 * (subdomain[1] + 2) * (subdomain[2] + 2), MPI_DOUBLE, rank - 1, 0,
                 MPI_COMM_WORLD, &status);*/
        
        MPI_Sendrecv(&sendBuffer[0][0], 5 * (subdomain[1] + 2) * (subdomain[2]+ 2), MPI_DOUBLE, rank - 1, 0,
                     &readBuffer[0][0], 5 * (subdomain[1] + 2) * (subdomain[2] + 2), MPI_DOUBLE, rank - 1, 0,
                     MPI_COMM_WORLD, &status);
    }

    
    // Right
    if (flagField[compute_index(x, y/2, z/2, subdomain)] == PARALLEL_BOUNDARY) {
        /*MPI_Send(&sendBuffer[1][0], 5 * (subdomain[1] + 2) * (subdomain[2]+ 2), MPI_DOUBLE, rank + 1, 0,
                 MPI_COMM_WORLD);
        
        MPI_Recv(&readBuffer[1][0], 5 * (subdomain[1] + 2) * (subdomain[2] + 2), MPI_DOUBLE, rank + 1, 0,
                 MPI_COMM_WORLD, &status);*/
        
        
        MPI_Sendrecv(&sendBuffer[1][0], 5 * (subdomain[1] + 2) * (subdomain[2]+ 2), MPI_DOUBLE, rank + 1, 0,
                     &readBuffer[1][0], 5 * (subdomain[1] + 2) * (subdomain[2] + 2), MPI_DOUBLE, rank + 1, 0,
                     MPI_COMM_WORLD, &status);
    }

    
    // Top
    if (flagField[compute_index(x/2, y/2, z, subdomain] == PARALLEL_BOUNDARY) {
        /*MPI_Send(&sendBuffer[2][0], 5 * (subdomain[0] + 2) * (subdomain[1]+ 2), MPI_DOUBLE, rank + proc[0]*proc[1], 0, MPI_COMM_WORLD);
        
        MPI_Recv(&readBuffer[2][0], 5 * (subdomain[0] + 2) * (subdomain[1] + 2), MPI_DOUBLE, rank + proc[0]*proc[1], 0, MPI_COMM_WORLD, &status);*/
                  
        MPI_Sendrecv(&sendBuffer[2][0], 5 * (subdomain[0] + 2) * (subdomain[1]+ 2), MPI_DOUBLE, rank + proc[0]*proc[1], 0,&readBuffer[2][0], 5 * (subdomain[0] + 2) * (subdomain[1] + 2), MPI_DOUBLE, rank + proc[0]*proc[1], 0, MPI_COMM_WORLD, &status);
                  
    }

    
    // Bottom
    if(flagField[compute_index(x/2, y/2, 0, subdomain] == PARALLEL_BOUNDARY){
       /*MPI_Send(&sendBuffer[3][0], 5 * (subdomain[0] + 2) * (subdomain[1]+ 2), MPI_DOUBLE, rank - proc[0]*proc[1], 0, MPI_COMM_WORLD);
                  
       MPI_Recv(&readBuffer[3][0], 5 * (subdomain[0] + 2) * (subdomain[1] + 2), MPI_DOUBLE, rank - proc[0]*proc[1], 0, MPI_COMM_WORLD, &status);*/
                 
                 
                 
       MPI_Sendrecv(&sendBuffer[3][0], 5 * (subdomain[0] + 2) * (subdomain[1]+ 2), MPI_DOUBLE, rank - proc[0]*proc[1], 0, &readBuffer[3][0], 5 * (subdomain[0] + 2) * (subdomain[1] + 2), MPI_DOUBLE, rank - proc[0]*proc[1], 0, MPI_COMM_WORLD, &status);
                  
    }
            


    // Front
    if(flagField[compute_index(x/2, 0, z/2, subdomain] == PARALLEL_BOUNDARY){
                 
       /*MPI_Send(&sendBuffer[4][0], 5 * (subdomain[0] + 2) * (subdomain[2]+ 2), MPI_DOUBLE, rank - proc[1], 0,MPI_COMM_WORLD);
                 
       MPI_Recv(&readBuffer[4][0], 5 * (subdomain[0] + 2) * (subdomain[2] + 2), MPI_DOUBLE, rank - proc[1], 0,MPI_COMM_WORLD, &status);*/
                 
                 
       MPI_Sendrecv(&sendBuffer[4][0], 5 * (subdomain[0] + 2) * (subdomain[2]+ 2), MPI_DOUBLE, rank - proc[1], 0,&readBuffer[4][0], 5 * (subdomain[0] + 2) * (subdomain[2] + 2), MPI_DOUBLE, rank - proc[1], 0,MPI_COMM_WORLD, &status);
                              
                              
    }
                 
    
    //Back
                 
    if(flagField[compute_index(x/2, y, z/2, subdomain] == PARALLEL_BOUNDARY){
        
      /* MPI_Send(&sendBuffer[5][0], 5 * (subdomain[0] + 2) * (subdomain[2]+ 2), MPI_DOUBLE, rank + proc[1], 0, MPI_COMM_WORLD);
                              
       MPI_Recv(&readBuffer[5][0], 5 * (subdomain[0] + 2) * (subdomain[2] + 2), MPI_DOUBLE, rank + proc[1], 0,MPI_COMM_WORLD, &status);*/
                 
                 
       MPI_Sendrecv(&sendBuffer[5][0], 5 * (subdomain[0] + 2) * (subdomain[2]+ 2), MPI_DOUBLE, rank + proc[1], 0, &readBuffer[5][0], 5 * (subdomain[0] + 2) * (subdomain[2] + 2), MPI_DOUBLE, rank + proc[1], 0,MPI_COMM_WORLD, &status);
                 
    }
                 

    
    


}