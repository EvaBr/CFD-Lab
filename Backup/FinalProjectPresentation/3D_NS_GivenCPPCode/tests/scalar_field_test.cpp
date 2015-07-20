#include <iostream>
#include "data_structures.h"

#define size_x 10
#define size_y 10
#define size_z 10

using namespace std;

int main(){
    cout << endl << "Running scalar field test" << endl;

    ScalarField sfield2D ( size_x, size_y );
    ScalarField sfield3D ( size_x, size_y, size_z );

    FLOAT d2counter = 1.0;
    FLOAT d3counter = 1.0;

    cout << sfield2D.getNx() << endl;

    // Fill the fields completely with stuff
    for (int i = 0; i < size_x; i++){
        for (int j = 0; j < size_y; j++){

            sfield2D.getScalar(i,j) = d2counter;
            d2counter += 0.5;

            for ( int k = 0; k < size_z; k++ ){
                sfield3D.getScalar (i, j, k) = d3counter;
                d3counter += 0.5;
            }
        }
    }

    // Now read and see if the same values come out

    d2counter = 1.0;
    d3counter = 1.0;

    for (int i = 0; i < size_x; i++){
        for (int j = 0; j < size_y; j++){

            if ( sfield2D.getScalar ( i, j ) != d2counter ){
                std::cerr << "Error while reading the scalar values for 2D"
                    " field" << std::endl;
                return -1;
            }
            d2counter += 0.5;

            for ( int k = 0; k < size_z; k++ ){
                if ( sfield3D.getScalar ( i, j, k ) != d3counter ){
                    std::cerr << "Error while reading the scalar values for"
                       " 3D field" << std::endl;
                    return -1;
                }
                d3counter += 0.5;
            }
        }
    }

    std::cout << std::endl << "Test for scalar fields completed successfully" << std::endl << std::endl;

    return 0;
}
