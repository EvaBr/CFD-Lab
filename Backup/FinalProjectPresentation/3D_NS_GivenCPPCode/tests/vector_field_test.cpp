#include <iostream>
#include "data_structures.h"

#define size_x 10
#define size_y 10
#define size_z 10

using namespace std;

/** Vector comparison
 *
 * Reads up to three elements of the vector and returns 0 if all element
 * match between them.
 *
 * @param v1 First vector
 * @param v2 Other vector
 * @param dim Dimension of the vector / Number of elements to compare
 */
bool compare_vectors_fails ( FLOAT* v1, FLOAT* v2, int dim = 2 );

void setEntry ( FLOAT* entry, FLOAT val, int dim = 2 ){
    entry [1] = val;
    entry [0] = val - 1.0 / 3;
    if ( dim == 3 ) entry [2] = val + 1.0 / 3;
}

int main(){
    cout << endl << "Testing vector field" << endl;

    FLOAT entry[3];

    VectorField vfield2D ( size_x, size_y );
    VectorField vfield3D ( size_x, size_y, size_z );

    FLOAT d2counter = 1.0;
    FLOAT d3counter = 1.0;

    // Fill the arrays with data
    for (int i = 0; i < size_x; i++){
        for (int j = 0; j < size_y; j++){
            setEntry( entry, d2counter );
            vfield2D.setVector ( entry, i, j );
            d2counter ++;
            for (int k = 0; k < size_z; k++){
                setEntry( entry, d3counter , 3 );
                vfield3D.setVector ( entry, i, j, k );
                d3counter ++;
            }
        }
    }

    d2counter = 1.0;
    d3counter = 1.0;

    // Read the data out
    for (int i = 0; i < size_x; i++){
        for (int j = 0; j < size_y; j++){

            setEntry( entry, d2counter );
            if ( compare_vectors_fails ( entry, vfield2D.getVector (i,j) ) ){
                cerr << "Test for 2D array failed" << endl;
                exit (-1);
            }
            d2counter ++;

            for ( int k = 0; k < size_z; k++ ){

                setEntry( entry, d3counter , 3 );
                if ( compare_vectors_fails ( entry,
                            vfield3D.getVector (i,j,k), 3 ) ){
                    cerr << "Test for 2D array failed" << endl;
                    exit (-1);
                }
                d3counter ++;

            }
        }
    }

    cout << endl << "Vector field test finished successfully" << endl << endl;

    return 0;
}

bool compare_vectors_fails ( FLOAT* v1, FLOAT* v2, int dim ) {
    assertion ( ( dim == 2 ) || ( dim == 3 ) );
    for ( int i = 0; i < dim; i++ ){
        if ( v1 [i] != v2 [i] ) return 1;
    }
    return 0;
}
