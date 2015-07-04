#include <time.h>
#include <iostream>
#include "../flow_field.h"
#include "../iterators.h"
#include "../stencils/vtk_stencil.h"
#include "../parameters.h"


int main () {
    std::cout << "Starting VTK test" << std::endl;
    FlowField flowField ( 10, 10, 10 );

    clock_t start = clock();

    FLOAT velocity [3] = {1,1,1};

    for (int k = 0; k < flowField.getNz() + 3; k++ ){
        for (int j = 0; j < flowField.getNy() + 3; j++ ){
            for (int i = 0; i < flowField.getNx() + 3; i++ ){
                flowField.getPressure().getScalar(i,j,k) = (double) k;
                flowField.getVelocity().setVector(velocity, i,j,k);
            }
        }
    }

    std::cout << "Initialization time: " << (double) (clock() - start) / CLOCKS_PER_SEC
        << std::endl;
    start = clock();

    Parameters parameters;

    parameters.dx = 1;
    parameters.dy = 1;
    parameters.dz = 1;

    VTKStencil stencil( "/tmp/some_file", parameters );

    std::cout << "Stencil creation time: " << (double) (clock() - start) / CLOCKS_PER_SEC << std::endl;
    start = clock();

    stencil.openFile ( flowField, 5.0/3 );

    std::cout << "File-openning and grid data writing time: " << (double) (clock() - start) / CLOCKS_PER_SEC << std::endl;
    start = clock();

    FieldIterator iterator( flowField, stencil );
    iterator.iterateInnerCells();

    std::cout << "Iteration time: " << (double) (clock() - start) / CLOCKS_PER_SEC << std::endl;
    start = clock();

    stencil.write( flowField );
    std::cout << "Writing time: " << (double) (clock() - start) / CLOCKS_PER_SEC << std::endl;

    stencil.closeFile();

}
