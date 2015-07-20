#include "FlowField.h"
#include "DataStructures.h"
#include <iostream>

#define size_x 20
#define size_y 25

// Some basic tests for the flow field. Namely, that the members can be read
// and written at some points.

int main ( int argc, char * argv[] ){

    FlowField field ( size_x, size_y );
    std::cout << field.getPressure().getScalar ( 10,10 ) << std::endl;
    field.getPressure().getScalar ( 10,10 ) = 10;
    std::cout << field.getPressure().getScalar ( 10,10 ) / 2 << std::endl;

    std::cout << field.getFlags().getNx () << std::endl;
    std::cout << field.getFlags().getNy () << std::endl;
    field.getFlags().getValue ( 10,10 ) = 7;
    std::cout << field.getFlags().getValue ( 10,10 ) / 3 << std::endl;

}
