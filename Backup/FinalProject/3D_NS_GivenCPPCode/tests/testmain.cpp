#include <iostream>
#include <vector>
#include "../data_structures.h"
#include "../definitions.h"
#include "../parameters.h"
#include "../flow_field.h"
#include "../stencils/vtk_stencil.h"
#include "../stencils/moving_wall_stencil.h"
#include "../iterators.h"
#include "../stencils/derivatives.h"

#define size_x 20
#define size_y 25
#define size_z 15

bool compare_vectors_fails ( FLOAT* v1, FLOAT* v2, int dim = 2 ) {
    assertion ( ( dim == 2 ) || ( dim == 3 ) );
    for ( int i = 0; i < dim; i++ ){
        if ( v1 [i] != v2 [i] ) return 1;
    }
    return 0;
}

void setEntry ( FLOAT* entry, FLOAT val, int dim = 2 ){
    entry [1] = val;
    entry [0] = val - 1.0 / 3;
    if ( dim == 3 ) entry [2] = val + 1.0 / 3;
}

class Test {
  public:
    virtual void validate() = 0;
    virtual std::string getTestName() const = 0;
};

class FlowFieldTest: public Test {
  public:
    virtual void validate(){
        FlowField field ( size_x, size_y );
        std::cout << field.getPressure().getScalar ( 10,10 ) << std::endl;
        field.getPressure().getScalar ( 10,10 ) = 10;
        std::cout << field.getPressure().getScalar ( 10,10 ) / 2 << std::endl;

        std::cout << field.getFlags().getNx () << std::endl;
        std::cout << field.getFlags().getNy () << std::endl;
        field.getFlags().getValue ( 10,10 ) = 7;
        std::cout << field.getFlags().getValue ( 10,10 ) / 3 << std::endl;
    }
    virtual std::string getTestName() const {return "flow-field-test";}
};

class ScalarFieldTest : public Test {
    public:
        virtual void validate(){
            std::cout << std::endl << "Running scalar field test" << std::endl;

            ScalarField sfield2D ( size_x, size_y );
            ScalarField sfield3D ( size_x, size_y, size_z );

            FLOAT d2counter = 1.0;
            FLOAT d3counter = 1.0;

            std::cout << sfield2D.getNx() << std::endl;

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
                        // return -1;
                    }
                    d2counter += 0.5;

                    for ( int k = 0; k < size_z; k++ ){
                        if ( sfield3D.getScalar ( i, j, k ) != d3counter ){
                            std::cerr << "Error while reading the scalar values for"
                                " 3D field" << std::endl;
                            // return -1;
                        }
                        d3counter += 0.5;
                    }
                }
            }

            std::cout << std::endl << "Test for scalar fields completed successfully" << std::endl << std::endl;
        }

        virtual std::string getTestName() const {return "scalar field test";}
};


class VectorFieldTest : public Test{
    public:

        virtual void validate () {
            std::cout << std::endl << "Testing vector field" << std::endl;

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
                        std::cerr << "Test for 2D array failed" << std::endl;
                        exit (-1);
                    }
                    d2counter ++;

                    for ( int k = 0; k < size_z; k++ ){

                        setEntry( entry, d3counter , 3 );
                        if ( compare_vectors_fails ( entry,
                                    vfield3D.getVector (i,j,k), 3 ) ){
                            std::cerr << "Test for 2D array failed" << std::endl;
                            exit (-1);
                        }
                        d3counter ++;

                    }
                }
            }

            std::cout << std::endl << "Vector field test finished successfully" << std::endl << std::endl;
        }

        virtual std::string getTestName() const {return "Vector field test";}
};


class VTKTest : public Test {
    public:

        virtual void validate () {
            std::cout << "Starting VTK test" << std::endl;
            FlowField flowField ( 50, 50, 50 );

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
            iterator.iterate();

            std::cout << "Iteration time: " << (double) (clock() - start) / CLOCKS_PER_SEC << std::endl;
            start = clock();

            stencil.write( flowField );
            std::cout << "Writing time: " << (double) (clock() - start) / CLOCKS_PER_SEC << std::endl;

            stencil.closeFile();
        }

        virtual std::string getTestName() const {return "VTK test";}
};


class BoundaryTest : public Test {

    public:

        virtual void validate () {

            std::cout << "Starting boundary test" << std::endl;

            FlowField flowField2D (5,5);
            FlowField flowField3D (5,5,5);

            Parameters parameters;
            parameters.dx = 1;
            parameters.dy = 1;
            parameters.dz = 1;
            parameters.vector_top[0] = 1;
            parameters.vector_top[1] = 2;

            parameters.vector_bottom[0] = -1;
            parameters.vector_bottom[1] = 0;

            parameters.vector_right[0] = -1;
            parameters.vector_right[1] = 1.5;

            parameters.vector_front[0] = 4;
            parameters.vector_front[1] = 2.5;

            MovingWallStencil stencil (parameters);
            GlobalBoundaryIterator iterator2D (flowField2D, stencil);
            GlobalBoundaryIterator iterator3D (flowField3D, stencil);
            iterator2D.iterate();
            iterator3D.iterate();

            for (int j = flowField2D.getNy()+2; j >= 1; j--){
                for (int i = 1; i <= flowField2D.getNx()+2; i++){
                    std::cout << flowField2D.getVelocity().getVector(i,j)[0] << "\t";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
            for (int j = flowField2D.getNy()+2; j >= 1; j--){
                for (int i = 1; i <= flowField2D.getNx()+2; i++){
                    std::cout << flowField2D.getVelocity().getVector(i,j)[1] << "\t";
                }
                std::cout << std::endl;
            }

            std::cout << std::endl;
            for (int j = flowField3D.getNy()+2; j >= 1; j--){
                for (int i = 1; i <= flowField3D.getNx()+2; i++){
                    std::cout << flowField3D.getVelocity().getVector(i,j,1)[0] << "\t";
                }
                std::cout << std::endl;
            }
        }
        virtual std::string getTestName() const {return "Boundary test";}
};

class DerivativesTest : public Test {

    public:

        DerivativesTest(){
            eps = 1e-10;
            parameters.dx = 1;
            parameters.dy = 1;
            parameters.dz = 1;
            parameters.dt = 1;
        }

        void setAllValues(double val){
            for (int i = 0; i < 81; i++){
                lv[i] = val;
            }
        }

        void printResults(){
            std::cout <<
             "dudx  : " << dudx  (lv, parameters) << std::endl <<
             "dvdy  : " << dvdy  (lv, parameters) << std::endl <<
             "dwdz  : " << dwdz  (lv, parameters) << std::endl <<

             "d2udx2: " << d2udx2(lv, parameters) << std::endl <<
             "d2udy2: " << d2udy2(lv, parameters) << std::endl <<
             "d2udz2: " << d2udz2(lv, parameters) << std::endl <<
             "d2vdx2: " << d2vdx2(lv, parameters) << std::endl <<
             "d2vdy2: " << d2vdy2(lv, parameters) << std::endl <<
             "d2vdz2: " << d2vdz2(lv, parameters) << std::endl <<
             "d2wdx2: " << d2wdx2(lv, parameters) << std::endl <<
             "d2wdy2: " << d2wdy2(lv, parameters) << std::endl <<
             "d2wdz2: " << d2wdz2(lv, parameters) << std::endl <<

             "duvdx : " << duvdx (lv, parameters) << std::endl <<
             "duvdy : " << duvdy (lv, parameters) << std::endl <<
             "duwdx : " << duwdx (lv, parameters) << std::endl <<
             "duwdz : " << duwdz (lv, parameters) << std::endl <<
             "dvwdy : " << dvwdy (lv, parameters) << std::endl <<
             "dvwdz : " << dvwdz (lv, parameters) << std::endl <<

             "du2dx : " << du2dx (lv, parameters) << std::endl <<
             "dv2dy : " << dv2dy (lv, parameters) << std::endl <<
             "dw2dz : " << dw2dz (lv, parameters) << std::endl <<

            std::endl;
        }

        virtual std::string getTestName() const {return "derivatives-test";}

        virtual void validate(){
            setAllValues(0);
            for (int i = -1; i < 2; i++){
                for (int j = -1; j < 2; j++){
                    for (int k = -1; k < 2; k++){
                        for (int c = 0; c < 3; c++){
                            lv[mapd(i,j,k,c)] = 1;
                        }
                    }
                }
            }

            for (int i = 0; i < 81; i++){
                if (!lv[i]) {
                    std::cerr << "Failed accessing all elements from coordinates and components"
                        << std::endl;
                }
            }

            lv[mapd(1,0,0,0)]  = 3;
            lv[mapd(0,0,0,0)]  = 2;
            lv[mapd(-1,0,0,0)] = 1;

            // TODO and here it does something odd.

            printResults();
        }

    private:

        FLOAT eps;
        Parameters parameters;
        FLOAT lv[81];     // Local velocity
};


int main(int argc, char * argv[]){
    std::vector<Test*> tests;

    tests.push_back(new ScalarFieldTest());
    tests.push_back(new VectorFieldTest());
    // tests.push_back(new VTKTest());
    tests.push_back(new FlowFieldTest());
    tests.push_back(new BoundaryTest());
    tests.push_back(new DerivativesTest());

    for (unsigned int i = 0; i<tests.size(); i++){
        std::cout << "Running " << tests[i] -> getTestName() << std::endl;
        tests[i]->validate();
        std::cout << std::endl;
    }
}
