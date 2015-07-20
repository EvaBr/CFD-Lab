#include "BFInputStencils.h"

FLOAT computeVelocity3D (FlowField & flowField, int i, int j, int k, int stepSize,
                         const Parameters & parameters){
    if (j + parameters.parallel.firstCorner[1] >= stepSize + 2) {
        // Get the size of the inlet in Y. A 3 is subtracted because of the boundary cells
        const int inletYSize = parameters.geometry.sizeY - stepSize;
        const int inletZSize = parameters.geometry.sizeZ;

        const FLOAT y = (double)j + parameters.parallel.firstCorner[1] - 1.5 - (double)stepSize;
        const FLOAT z = (double)k + parameters.parallel.firstCorner[2] - 1.5;

        return 36.0 * parameters.walls.vectorLeft[0] /
                      (inletZSize * inletZSize * inletYSize * inletYSize) *
                      y * (y - inletYSize) * z * (z - inletZSize);
    } else {
        return 0.0;
    }
}

FLOAT computeVelocity2D (FlowField & flowField, int i, int j, int stepSize,
                         const Parameters & parameters){
    if (j + parameters.parallel.firstCorner[1] >= stepSize + 2) {
        // Get the size of the inlet in Y. A 3 is subtracted because of the boundary cells
        const int inletYSize = parameters.geometry.sizeY - stepSize;

        const FLOAT y = (double)j + parameters.parallel.firstCorner[1] - 1.5 - (double)stepSize;

        // DMITRIIS VERSION: for turbulence, please use: return parameters.walls.vectorLeft[0];
        return 6.0 * parameters.walls.vectorLeft[0] /
                     (inletYSize * inletYSize) * y * (inletYSize - y);
    } else {
        return 0.0;
    }
}

BFInputVelocityStencil::BFInputVelocityStencil (const Parameters & parameters) :
    BoundaryStencil<FlowField> (parameters),
    // Here, the obstacle size is set to zero if it was set as negative at the configuration
    _stepSize (parameters.bfStep.height > 0 ? parameters.bfStep.height : 0)
{}

// Most of the functions are empty, and they shouldn't be called, assuming that the input is always
// located at the left.

void BFInputVelocityStencil::applyLeftWall   ( FlowField & flowField, int i, int j ){
    flowField.getVelocity().getVector(i,j)[0] =
        computeVelocity2D(flowField, i, j, _stepSize, _parameters);
    flowField.getVelocity().getVector(i,j)[1] = -flowField.getVelocity().getVector(i+1,j)[1];
}

void BFInputVelocityStencil::applyRightWall  ( FlowField & flowField, int i, int j ){}
void BFInputVelocityStencil::applyBottomWall ( FlowField & flowField, int i, int j ){}
void BFInputVelocityStencil::applyTopWall    ( FlowField & flowField, int i, int j ){}

void BFInputVelocityStencil::applyLeftWall   ( FlowField & flowField, int i, int j, int k ){
    flowField.getVelocity().getVector(i,j,k)[0] =
        computeVelocity3D(flowField, i, j, k, _stepSize, _parameters);
    flowField.getVelocity().getVector(i,j,k)[1] = -flowField.getVelocity().getVector(i+1,j,k)[1];
    flowField.getVelocity().getVector(i,j,k)[2] = -flowField.getVelocity().getVector(i+1,j,k)[2];
}

void BFInputVelocityStencil::applyRightWall  ( FlowField & flowField, int i, int j, int k ){}
void BFInputVelocityStencil::applyBottomWall ( FlowField & flowField, int i, int j, int k ){}
void BFInputVelocityStencil::applyTopWall    ( FlowField & flowField, int i, int j, int k ){}
void BFInputVelocityStencil::applyFrontWall  ( FlowField & flowField, int i, int j, int k ){}
void BFInputVelocityStencil::applyBackWall   ( FlowField & flowField, int i, int j, int k ){}


BFInputFGHStencil::BFInputFGHStencil(const Parameters & parameters) :
    BoundaryStencil<FlowField> (parameters),
    _stepSize (parameters.bfStep.height > 0 ? parameters.bfStep.height : 0)
{}

void BFInputFGHStencil::applyLeftWall   ( FlowField & flowField, int i, int j ){
    flowField.getFGH().getVector(i,j)[0] =
        computeVelocity2D(flowField, i, j, _stepSize, _parameters);
}

void BFInputFGHStencil::applyRightWall  ( FlowField & flowField, int i, int j ){}
void BFInputFGHStencil::applyBottomWall ( FlowField & flowField, int i, int j ){}
void BFInputFGHStencil::applyTopWall    ( FlowField & flowField, int i, int j ){}

void BFInputFGHStencil::applyLeftWall  ( FlowField & flowField, int i, int j, int k ){
    flowField.getFGH().getVector(i,j,k)[0] =
        computeVelocity3D (flowField, i, j, k, _stepSize, _parameters);
}

void BFInputFGHStencil::applyRightWall  ( FlowField & flowField, int i, int j, int k ){}
void BFInputFGHStencil::applyBottomWall ( FlowField & flowField, int i, int j, int k ){}
void BFInputFGHStencil::applyTopWall    ( FlowField & flowField, int i, int j, int k ){}
void BFInputFGHStencil::applyFrontWall  ( FlowField & flowField, int i, int j, int k ){}
void BFInputFGHStencil::applyBackWall   ( FlowField & flowField, int i, int j, int k ){}
