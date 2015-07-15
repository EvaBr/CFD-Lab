#include "VelocityStencil.h"

VelocityStencil::VelocityStencil ( const Parameters & parameters ) : FieldStencil<FlowField> ( parameters ) {}


void VelocityStencil::apply ( FlowField & flowField, int i, int j ){

    const FLOAT dt = _parameters.timestep.dt;
    const FLOAT dx = _parameters.geometry.dx, dy = _parameters.geometry.dy;
    const int obstacle = flowField.getFlags().getValue(i, j);
    VectorField & velocity = flowField.getVelocity();

    if ((obstacle & OBSTACLE_SELF) == 0){ // If this is a fluid cell
        if ((obstacle & OBSTACLE_RIGHT) == 0){  // Check whether the neighbor is also fluid
            flowField.getVelocity().getVector(i,j)[0] = flowField.getFGH().getVector(i,j)[0] - dt/dx *
                (flowField.getPressure().getScalar(i+1,j) - flowField.getPressure().getScalar(i,j));
        } else {    // Otherwise, set to zero
            flowField.getVelocity().getVector(i,j)[0] = 0;
        }   // Note that we only set one direction per cell. The neighbor at the left is
            // responsible for the other side
        if ((obstacle & OBSTACLE_TOP) == 0){
            flowField.getVelocity().getVector(i,j)[1] = flowField.getFGH().getVector(i,j)[1] - dt/dy *
                (flowField.getPressure().getScalar(i,j+1) - flowField.getPressure().getScalar(i,j));
        } else {
            flowField.getVelocity().getVector(i,j)[1] = 0;
        }
    } else {    // Otherwise we have an obstacle cell
        if ((obstacle & OBSTACLE_LEFT) == 0){   // If the left neighbor is fluid
            velocity.getVector(i,j)[1] = -velocity.getVector(i-1,j)[1];
            // Here, we don't set the normal velocity to zero, since it has already been taken care
            // of in the previous part of the branch
        }
        if ((obstacle & OBSTACLE_RIGHT) == 0){
            velocity.getVector(i,j)[0] = 0.0;
            velocity.getVector(i,j)[1] = -velocity.getVector(i+1,j)[1];
        }
        if ((obstacle & OBSTACLE_BOTTOM) == 0){
            velocity.getVector(i,j)[0] = -velocity.getVector(i,j-1)[0];
        }
        if ((obstacle & OBSTACLE_TOP) == 0){
            velocity.getVector(i,j)[0] = -velocity.getVector(i,j+1)[0];
            velocity.getVector(i,j)[1] = 0.0;
        }
    }
}


void VelocityStencil::apply ( FlowField & flowField, int i, int j, int k ){

    const FLOAT dt = _parameters.timestep.dt;
    const FLOAT dx = _parameters.geometry.dx, dy = _parameters.geometry.dy,
                dz = _parameters.geometry.dz;
    const int obstacle = flowField.getFlags().getValue(i, j, k);
    VectorField & velocity = flowField.getVelocity();

    if ((obstacle & OBSTACLE_SELF) == 0) {
        if ((obstacle & OBSTACLE_RIGHT) == 0) {
            flowField.getVelocity().getVector(i,j,k)[0] = flowField.getFGH().getVector(i,j,k)[0] - dt/dx *
                (flowField.getPressure().getScalar(i+1,j,k) - flowField.getPressure().getScalar(i,j,k));
        } else {
            flowField.getVelocity().getVector(i, j, k)[0] = 0.0;
        }
        if ((obstacle & OBSTACLE_TOP) == 0) {
            flowField.getVelocity().getVector(i,j,k)[1] = flowField.getFGH().getVector(i,j,k)[1] - dt/dy *
                (flowField.getPressure().getScalar(i,j+1,k) - flowField.getPressure().getScalar(i,j,k));
        } else {
            flowField.getVelocity().getVector(i, j, k)[1] = 0.0;
        }
        if ((obstacle & OBSTACLE_BACK) == 0) {
            flowField.getVelocity().getVector(i,j,k)[2] = flowField.getFGH().getVector(i,j,k)[2] - dt/dz *
                (flowField.getPressure().getScalar(i,j,k+1) - flowField.getPressure().getScalar(i,j,k));
        } else {
            flowField.getVelocity().getVector(i, j, k)[2] = 0.0;
        }
    } else {    // It this is an obstacle cell
        if ((obstacle & OBSTACLE_LEFT) == 0) { // If the left neighbour is a fluid cell
            velocity.getVector(i, j, k)[1] = -velocity.getVector(i-1, j, k)[1];
            velocity.getVector(i, j, k)[2] = -velocity.getVector(i-1, j, k)[2];
        }
        if ((obstacle & OBSTACLE_RIGHT) == 0) {
            velocity.getVector(i, j, k)[0] = 0.0;
            velocity.getVector(i, j, k)[1] = -velocity.getVector(i+1, j, k)[1];
            velocity.getVector(i, j, k)[2] = -velocity.getVector(i+1, j, k)[2];
        }
        if ((obstacle & OBSTACLE_BOTTOM) == 0) {
            velocity.getVector(i, j, k)[0] = -velocity.getVector(i, j-1, k)[0];
            velocity.getVector(i, j, k)[2] = -velocity.getVector(i, j-1, k)[2];
        }
        if ((obstacle & OBSTACLE_TOP) == 0) {
            velocity.getVector(i, j, k)[0] = -velocity.getVector(i, j+1, k)[0];
            velocity.getVector(i, j, k)[1] = 0.0;
            velocity.getVector(i, j, k)[2] = -velocity.getVector(i, j+1, k)[2];
        }
        if ((obstacle & OBSTACLE_FRONT) == 0) {
            velocity.getVector(i, j, k)[0] = -velocity.getVector(i, j, k-1)[0];
            velocity.getVector(i, j, k)[1] = -velocity.getVector(i, j, k-1)[1];
        }
        if ((obstacle & OBSTACLE_BACK) == 0) {
            velocity.getVector(i, j, k)[0] = -velocity.getVector(i, j, k+1)[0];
            velocity.getVector(i, j, k)[1] = -velocity.getVector(i, j, k+1)[1];
            velocity.getVector(i, j, k)[2] = 0.0;
        }
    }
}
