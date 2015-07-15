#ifndef _DERIVATIVES_H_
#define _DERIVATIVES_H_

#include <math.h>
#include "../Definitions.h"
#include "../Parameters.h"

// Load the local velocity cube with relevant velocities of the 2D plane
inline void loadLocalVelocity2D(FlowField & flowField, FLOAT * const localVelocity, int i, int j){
    for (int row = -1; row <= 1; row++ ){
        for ( int column = -1; column <= 1; column ++ ){
            const FLOAT * const point = flowField.getVelocity().getVector(i + column, j + row);
            localVelocity[39 + 9*row + 3*column]     = point[0]; // x-component
            localVelocity[39 + 9*row + 3*column + 1] = point[1]; // y-component
        }
    }
}

// Load the local velocity cube with surrounding velocities
inline void loadLocalVelocity3D(FlowField & flowField, FLOAT * const localVelocity, int i, int j, int k){
    for ( int layer = -1; layer <= 1; layer ++ ){
        for ( int row = -1; row <= 1; row++ ){
            for ( int column = -1; column <= 1; column ++ ){
                const FLOAT * const point = flowField.getVelocity().getVector(i + column, j + row, k + layer);
                localVelocity[39 + 27*layer + 9*row + 3*column    ] = point[0]; // x-component
                localVelocity[39 + 27*layer + 9*row + 3*column + 1] = point[1]; // y-component
                localVelocity[39 + 27*layer + 9*row + 3*column + 2] = point[2]; // z-component
            }
        }
    }
}


// Maps an index and a component to the corresponding value in the cube.
inline int mapd (int i, int j, int k, int component){
    return 39 + 27*k + 9*j + 3*i + component;
}

// Derivative functions. They are applied to a cube of 3x3x3 cells.
inline FLOAT dudx ( const FLOAT * const lv, const Parameters & parameters ) {
    return ( lv [mapd(0,0,0,0)] - lv [mapd(-1,0,0,0)] ) / parameters.geometry.dx;
}

inline FLOAT dvdy ( const FLOAT * const lv, const Parameters & parameters ) {
    return ( lv [mapd(0,0,0,1)] - lv [mapd(0,-1,0,1)] ) / parameters.geometry.dy;
}

inline FLOAT dwdz ( const FLOAT * const lv, const Parameters & parameters ) {
    return ( lv [mapd(0,0,0,2)] - lv [mapd(0,0,-1,2)] ) / parameters.geometry.dz;
}



inline FLOAT d2udx2 ( const FLOAT * const lv, const Parameters & parameters ) {
    return ( lv[mapd(1,0,0,0)] - 2*lv[mapd(0,0,0,0)] + lv[mapd(-1,0,0,0)] )
        / ( parameters.geometry.dx * parameters.geometry.dx );
}

inline FLOAT d2udy2 ( const FLOAT * const lv, const Parameters & parameters ) {
    return ( lv[mapd(0,1,0,0)] - 2*lv[mapd(0,0,0,0)] + lv[mapd(0,-1,0,0)] )
        / ( parameters.geometry.dy * parameters.geometry.dy );
}

inline FLOAT d2udz2 ( const FLOAT * const lv, const Parameters & parameters ) {
    return ( lv[mapd(0,0,1,0)] - 2*lv[mapd(0,0,0,0)] + lv[mapd(0,0,-1,0)] )
        / ( parameters.geometry.dz * parameters.geometry.dz );
}

inline FLOAT d2vdx2 ( const FLOAT * const lv, const Parameters & parameters ) {
    return ( lv[mapd(1,0,0,1)] - 2*lv[mapd(0,0,0,1)] + lv[mapd(-1,0,0,1)] )
        / ( parameters.geometry.dx * parameters.geometry.dx );
}

inline FLOAT d2vdy2 ( const FLOAT * const lv, const Parameters & parameters ) {
    return ( lv[mapd(0,1,0,1)] - 2*lv[mapd(0,0,0,1)] + lv[mapd(0,-1,0,1)] )
        / ( parameters.geometry.dy * parameters.geometry.dy );
}

inline FLOAT d2vdz2 ( const FLOAT * const lv, const Parameters & parameters ) {
    return ( lv[mapd(0,0,1,1)] - 2*lv[mapd(0,0,0,1)] + lv[mapd(0,0,-1,1)] )
        / ( parameters.geometry.dz * parameters.geometry.dz );
}

inline FLOAT d2wdx2 ( const FLOAT * const lv, const Parameters & parameters ) {
    return ( lv[mapd(1,0,0,2)] - 2*lv[mapd(0,0,0,2)] + lv[mapd(-1,0,0,2)] )
        / ( parameters.geometry.dx * parameters.geometry.dx );
}

inline FLOAT d2wdy2 ( const FLOAT * const lv, const Parameters & parameters ) {
    return ( lv[mapd(0,1,0,2)] - 2*lv[mapd(0,0,0,2)] + lv[mapd(0,-1,0,2)] )
        / ( parameters.geometry.dy * parameters.geometry.dy );
}

inline FLOAT d2wdz2 ( const FLOAT * const lv, const Parameters & parameters ) {
    return ( lv[mapd(0,0,1,2)] - 2*lv[mapd(0,0,0,2)] + lv[mapd(0,0,-1,2)] )
        / ( parameters.geometry.dz * parameters.geometry.dz );
}



inline FLOAT duvdx ( const FLOAT * const lv, const Parameters & parameters ) {
    return 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,0)] + lv [mapd(0,1,0,0)] ) *
                         ( lv [mapd(0,0,0,1)] + lv [mapd(1,0,0,1)] ) ) -
                       ( ( lv [mapd(-1,0,0,0)] + lv [mapd(-1,1,0,0)] ) *
                         ( lv [mapd(-1,0,0,1)] + lv [mapd(0,0,0,1)] ) ) ) +
      parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,0)] + lv [mapd(0,1,0,0)] ) *
                              ( lv [mapd(0,0,0,1)] - lv [mapd(1,0,0,1)] ) ) -
                       ( fabs ( lv [mapd(-1,0,0,0)] + lv [mapd(-1,1,0,0)] ) *
                              ( lv [mapd(-1,0,0,1)] - lv [mapd(0,0,0,1)] ) ) ) ) /
                       parameters.geometry.dx;
}

inline FLOAT duvdy ( const FLOAT * const lv, const Parameters & parameters ) {
    return 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,1)] + lv [mapd(1,0,0,1)] ) *
                         ( lv [mapd(0,0,0,0)] + lv [mapd(0,1,0,0)] ) ) -
                       ( ( lv [mapd(0,-1,0,1)] + lv [mapd(1,-1,0,1)] ) *
                         ( lv [mapd(0,-1,0,0)] + lv [mapd(0,0,0,0)] ) ) ) +
      parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,1)] + lv [mapd(1,0,0,1)] ) *
                              ( lv [mapd(0,0,0,0)] - lv [mapd(0,1,0,0)] ) ) -
                       ( fabs ( lv [mapd(0,-1,0,1)] + lv [mapd(1,-1,0,1)] ) *
                              ( lv [mapd(0,-1,0,0)] - lv [mapd(0,0,0,0)] ) ) ) ) /
                       parameters.geometry.dy;
}

inline FLOAT duwdx ( const FLOAT * const lv, const Parameters & parameters ) {
    return 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,0)] + lv [mapd(0,0,1,0)] ) *
                         ( lv [mapd(0,0,0,2)] + lv [mapd(1,0,0,2)] ) ) -
                       ( ( lv [mapd(-1,0,0,0)] + lv [mapd(-1,0,1,0)] ) *
                         ( lv [mapd(-1,0,0,2)] + lv [mapd(0,0,0,2)] ) ) ) +
      parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,0)] + lv [mapd(0,0,1,0)] ) *
                              ( lv [mapd(0,0,0,2)] - lv [mapd(1,0,0,2)] ) ) -
                       ( fabs ( lv [mapd(-1,0,0,0)] + lv [mapd(-1,0,1,0)] ) *
                              ( lv [mapd(-1,0,0,2)] - lv [mapd(0,0,0,2)] ) ) ) ) /
                       parameters.geometry.dx;
}

inline FLOAT duwdz ( const FLOAT * const lv, const Parameters & parameters ) {
    return 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,2)] + lv [mapd(1,0,0,2)] ) *
                         ( lv [mapd(0,0,0,0)] + lv [mapd(0,0,1,0)] ) ) -
                       ( ( lv [mapd(0,0,-1,2)] + lv [mapd(1,0,-1,2)] ) *
                         ( lv [mapd(0,0,-1,0)] + lv [mapd(0,0,0,0)] ) ) ) +
      parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,2)] + lv [mapd(1,0,0,2)] ) *
                              ( lv [mapd(0,0,0,0)] - lv [mapd(0,0,1,0)] ) ) -
                       ( fabs ( lv [mapd(0,0,-1,2)] + lv [mapd(1,0,-1,2)] ) *
                              ( lv [mapd(0,0,-1,0)] - lv [mapd(0,0,0,0)] ) ) ) ) /
                       parameters.geometry.dz;
}

inline FLOAT dvwdy ( const FLOAT * const lv, const Parameters & parameters ) {
    return 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,1)] + lv [mapd(0,0,1,1)] ) *
                         ( lv [mapd(0,0,0,2)] + lv [mapd(0,1,0,2)] ) ) -
                       ( ( lv [mapd(0,-1,0,1)] + lv [mapd(0,-1,1,1)] ) *
                         ( lv [mapd(0,-1,0,2)] + lv [mapd(0,0,0,2)] ) ) ) +
      parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,1)] + lv [mapd(0,0,1,1)] ) *
                              ( lv [mapd(0,0,0,2)] - lv [mapd(0,1,0,2)] ) ) -
                       ( fabs ( lv [mapd(0,-1,0,1)] + lv [mapd(0,-1,1,1)] ) *
                              ( lv [mapd(0,-1,0,2)] - lv [mapd(0,0,0,2)] ) ) ) ) /
                       parameters.geometry.dy;
}

inline FLOAT dvwdz ( const FLOAT * const lv, const Parameters & parameters ) {
    return 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,2)] + lv [mapd(0,1,0,2)] ) *
                         ( lv [mapd(0,0,0,1)] + lv [mapd(0,0,1,1)] ) ) -
                       ( ( lv [mapd(0,0,-1,2)] + lv [mapd(0,1,-1,2)] ) *
                         ( lv [mapd(0,0,-1,1)] + lv [mapd(0,0,0,1)] ) ) ) +
      parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,2)] + lv [mapd(1,0,0,2)] ) *
                              ( lv [mapd(0,0,0,1)] - lv [mapd(0,0,1,1)] ) ) -
                       ( fabs ( lv [mapd(0,0,-1,2)] + lv [mapd(0,1,-1,2)] ) *
                              ( lv [mapd(0,0,-1,1)] - lv [mapd(0,0,0,1)] ) ) ) ) /
                       parameters.geometry.dz;
}

inline FLOAT du2dx ( const FLOAT * const lv, const Parameters & parameters ) {
    return 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,0)] + lv [mapd(1,0,0,0)] ) *
                         ( lv [mapd(0,0,0,0)] + lv [mapd(1,0,0,0)] ) ) -
                       ( ( lv [mapd(-1,0,0,0)] + lv [mapd(0,0,0,0)] ) *
                         ( lv [mapd(-1,0,0,0)] + lv [mapd(0,0,0,0)] ) ) ) +
      parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,0)] + lv [mapd(1,0,0,0)] ) *
                              ( lv [mapd(0,0,0,0)] - lv [mapd(1,0,0,0)] ) ) -
                       ( fabs ( lv [mapd(-1,0,0,0)] + lv [mapd(0,0,0,0)] ) *
                              ( lv [mapd(-1,0,0,0)] - lv [mapd(0,0,0,0)] ) ) ) ) /
                       parameters.geometry.dx;
}

inline FLOAT dv2dy ( const FLOAT * const lv, const Parameters & parameters ) {
    return 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,1)] + lv [mapd(0,1,0,1)] ) *
                         ( lv [mapd(0,0,0,1)] + lv [mapd(0,1,0,1)] ) ) -
                       ( ( lv [mapd(0,-1,0,1)] + lv [mapd(0,0,0,1)] ) *
                         ( lv [mapd(0,-1,0,1)] + lv [mapd(0,0,0,1)] ) ) ) +
      parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,1)] + lv [mapd(0,1,0,1)] ) *
                              ( lv [mapd(0,0,0,1)] - lv [mapd(0,1,0,1)] ) ) -
                       ( fabs ( lv [mapd(0,-1,0,1)] + lv [mapd(0,0,0,1)] ) *
                              ( lv [mapd(0,-1,0,1)] - lv [mapd(0,0,0,1)] ) ) ) ) /
                       parameters.geometry.dy;
}

inline FLOAT dw2dz ( const FLOAT * const lv, const Parameters & parameters ) {
    return 1.0 /4.0 * ( ( ( ( lv [mapd(0,0,0,2)] + lv [mapd(0,0,1,2)] ) *
                         ( lv [mapd(0,0,0,2)] + lv [mapd(0,0,1,2)] ) ) -
                       ( ( lv [mapd(0,0,-1,2)] + lv [mapd(0,0,0,2)] ) *
                         ( lv [mapd(0,0,-1,2)] + lv [mapd(0,0,0,2)] ) ) ) +
      parameters.solver.gamma * ( ( fabs ( lv [mapd(0,0,0,2)] + lv [mapd(0,0,1,2)] ) *
                              ( lv [mapd(0,0,0,2)] - lv [mapd(0,0,1,2)] ) ) -
                       ( fabs ( lv [mapd(0,0,-1,2)] + lv [mapd(0,0,0,2)] ) *
                              ( lv [mapd(0,0,-1,2)] - lv [mapd(0,0,0,2)] ) ) ) ) /
                       parameters.geometry.dz;
}


inline FLOAT computeF2D(const FLOAT * const localVelocity, const Parameters & parameters, FLOAT dt){
    return localVelocity [mapd(0,0,0,0)]
        + dt * ( 1 / parameters.flow.Re * ( d2udx2 ( localVelocity, parameters )
                    + d2udy2(localVelocity, parameters)) - du2dx (localVelocity, parameters)
                    - duvdy (localVelocity, parameters) + parameters.environment.gx);
}

inline FLOAT computeG2D(const FLOAT * const localVelocity, const Parameters & parameters, FLOAT dt){
    return localVelocity [mapd(0,0,0,1)] 
        + dt * ( 1 / parameters.flow.Re * ( d2vdx2 ( localVelocity, parameters )
                    + d2vdy2(localVelocity, parameters)) - duvdx (localVelocity, parameters)
                    - dv2dy (localVelocity, parameters) + parameters.environment.gy);
}


inline FLOAT computeF3D(const FLOAT * const localVelocity, const Parameters & parameters, FLOAT dt){
    return localVelocity [mapd(0,0,0,0)]
                +  dt * ( 1 / parameters.flow.Re * ( d2udx2 ( localVelocity, parameters )
                + d2udy2 ( localVelocity, parameters ) + d2udz2 ( localVelocity, parameters ) )
                - du2dx ( localVelocity, parameters ) - duvdy ( localVelocity, parameters )
                - duwdz ( localVelocity, parameters ) + parameters.environment.gx );
}


inline FLOAT computeG3D(const FLOAT * const localVelocity, const Parameters & parameters, FLOAT dt){
    return localVelocity [mapd(0,0,0,1)]
                +  dt * ( 1 / parameters.flow.Re * ( d2vdx2 ( localVelocity, parameters )
                + d2vdy2 ( localVelocity, parameters ) + d2vdz2 ( localVelocity, parameters ) )
                - dv2dy ( localVelocity, parameters ) - duvdx ( localVelocity, parameters )
                - dvwdz ( localVelocity, parameters ) + parameters.environment.gy );
}


inline FLOAT computeH3D(const FLOAT * const localVelocity, const Parameters & parameters, FLOAT dt){
    return localVelocity [mapd(0,0,0,2)]
                +  dt * ( 1 / parameters.flow.Re * ( d2wdx2 ( localVelocity, parameters )
                + d2wdy2 ( localVelocity, parameters ) + d2wdz2 ( localVelocity, parameters ) )
                - dw2dz ( localVelocity, parameters ) - duwdx ( localVelocity, parameters )
                - dvwdy ( localVelocity, parameters ) + parameters.environment.gz );
}

#endif
