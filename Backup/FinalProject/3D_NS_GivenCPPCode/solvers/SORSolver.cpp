#include <float.h> // To get the max double
#include <math.h>
#include "SORSolver.h"

SORSolver::SORSolver(FlowField & flowField, const Parameters & parameters):
    LinearSolver(flowField, parameters){}

void SORSolver::solve() {
    // Solves the heat equation
    FLOAT resnorm = DBL_MAX, tol = 1e-4;

    double omg = 1.7;
    int iterations = -1;
    int it = 0;

    int nx = _flowField.getNx(), ny = _flowField.getNy(), nz = _flowField.getNz();
    FLOAT dx = _parameters.geometry.dx, dy = _parameters.geometry.dy, dz = _parameters.geometry.dz;
    ScalarField & P = _flowField.getPressure();
    if (_parameters.geometry.dim == 3){
        do {
            for (int k = 2; k < nz + 2; k++){
                for (int j = 2; j < ny + 2; j++){
                    for (int i = 2; i < nx + 2; i++){
                        P.getScalar(i,j,k) =
                            (((P.getScalar(i-1,j,k) + P.getScalar(i+1,j,k)) / (dx*dx)) +
                             ((P.getScalar(i,j-1,k) + P.getScalar(i,j+1,k)) / (dy*dy)) +
                             ((P.getScalar(i,j,k-1) + P.getScalar(i,j,k+1)) / (dz*dz)) -
                             _flowField.getRHS().getScalar(i,j,k))
                            * omg / (2 * (1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz))) + (1 - omg) *
                            P.getScalar(i,j,k);
                    }
                }
            }

            for ( int j = 2; j < ny + 2; j++ ) {
                for (int k = 2; k < nz + 2; k++ ) {
                    P.getScalar(1,j,k) = P.getScalar(2,j,k);
                    P.getScalar(nx+2,j,k) = P.getScalar(nx+1,j,k);
                }
            }

            for ( int i = 2; i < nx + 2; i++ ) {
                for (int k = 2; k < nz + 2; k++ ) {
                    P.getScalar(i,1,k) = P.getScalar(i,2,k);
                    P.getScalar(i,ny+2,k) = P.getScalar(i,ny+1,k);
                }
            }

            for ( int i = 2; i < nx + 2; i++ ) {
                for (int j = 2; j < ny + 2; j++) {
                    P.getScalar(i,j,1) = P.getScalar(i,j,2);
                    P.getScalar(i,j,nz+2) = P.getScalar(i,j,nz+1);
                }
            }

            resnorm = 0;
            for (int k = 2; k < nz + 2; k++){
                for (int j = 2; j < ny + 2; j++){
                    for (int i = 2; i < nx + 2; i++){
                        resnorm += pow(fabs(
                                    (P.getScalar(i-1,j,k) - 2*P.getScalar(i,j,k) +
                                     P.getScalar(i+1,j,k)) / (dx * dx) +
                                    (P.getScalar(i,j-1,k) - 2*P.getScalar(i,j,k) +
                                     P.getScalar(i,j+1,k)) / (dy * dy) +
                                    (P.getScalar(i,j,k-1) - 2*P.getScalar(i,j,k) +
                                     P.getScalar(i,j,k+1)) / (dz * dz) -
                                    _flowField.getRHS().getScalar(i,j,k)
                                    ),2);
                    }
                }
            }
            resnorm = sqrt (resnorm / (nx * ny * nz));
            // std::cout << "Residual norm: " << resnorm << std::endl;

            it++;
            iterations--;
        } while ( resnorm > tol && iterations);

    } if (_parameters.geometry.dim == 2){
        do {
            for (int j = 2; j < ny + 2; j++){
                for (int i = 2; i < nx + 2; i++){
                    P.getScalar(i,j) =
                        (((P.getScalar(i-1,j) + P.getScalar(i+1,j)) / (dx*dx)) +
                         ((P.getScalar(i,j-1) + P.getScalar(i,j+1)) / (dy*dy)) -
                         _flowField.getRHS().getScalar(i,j))
                        * omg / (2 * (1/(dx*dx) + 1/(dy*dy))) + (1 - omg) * P.getScalar(i,j);
                }
            }

            resnorm = 0;
            for (int j = 2; j < ny + 2; j++){
                for (int i = 2; i < nx + 2; i++){
                    resnorm += ((P.getScalar(i-1,j) - 2*P.getScalar(i,j) +
                                 P.getScalar(i+1,j)) / (dx * dx) +
                                (P.getScalar(i,j-1) - 2*P.getScalar(i,j) +
                                 P.getScalar(i,j+1)) / (dy * dy) -
                                _flowField.getRHS().getScalar(i,j)) *
                               ((P.getScalar(i-1,j) - 2*P.getScalar(i,j) +
                                 P.getScalar(i+1,j)) / (dx * dx) +
                                (P.getScalar(i,j-1) - 2*P.getScalar(i,j) +
                                 P.getScalar(i,j+1)) / (dy * dy) -
                                _flowField.getRHS().getScalar(i,j));
                }
            }
            resnorm = sqrt (resnorm / (nx * ny));

            for ( int j = 2; j < ny + 2; j++ ) {
                P.getScalar(1,j) = P.getScalar(2,j);
                P.getScalar(nx+2,j) = P.getScalar(nx+1,j);
            }

            for ( int i = 2; i < nx + 2; i++ ) {
                P.getScalar(i,1) = P.getScalar(i,2);
                P.getScalar(i,ny+2) = P.getScalar(i,ny+1);
            }

            // std::cout << "Residual norm: " << resnorm << std::endl;
            iterations--;
            it++;
        } while ( resnorm > tol && iterations);
    }
}
