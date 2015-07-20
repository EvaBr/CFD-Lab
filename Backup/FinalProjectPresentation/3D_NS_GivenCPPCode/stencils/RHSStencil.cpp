#include "RHSStencil.h"

RHSStencil::RHSStencil ( const Parameters & parameters ) : FieldStencil<FlowField> ( parameters ) {}


void RHSStencil::apply ( FlowField & flowField, int i, int j ) {
    flowField.getRHS().getScalar (i, j) = 1.0 / _parameters.timestep.dt *
        ( ( flowField.getFGH().getVector(i, j)[0] - flowField.getFGH().getVector(i-1, j)[0] )
          / _parameters.geometry.dx
        + ( flowField.getFGH().getVector(i, j)[1] - flowField.getFGH().getVector(i, j-1)[1] )
          / _parameters.geometry.dy );
}


void RHSStencil::apply ( FlowField & flowField, int i, int j, int k ) {
    flowField.getRHS().getScalar (i, j, k) = 1.0 / _parameters.timestep.dt *
        ( (flowField.getFGH().getVector(i, j, k)[0] - flowField.getFGH().getVector(i-1, j, k)[0])
          / _parameters.geometry.dx
        + (flowField.getFGH().getVector(i, j, k)[1] - flowField.getFGH().getVector(i, j-1, k)[1])
          / _parameters.geometry.dy
        + (flowField.getFGH().getVector(i, j, k)[2] - flowField.getFGH().getVector(i, j, k-1)[2])
          / _parameters.geometry.dz );
}
