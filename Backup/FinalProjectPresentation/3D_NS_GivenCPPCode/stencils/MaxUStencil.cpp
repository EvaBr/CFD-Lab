#include "MaxUStencil.h"
#include <algorithm>
#include <math.h>


MaxUStencil::MaxUStencil (const Parameters & parameters) :
    FieldStencil<FlowField> (parameters), BoundaryStencil<FlowField> (parameters) {
    reset();
}

void MaxUStencil::apply (FlowField & flowField, int i, int j){
    cellMaxValue(flowField, i, j);
}

void MaxUStencil::apply (FlowField & flowField, int i, int j, int k){
    cellMaxValue(flowField, i, j, k);
}

void MaxUStencil::applyLeftWall   ( FlowField & flowField, int i, int j ){
    cellMaxValue(flowField, i, j);
}

void MaxUStencil::applyRightWall  ( FlowField & flowField, int i, int j ){
    cellMaxValue(flowField, i, j);
}

void MaxUStencil::applyBottomWall ( FlowField & flowField, int i, int j ){
    cellMaxValue(flowField, i, j);
}

void MaxUStencil::applyTopWall    ( FlowField & flowField, int i, int j ){
    cellMaxValue(flowField, i, j);
}

void MaxUStencil::applyLeftWall   ( FlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}

void MaxUStencil::applyRightWall  ( FlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}

void MaxUStencil::applyBottomWall ( FlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}

void MaxUStencil::applyTopWall    ( FlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}

void MaxUStencil::applyFrontWall  ( FlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}

void MaxUStencil::applyBackWall   ( FlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}


void MaxUStencil::cellMaxValue(FlowField & flowField, int i, int j){
    FLOAT * velocity = flowField.getVelocity().getVector(i, j);
    if (fabs(velocity[0]) > _maxValues[0]){
        _maxValues[0] = fabs(velocity[0]);
    }
    if (fabs(velocity[1]) > _maxValues[1]){
        _maxValues[1] = fabs(velocity[1]);
    }
}

void MaxUStencil::cellMaxValue(FlowField & flowField, int i, int j, int k){
    FLOAT * velocity = flowField.getVelocity().getVector(i, j, k);
    if (fabs(velocity[0]) > _maxValues[0]){
        _maxValues[0] = fabs(velocity[0]);
    }
    if (fabs(velocity[1]) > _maxValues[1]){
        _maxValues[1] = fabs(velocity[1]);
    }
    if (fabs(velocity[2]) > _maxValues[2]){
        _maxValues[2] = fabs(velocity[2]);
    }
}

void MaxUStencil::reset () {
    _maxValues[0] = 0;
    _maxValues[1] = 0;
    _maxValues[2] = 0;
}

const FLOAT * MaxUStencil::getMaxValues() const{
    return _maxValues;
}
