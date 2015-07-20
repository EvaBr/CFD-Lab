#include "BFStepInitStencil.h"

BFStepInitStencil::BFStepInitStencil (const Parameters & parameters) :
    FieldStencil<FlowField> (parameters), xLimit (parameters.bfStep.width), yLimit (parameters.bfStep.height)
{}

void BFStepInitStencil::apply(FlowField & flowField, int i, int j){
    IntScalarField & flags = flowField.getFlags();
    // DMITRIIS VERSION : Comment all subsequent checks and set flag to zero
    if (i + _parameters.parallel.firstCorner[0] - 2 < xLimit &&
        j + _parameters.parallel.firstCorner[1] - 2 < yLimit){
        flags.getValue(i, j) = OBSTACLE_SELF;
    }
    if (i + _parameters.parallel.firstCorner[0] - 2-1 < xLimit &&
        j + _parameters.parallel.firstCorner[1] - 2 < yLimit){
        flags.getValue(i, j) += OBSTACLE_LEFT;
    }
    if (i + _parameters.parallel.firstCorner[0] - 2+1 < xLimit &&
        j + _parameters.parallel.firstCorner[1] - 2 < yLimit){
        flags.getValue(i, j) += OBSTACLE_RIGHT;
    }
    if (i + _parameters.parallel.firstCorner[0] - 2 < xLimit &&
        j + _parameters.parallel.firstCorner[1] - 2-1 < yLimit){
        flags.getValue(i, j) += OBSTACLE_BOTTOM;
    }
    if (i + _parameters.parallel.firstCorner[0] - 2 < xLimit &&
        j + _parameters.parallel.firstCorner[1] - 2+1 < yLimit){
        flags.getValue(i, j) += OBSTACLE_TOP;
    }
}

void BFStepInitStencil::apply(FlowField & flowField, int i, int j, int k){
    IntScalarField & flags = flowField.getFlags();
    if (i + _parameters.parallel.firstCorner[0] - 2 < xLimit &&
        j + _parameters.parallel.firstCorner[1] - 2 < yLimit){
        flags.getValue(i, j, k) = OBSTACLE_SELF;

        // The obstacle is 2D, so we can say the following
        flags.getValue(i, j, k) += OBSTACLE_FRONT;
        flags.getValue(i, j, k) += OBSTACLE_BACK;
    }
    if (i + _parameters.parallel.firstCorner[0] - 2-1 < xLimit &&
        j + _parameters.parallel.firstCorner[1] - 2 < yLimit){
        flags.getValue(i, j, k) += OBSTACLE_LEFT;
    }
    if (i + _parameters.parallel.firstCorner[0] - 2+1 < xLimit &&
        j + _parameters.parallel.firstCorner[1] - 2 < yLimit){
        flags.getValue(i, j, k) += OBSTACLE_RIGHT;
    }
    if (i + _parameters.parallel.firstCorner[0] - 2 < xLimit &&
        j + _parameters.parallel.firstCorner[1] - 2-1 < yLimit){
        flags.getValue(i, j, k) += OBSTACLE_BOTTOM;
    }
    if (i + _parameters.parallel.firstCorner[0] - 2 < xLimit &&
        j + _parameters.parallel.firstCorner[1] - 2+1 < yLimit){
        flags.getValue(i, j, k) += OBSTACLE_TOP;
    }
}
