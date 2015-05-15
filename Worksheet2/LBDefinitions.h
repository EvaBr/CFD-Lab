#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

#include <math.h>

  static const int LATTICEVELOCITIES[19][3] = {{0, -1.0, -1.0}, {-1.0, 0, -1.0}, {0, 0, -1.0}, {1.0, 0, -1.0}, {0, 1.0, -1.0}, {-1.0, -1.0, 0}, {0, -1.0, 0},
						{1.0, -1.0, 0}, {-1.0, 0, 0}, {0, 0, 0}, {1.0, 0, 0}, {-1.0, 1.0, 0}, {0, 1.0, 0}, {1.0, 1.0, 0},
						{0, -1.0, 1.0}, {-1.0, 0, 1.0}, {0, 0, 1.0}, {1.0, 0, 1.0}, {0, 1.0, 1.0}};
  static const double LATTICEWEIGHTS[19] = {0.027778, 0.027778, 0.055556, 0.027778, 0.027778, 0.027778, 0.055556, 0.027778, 0.055556, 0.33333,
						0.055556, 0.027778, 0.055556, 0.027778, 0.027778, 0.027778, 0.055556, 0.027778, 0.027778};
  static const double C_S = 0.57735;



  static const int FLUID = 0;
  static const int NO_SLIP = 1;
  static const int MOVING_WALL = 2;

#endif

