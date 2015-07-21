#ifndef __SOR_H_
#define __SOR_H_


#include "helper.h"

/**
 * One GS iteration for the pressure Poisson equation. Besides, the routine must
 * also set the boundary values for P according to the specification. The
 * residual for the termination criteria has to be stored in res.
 *
 * An \omega = 1 GS - implementation is given within sor.c.
 */
int sor(
  double omg,
  double dx,
  double dy,
  double dz,
  int    imax,
  int    jmax,
  int    kmax,
  double ***P,
  double ***RS,
  double *res,
  int ***Flag,
  struct p_pointer *PP1,
	int FluidCells
/*  double presLeft,
  double presRight*/
);


#endif
