#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/*
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V,
  double **P,
  int wl,
  int wr,
  int wt,
  int wb,
  double **F,
  double **G,
  char *problem,
  int **Flag,
  double UI
);

void spec_boundary_val(
  char *problem,
  int imax,
  int jmax,
  double **U,
  double **V,
  int **Flag,
  double UI
);

#endif
