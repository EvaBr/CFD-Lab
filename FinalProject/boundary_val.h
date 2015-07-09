#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/*
 * The boundary values of the problem are set.
 */

void boundaryvalues_no_slip(
  int i,
  int j,
  int k,
  double ***U,
  double ***V,
  double ***W,
  int ***Flag
);
void boundaryvalues_free_slip(
  int i,
  int j,
  int k,
  double ***U,
  double ***V,
  double ***W,
  int ***Flag
);


void boundaryvalues_outflow(
    int i,
    int j,
    int k,
    double ***U,
    double ***V,
    double ***W,
    int ***Flag
);

void boundaryvalues(
  int imax,
  int jmax,
  int kmax,
  double ***U,
  double ***V,
  double ***W,
  double ***P,
  int wl,
  int wr,
  int wf,
  int wh,
  int wt,
  int wb,
  double ***F,
  double ***G,
  double ***H,
  char *problem,
  int ***Flag,
  double velIN,
  double velMW
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
