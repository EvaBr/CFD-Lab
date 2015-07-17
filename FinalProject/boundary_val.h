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
void boundaryvalues_moving_wall(
		int i,
		int j,
		int k,
		double ***U,
		double ***V,
		double ***W,
		int ***Flag,
		double *velMW
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
void boundaryvalues_inflow(
		int i,
		int j,
		int k,
		double ***U,
		double ***V,
		double ***W,
		int ***Flag,
		double velIN
);

void boundaryvalues_pressure(
		double ***P,
		int ***Flag,
		int imax,
		int jmax,
		int kmax
);


void boundaryvalues(
		int imax,
		int jmax,
		int kmax,
		double ***U,
		double ***V,
		double ***W,
		double ***P,
		double ***F,
		double ***G,
		double ***H,
		char *problem,
		int ***Flag,
		double velIN,
		double *velMW
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
