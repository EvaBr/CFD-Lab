#include "sor.h"
#include <math.h>
#include "boundary_val.h"
#include "helper.h"




void sor(
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
		int ***Flag
		/*double presLeft,
  double presRight*/
) {
	int i,j,k;
	double rloc;
	double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)+1.0/(dz*dz)));

	double FluidCells = 0.0;
	double tmp = 0;
	double dx2 = dx*dx;
	double dy2 = dy*dy;
	double dz2 = dz*dz;

	double omg2 = 1-omg;

	/* SOR iteration */

	for(i = 1; i <= imax; i++) {
		for(j = 1; j<=jmax; j++) {
			for(k = 1; k<=kmax; k++) {
				if(isfluid(Flag[i][j][k]) && !emptyneighbor(Flag[i][j][k])){


					P[i][j][k] = omg2*P[i][j][k]+
					                coeff*(
					                      (P[i+1][j  ][k  ]+P[i-1][j  ][k  ])/dx2 +
					                      (P[i  ][j+1][k  ]+P[i  ][j-1][k  ])/dy2 +
					                      (P[i  ][j  ][k+1]+P[i  ][j  ][k-1])/dz2

					                    - RS[i  ][j  ][k  ]
					                       );


					FluidCells++;
				}
			}
		}
	}


	/* compute the residual */
	rloc = 0;
	for(i = 1; i <= imax; i++) {
		for(j = 1; j <= jmax; j++) {
			for(k = 1; k <= kmax; k++){
				if(isfluid(Flag[i][j][k])&& !emptyneighbor(Flag[i][j][k])){
					tmp =  (P[i+1][j][k]-2.0*P[i][j][k]+P[i-1][j][k])/(dx2) + (P[i][j+1][k]-2.0*P[i][j][k]+P[i][j-1][k])/(dy2) + (P[i][j][k+1]-2.0*P[i][j][k]+P[i][j][k-1])/(dz2) - RS[i][j][k];
					rloc += tmp*tmp;

				}
			}
		}
	}
	rloc = rloc/FluidCells;
	rloc = sqrt(rloc);
	/* set residual */
	*res = rloc;



	boundaryvalues_pressure(P,Flag,imax,jmax,kmax);


	//check if given pressure, set Dirichlet or Neuman BC according to that, for left and right wall
	/*for (j=0; j<=jmax+1; j++) {
	  if((Flag[0][j] & 32) != 0) {//pressure given -> overwrite with Dirichlet BC
	        P[0][j] = presLeft*2.0 - P[1][j];
	  }
	  if((Flag[imax+1][j] & 32) != 0) {//pressure given -> overwrite with Dirichlet BC
 	        P[imax+1][j] = 2.0*presRight - P[imax][j];
          }
  }*/
}
