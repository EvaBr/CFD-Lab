#include "sor.h"
#include <math.h>
#include "boundary_val.h"
#include "helper.h"
#include <omp.h>

#define SOR_ITER_MAX 3


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
		/*double presLeft,
  double presRight*/
) {
	int i,j,k;
	int s;
	double rloc;
	double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)+1.0/(dz*dz)));

	double tmp = 0;
	double dx2 = 1.0/(dx*dx);
	double dy2 = 1.0/(dy*dy);
	double dz2 = 1.0/(dz*dz);

	double omg2 = 1-omg;
	int sor_iter = 0;
	/* SOR iteration */


	for(sor_iter = 1; sor_iter <= SOR_ITER_MAX; sor_iter++) {
		#pragma omp parallel for private(i,j,k)
		for(s=0;s<FluidCells;s++){
			struct p_pointer *pp = &PP1[s];
			if(pp->p){
				i = pp->i;
				j = pp->j;
				k = pp->k;
				P[i][j][k] = omg2*P[i][j][k]+
						coeff*(
								(P[i+1][j  ][k  ]+P[i-1][j  ][k  ])*dx2 +
								(P[i  ][j+1][k  ]+P[i  ][j-1][k  ])*dy2 +
								(P[i  ][j  ][k+1]+P[i  ][j  ][k-1])*dz2

								- RS[i  ][j  ][k  ]
						);
			}
		}


		if(sor_iter==SOR_ITER_MAX){
			rloc = 0;
		#pragma omp parallel for private(i,j,k,tmp)
			for(s=0;s<FluidCells;s++){
				struct p_pointer *pp = &PP1[s];
				i = pp->i;
				j = pp->j;
				k = pp->k;
				if(isfluid(Flag[i][j][k]) && !emptyneighbor(Flag[i][j][k])){
					tmp =  (P[i+1][j][k]-2.0*P[i][j][k]+P[i-1][j][k])*dx2 + (P[i][j+1][k]-2.0*P[i][j][k]+P[i][j-1][k])*dy2 + (P[i][j][k+1]-2.0*P[i][j][k]+P[i][j][k-1])*dz2 - RS[i][j][k];
					tmp = tmp*tmp;
					rloc += tmp;
				}

			}

			rloc = rloc/FluidCells;
			rloc = sqrt(rloc);
			/* set residual */
			*res = rloc;
		}

		boundaryvalues_pressure(P,Flag,imax,jmax,kmax);

	}


	/* compute the residual */





	return SOR_ITER_MAX;

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
