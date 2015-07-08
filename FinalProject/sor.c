#include "sor.h"
#include <math.h>
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

  /* SOR iteration */
  for(i = 1; i <= imax; i++) {
    for(j = 1; j<=jmax; j++) {
      for(k = 1; k<=kmax; k++) {
        if(isfluid(Flag[i][j][k])){
	         P[i][j][k] = (1.0-omg)*P[i][j][k]
        	      + coeff*( (P[i+1][j][k]+P[i-1][j][k])/(dx*dx) + (P[i][j+1][k]+P[i][j-1][k])/(dy*dy) + (P[i][j][k+1]-P[i][j][k])/(dz*dz) - RS[i][j][k]);
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
        if(isfluid(Flag[i][j][k])){
	         rloc += ( (P[i+1][j][k]-2.0*P[i][j][k]+P[i-1][j][k])/(dx*dx) + (P[i][j+1][k]-2.0*P[i][j][k]+P[i][j-1][k])/(dy*dy) + (P[i][j][k+1]-2.0*P[i][j][k]+P[i][j][k-1])/(dz*dz) - RS[i][j][k])*
        	      ( (P[i+1][j][k]-2.0*P[i][j][k]+P[i-1][j][k])/(dx*dx) + (P[i][j+1][k]-2.0*P[i][j][k]+P[i][j-1][k])/(dy*dy) + (P[i][j][k+1]-2.0*P[i][j][k]+P[i][j][k-1])/(dz*dz) - RS[i][j][k]);
        }
      }
    }
  }
  rloc = rloc/FluidCells;
  rloc = sqrt(rloc);
  /* set residual */
  *res = rloc;


  /* set boundary values, here just for the 'real' boundaries - no air included yet (if even needed?) */
  for(i = 0; i <= imax+1; i++) {
    //P[i][0] = P[i][1];		//boundary cond at lower and upper wall
    //P[i][jmax+1] = P[i][jmax];
    for(j = 0; j <= jmax+1; j++) {
      for(k = 0; k <= kmax+1; k++) {
	       switch(getcelltype(Flag[i][j][k])){
		        case B_N: P[i][j][k] = P[i][j+1][k]; break;
		        case B_O: P[i][j][k] = P[i+1][j][k]; break;
		        case B_S: P[i][j][k] = P[i][j-1][k]; break;
		        case B_W: P[i][j][k] = P[i-1][j][k]; break;
            case B_U: P[i][j][k] = P[i][j][k+1]; break;
		        case B_D: P[i][j][k] = P[i][j][k-1]; break;

		        case B_NO: P[i][j][k] = (P[i+1][j][k] + P[i][j+1][k])*0.5; break;
		        case B_NW: P[i][j][k] = (P[i-1][j][k] + P[i][j+1][k])*0.5; break;
		        case B_SO: P[i][j][k] = (P[i+1][j][k] + P[i][j-1][k])*0.5; break;
		        case B_SW: P[i][j][k] = (P[i-1][j][k] + P[i][j-1][k])*0.5; break;
            case B_NU: P[i][j][k] = (P[i][j][k+1] + P[i][j+1][k])*0.5; break;
		        case B_ND: P[i][j][k] = (P[i][j][k-1] + P[i][j+1][k])*0.5; break;
		        case B_SU: P[i][j][k] = (P[i][j][k+1] + P[i][j-1][k])*0.5; break;
		        case B_SD: P[i][j][k] = (P[i][j][k-1] + P[i][j-1][k])*0.5; break;
            case B_OU: P[i][j][k] = (P[i+1][j][k] + P[i][j][k+1])*0.5; break;
		        case B_WU: P[i][j][k] = (P[i-1][j][k] + P[i][j][k+1])*0.5; break;
		        case B_OD: P[i][j][k] = (P[i+1][j][k] + P[i][j][k-1])*0.5; break;
		        case B_WD: P[i][j][k] = (P[i-1][j][k] + P[i][j][k-1])*0.5; break;

            case B_NOU: P[i][j][k] = (P[i+1][j][k] + P[i][j+1][k] + P[i][j][k+1])*1.0/3.0; break;
		        case B_NWU: P[i][j][k] = (P[i-1][j][k] + P[i][j+1][k] + P[i][j][k+1])*1.0/3.0; break;
		        case B_SOU: P[i][j][k] = (P[i+1][j][k] + P[i][j-1][k] + P[i][j][k+1])*1.0/3.0; break;
		        case B_SWU: P[i][j][k] = (P[i-1][j][k] + P[i][j-1][k] + P[i][j][k+1])*1.0/3.0; break;
            case B_NOD: P[i][j][k] = (P[i+1][j][k] + P[i][j+1][k] + P[i][j][k-1])*1.0/3.0; break;
		        case B_NWD: P[i][j][k] = (P[i-1][j][k] + P[i][j+1][k] + P[i][j][k-1])*1.0/3.0; break;
		        case B_SOD: P[i][j][k] = (P[i+1][j][k] + P[i][j-1][k] + P[i][j][k-1])*1.0/3.0; break;
		        case B_SWD: P[i][j][k] = (P[i-1][j][k] + P[i][j-1][k] + P[i][j][k-1])*1.0/3.0; break;

		        //case C_B: P[i][j] = 0; break;
	          default: break;
          }
      }
    }
  }


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
