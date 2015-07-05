#include "helper.h"
#include <math.h>

void calculate_fg(
  double Re,
  double GX,
  double GY,
  double GZ,
  double alpha,
  double dt,
  double dx,
  double dy,
  double dz,
  int imax,
  int jmax,
  int kmax,
  double ***U,
  double ***V,
  double ***W,
  double ***F,
  double ***G,
  double ***H,
  int ***Flag
){
/*	double uij, vij, uuj, udj, uiu, uid, viu, vid, vuj, vdj, udu, vud;
	double Dx = 1/dx, Dy = 1/dy, Dz = 1/dz;
	int i, j, k;
	for (i=1; i<imax+1; i++){
		for  (j=1; j<jmax+1; j++){
			switch (Flag[i][j][k]){
				case (C_F):
					uij = U[i][j];
					uuj = U[i+1][j];
					udj = U[i-1][j];
					uiu = U[i][j+1];
					uid = U[i][j-1];
					udu = U[i-1][j+1];

					vij = V[i][j];
					viu = V[i][j+1];
					vid = V[i][j-1];
					vuj = V[i+1][j];
					vdj = V[i-1][j];
					vud = V[i+1][j-1];

					if (Flag[i+1][j] == C_F){ // add this if-loop in double for-loop to ensure calculation is only on edges
							// separating two fluid cells.

						F[i][j] = uij + dt*(1/Re*((uuj-2*uij+udj)*Dx*Dx + (uiu-2*uij+uid)*Dy*Dy) -
							0.25*Dx*(pow((uij+uuj),2)-pow((udj+uij),2) + alpha*(fabs(uij+uuj)*(uij-uuj)-fabs(udj+uij)*(udj-uij))) -
							0.25*Dy*((vij+vuj)*(uij+uiu)-(vid+vud)*(uid+uij) + alpha*(fabs(vij+vuj)*(uij-uiu)-fabs(vid+vud)*(uid-uij))) +
							GX);
					}
					if (Flag[i][j+1] == C_F){
						G[i][j] = vij + dt*(1/Re*((vuj-2*vij+vdj)*Dx*Dx + (viu-2*vij+vid)*Dy*Dy) -
							0.25*Dy*(pow((vij+viu),2)-pow((vid+vij),2) + alpha*(fabs(vij+viu)*(vij-viu)-fabs(vid+vij)*(vid-vij))) -
							0.25*Dx*((uij+uiu)*(vij+vuj)-(udj+udu)*(vdj+vij) + alpha*(fabs(uij+uiu)*(vij-vuj)-fabs(udj+udu)*(vdj-vij))) +
							GY);
					break;
				case (B_N):
					G[i][j] = V[i][j]; break;
				case (B_W):
					F[i-1][j] = U[i-1][j]; break;
				case (B_O):
					F[i][j] = U[i][j]; break;
				case (B_S):
					G[i][j-1] = V[i][j-1]; break;
				case (B_NO):
					G[i][j] = V[i][j];
					F[i][j] = U[i][j];
					break;
				case (B_NW):
					G[i][j] = V[i][j];
					F[i-1][j] = U[i-1][j];
					break;
				case (B_SO):
					G[i][j-1] = V[i][j-1];
					F[i][j] = U[i][j];
					break;
				case (B_SW):
					F[i-1][j] = U[i-1][j];
					G[i][j-1] = V[i][j-1];
					break;
				}
			}
		}*/
		/* rewrite G(i,0) and G(i, jmax) with bound.cond. for G */
//		G[i][0] = V[i][0];
//		G[i][jmax] = V[i][jmax];
//	}
//	for (j=1; j<jmax+1; j++){
		/* rewrite F(0,j) and F(imax, j) with bound.cond. for F */
//		F[0][j] = U[0][j];
//		F[imax][j] = U[imax][j];
//	}
}

void calculate_rs(
  double dt,
  double dx,
  double dy,
  double dz,
  int imax,
  int jmax,
  int kmax,
  double ***F,
  double ***G,
  double ***H,
  double ***RS
){
	int i, j, k;
	/*range of indices {1:imax}x{1:jmax} for RS*/
	for(i=1; i<imax+1; i++){
		for(j=1; j<jmax+1; j++){
      for(k=0; k<kmax+1; k++){
			  RS[i][j][k] = ( (F[i][j][k] - F[i-1][j][k])/dx + (G[i][j][k] - G[i][j-1][k])/dy + (H[i][j][k] - H[i][j][k-1])/dz ) / dt;
		  }
	  }
  }
}

void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  double dz,
  int imax,
  int jmax,
  int kmax,
  double ***U,
  double ***V,
  double ***W
){
	double tmp, maxi1, maxi2, maxi3;
	/* we rewrite dt only if tau is nonnegative, else we do nothing */
	if(tau>=0){
		tmp = 0.5*Re*(dx*dx*dy*dy*dz*dz)/(dx*dx+dy*dy+dz*dz);
		maxi1 = tmax(U, imax, jmax, kmax);
		maxi2 = tmax(V, imax, jmax, kmax);
    maxi3 = tmax(W, imax, jmax, kmax);
		*dt = tau * fmin(tmp, fmin(dx/maxi1, fmin(dy/maxi2, dz/maxi3)));
 	}
}


void calculate_uv(
  double dt,
  double dx,
  double dy,
  double dz,
  int imax,
  int jmax,
  int kmax,
  double ***U,
  double ***V,
  double ***W,
  double ***F,
  double ***G,
  double ***H,
  double ***P,
  int ***Flag
){
	int i, j, k;
	/*range of indices {1:(imax-1)}x{1:jmax}x{1:kmax} for F; {1:imax}x{1:(jmax-1)}x{1:kmax} for G, etc*/

	/*we only calculate U, V and W at the edges between two fluid cells*/
	for(i=1; i<imax; i++){
		for(j=1; j<jmax+1; j++){
      for(k=1; k<kmax+1; k++){
			  if (isfluid(i,j,k, Flag) && isfluid(i+1,j,k, Flag)){
				   U[i][j][k] = F[i][j][k] - dt/dx * (P[i+1][j][k] - P[i][j][k]);
			  }
      }
		}
	}

	for(i=1; i<imax+1; i++){
		for(j=1; j<jmax; j++){
      for(k=1; k<kmax+1; k++){
			  if (isfluid(i,j,k, Flag) && isfluid(i,j+1,k, Flag)){
				  V[i][j][k] = G[i][j][k] - dt/dy * (P[i][j+1][k] - P[i][j][k]);
			  }
      }
		}
	}

  for(i=1; i<imax+1; i++){
    for(j=1; j<jmax+1; j++){
      for(k=1; k<kmax; k++){
        if (isfluid(i,j,k, Flag) && isfluid(i,j,k+1, Flag)){
          V[i][j][k] = H[i][j][k] - dt/dz * (P[i][j][k+1] - P[i][j][k]);
        }
      }
    }
  }

}
