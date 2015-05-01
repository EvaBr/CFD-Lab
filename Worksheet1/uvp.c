#include "helper.h"

void calculate_fg(
  double Re,
  double GX,
  double GY,
  double alpha,
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G
){
	double uij, vij, uuj, udj, uiu, uid, viu, vid, vuj, vdj, udu, vud;
	double Dx = 1/dx, Dy = 1/dy;
	int i, j;
	for (i=1; i<imax; i++){
		for  (j=1; j<jmax; j++){
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

			F[i][j] = uij + dt*(1/Re*((uuj-2*uij+udj)*Dx*Dx + (uiu-2*uij+uid)*Dy*Dy) -
				0.25*Dx*(pow((uij+uuj),2)-pow((udj+uij),2) + alpha*(abs(uij+uuj)*(uij-uuj)-abs(udj+uij)*(udj-uij))) -
				0.25*Dy*((vij+vuj)*(uij+uiu)-(vid+vud)*(uid+uij) + alpha*(abs(vij+vuj)*(uij-uiu)-abs(vid+vud)*(uid-uij))) +
				GX);
			G[i][j] = vij + dt*(1/Re*((vuj-2*vij+vdj)*Dx*Dx + (viu-2*vij+vid)*Dy*Dy) -
				0.25*Dy*(pow((vij+viu),2)-pow((vid+vij),2) + alpha*(abs(vij+viu)*(vij-viu)-abs(vid+vij)*(vid-vij))) -
				0.25*Dx*((uij+uiu)*(vij+vuj)-(udj+udu)*(vdj+vij) + alpha*(abs(uij+uiu)*(vij-vuj)-abs(udj+udu)*(vdj-vij))) +
				GY);

		}
		/* j=jmax, but we go through all i's */
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
		/* setting last F */
		F[i][j] = uij + dt*(1/Re*((uuj-2*uij+udj)*Dx*Dx + (uiu-2*uij+uid)*Dy*Dy) -
			0.25*Dx*(pow((uij+uuj),2)-pow((udj+uij),2) + alpha*(abs(uij+uuj)*(uij-uuj)-abs(udj+uij)*(udj-uij))) -
			0.25*Dy*((vij+vuj)*(uij+uiu)-(vid+vud)*(uid+uij) + alpha*(abs(vij+vuj)*(uij-uiu)-abs(vid+vud)*(uid-uij))) +
			GX);

		/* bound.cond. for G */
		G[i][0] = V[i][0];
		G[i][jmax] = V[i][jmax];
	}
	G[imax][0] = V[imax][0];
	G[imax][jmax] = V[imax][jmax];
	/* now i=imax, but we need to visit all j's */
	for (j=1; j<jmax+1; j++){
		/* setting the last G */
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

		G[i][j] = vij + dt*(1/Re*((vuj-2*vij+vdj)*Dx*Dx + (viu-2*vij+vid)*Dy*Dy) -
			0.25*Dy*(pow((vij+viu),2)-pow((vid+vij),2) + alpha*(abs(vij+viu)*(vij-viu)-abs(vid+vij)*(vid-vij))) -
			0.25*Dx*((uij+uiu)*(vij+vuj)-(udj+udu)*(vdj+vij) + alpha*(abs(uij+uiu)*(vij-vuj)-abs(udj+udu)*(vdj-vij))) +
			GY);

		/* bound.cond. for F */
		F[0][j] = U[0][j];
		F[imax][j] = U[imax][j];
	}
};

void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS
);

void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V
){
	/* we rewrite dt only if tau is nonnegative, else we do nothing */
	if(tau>=0){
		*dt = tau*fmin( fmin( pow(0.5*Re*(pow(dx,-2) + pow(dy,-2)),-1),  dx/mmax(**U, imax, jmax)),  dy/mmax(**V, imax, jmax));
 	}
};

void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P
);

#endif
