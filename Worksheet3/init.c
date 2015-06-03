#include "helper.h"
#include "init.h"
#include <math.h>

int read_parameters( const char *szFileName,       /* name of the file */
                    double *Re,                /* reynolds number   */
                    double *UI,                /* velocity x-direction */
                    double *VI,                /* velocity y-direction */
                    double *PI,                /* pressure */
                    double *GX,                /* gravitation x-direction */
                    double *GY,                /* gravitation y-direction */
                    double *t_end,             /* end time */
                    double *xlength,           /* length of the domain x-dir.*/
                    double *ylength,           /* length of the domain y-dir.*/
                    double *dt,                /* time step */
                    double *dx,                /* length of a cell x-dir. */
                    double *dy,                /* length of a cell y-dir. */
                    int  *imax,                /* number of cells x-direction*/
                    int  *jmax,                /* number of cells y-direction*/
                    double *alpha,             /* uppwind differencing factor*/
                    double *omg,               /* relaxation factor */
                    double *tau,               /* safety factor for time step*/
                    int  *itermax,             /* max. number of iterations  */
		                               /* for pressure per time step */
                    double *eps,               /* accuracy bound for pressure*/
		    double *dt_value,          /* time for output */

		    int *wl,			/*initial boundary for left wall*/
		    int *wr,			/*initial boundary for right wall*/
		    int *wt,			/*initial boundary for top wall*/
		    int *wb,			/*initial boundary for bottom wall*/
		    char *problem,		/*problem to solve*/

		    double *presLeft,		/*pressure at the left wall*/
		    double *presRight,		/*pressure at the right wall*/
		    double *presDelta,		/*pressure difference across the domain*/

		    double *vel)		/*velocity of inflow or wall (in U direction)*/
{
   READ_DOUBLE( szFileName, *xlength );
   READ_DOUBLE( szFileName, *ylength );

   READ_DOUBLE( szFileName, *Re    );
   READ_DOUBLE( szFileName, *t_end );
   READ_DOUBLE( szFileName, *dt    );

   READ_INT   ( szFileName, *imax );
   READ_INT   ( szFileName, *jmax );

   READ_DOUBLE( szFileName, *omg   );
   READ_DOUBLE( szFileName, *eps   );
   READ_DOUBLE( szFileName, *tau   );
   READ_DOUBLE( szFileName, *alpha );

   READ_INT   ( szFileName, *itermax );
   READ_DOUBLE( szFileName, *dt_value );

   READ_DOUBLE( szFileName, *UI );
   READ_DOUBLE( szFileName, *VI );
   READ_DOUBLE( szFileName, *GX );
   READ_DOUBLE( szFileName, *GY );
   READ_DOUBLE( szFileName, *PI );

   READ_INT   ( szFileName, *wl );
   READ_INT   ( szFileName, *wr );
   READ_INT   ( szFileName, *wt );
   READ_INT   ( szFileName, *wb );

   READ_STRING( szFileName, problem );

   READ_DOUBLE( szFileName, *presLeft);
   READ_DOUBLE( szFileName, *presRight);
   READ_DOUBLE( szFileName, *presDelta);

   READ_DOUBLE( szFileName, *vel   );

   //take care of (in)valid pressure input
   if (*presDelta<=0){
   //    if (fmin(presLeft, presRight)<0): we dont have pressure input
	if  (fmin(*presLeft, *presRight)>=0){
		*presDelta = *presLeft - *presRight;
        }
   } else { //deltaP is given
	if  (*presLeft< *presDelta){
		if (*presRight<0){
			*presLeft = *presDelta;
			*presRight = 0.0;
		} else {
			*presLeft = *presRight + *presDelta;
		}
	} else {//pressure on left wall is also given
		*presRight = *presLeft - *presDelta;
	}
   }
   if (*presDelta>0){ // if pressure given, left and right bound. set to outflow
	*wl = 3;
	*wr = 3;
   }

   *dx = *xlength / (double)(*imax);
   *dy = *ylength / (double)(*jmax);

   return 1;
}

/**
 * The arrays U,V and P are initialized to the constant values UI, VI and PI on
 * the whole domain.
 */
void init_uvp(
  double UI,
  double VI,
  double PI,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **P
) {
	init_matrix(U, 0, imax+1, 0, jmax+1, UI);
	init_matrix(V, 0, imax+1, 0, jmax+1, VI);
	init_matrix(P, 0, imax+1, 0, jmax+1, PI);
}

/**
 * The integer array Flag is initialized to constants C_F for fluid cells and C_B
 * for obstacle cells as specified by the parameter problem.
 */
void init_flag(
  char *problem,
  int imax,
  int jmax,
  double presDelta,
  int **Flag
) {
	int i,j;
	int temp;
	//read the geometry
	int **Pic = read_pgm(problem);

	//initialisation to C_F and C_B
	for (int i=1; i<imax+1; i++){
		for (int j=1; j<jmax+1; j++){
			temp = min(Pic[i][j]*pow(2,4) + Pic[i+1][j]*pow(2,3) + Pic[i-1][j]*pow(2,2) + Pic[i][j-1]*2 + Pic[i][j+1], 16);
			if (temp == 3 || temp ==7 || (temp > 10 && temp < 15)){
				ERROR("Invalid geometry! Forbidden boundary cell found.\n"); }
			Flag[i][j] = temp;
		}
	}


	//set outer boundary flags - upper and lower:
	for (i=0; i<=imax+1; i++){
		if (Flag[i][1] == C_F) {
			Flag[i][0] = B_N;
		} else {
			Flag[i][0] = C_B;
		}
		if (Flag[i][jmax] == C_F) {
			Flag[i][jmax+1] = B_S;
		} else {
			Flag[i][jmax+1] = C_B;
		}
	}
	//set the outer boundary flags - right and left:
	for (j=0; j<=jmax+1; j++){
		if (Flag[1][j]==C_F) {
			Flag[0][j] = B_O;
		} else {
			Flag[0][j] = C_B;
		}
		if (Flag[imax][j] == C_F) {
			Flag[imax+1][j] = B_W;
		} else {
			Flag[imax+1][j] = C_B;
		}
	}


	// take care of the case when pressure is given
	for (i=0; i<=imax+1; i++){
		if (presDelta) {
			Flag[i][0] += 32;
			Flag[i][jmax+1] += 32;
		}
	}
	for (j=0; j<=jmax+1; j++){
		if (presDelta) {
			Flag[0][j] += 32;
			Flag[imax+1][j] += 32;
		}
	}

}
