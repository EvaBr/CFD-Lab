#include "helper.h"
#include "init.h"

#define C_F 100 //C_F value for fluid cells are temporarily set as 100
#define C_B 0 //interior of the obstacle. flag 00000
#define B_N 1 //northern edge cell. flag 00001
#define B_O 8 //eastern edge cell. flag 01000
#define B_NO 9 //northeastern edge cell. flag 01001
/*more Flag values to be defined*/


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
		    int *wb)			/*initial boundary for bottom wall*/

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
  char *problem
  int imax,
  int jmax,
  int **Flag
) {
	int i,j;
	//initialisation to C_F and C_B
	init_matrix(Flag, 0, imax+1, 0, jmax+1, C_F); //C_F value for fluid cells are temporarily set as 100

	for (i=0; i<=imax+1; i++){
		Flag[i][0] = C_B;
		Flag[i][jmax+1] = C_B;
	}
	for (j=0; j<=jmax+1; j++){
		Flag[0][j] = C_B;
		Flag[imax+1][j] = C_B;
	}
	//looping over all cells where the boundary cells are marked with 
	//the appropriate flags B_xy
        switch(problem){ //
                case KARMAN :
			//TODO
                case PLANE_SHEAR :
			//TODO
                case STEP :
			for(i=1; i<jmax/2; i++){
				for(j=1; j<jmax/2; j++){
					Flag[i][j] = C_B; //interior of the obstacle. flag 00000
				}
			}
			for(i=1; i<jmax/2; i++){
				Flag[i][jmax/2] = B_N; //northern edge cell. flag 00001
			}
			for(j=1; j<jmax/2; j++){
				Flag[jmax/2][j] = B_O; //eastern edge cell. flag 01000
			}
			Flag[jmax/2][jmax/2] = B_NO; //northeastern edge cell. flag 01001
        }

}
