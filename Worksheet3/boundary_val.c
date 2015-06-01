#include "init.h"
#include "boundary_val.h"
#include "helper.h"

/**
 * The boundary values of the 4 domain walls are set.
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
		    double **Flag,
		    char *problem,
		    double UI
                    ){

	switch(wl){ //left wall indices:u(0,i), (v(0,i)+v(1,i))/2
		case NO_SLIP :
			for(int j=1; j<=jmax; j++){
				U[0][j] = 0.0;
				V[0][j] = -V[1][j];
			}
		case FREE_SLIP :
			for(int j=1; j<=jmax; j++){
				U[0][j] = 0.0;
				V[0][j] = V[1][j];
			}
		case OUTFLOW :
			for(int j=1; j<=jmax; j++){
				U[0][j] = U[1][j];
				V[0][j] = V[1][j];
			}
	}

	switch(wr){ //right wall indices: u(imax,i), (v(imax,i)+v(imax+1,i))/2
		case NO_SLIP :
			for(int j=1; j<=jmax; j++){
				U[imax][j] = 0.0;
			        V[imax+1][j] = -V[imax][j];
			}
		case FREE_SLIP :
			for(int j=1; j<=jmax; j++){
				U[imax][j] = 0.0;
				V[imax+1][j] = V[imax][j];
			}
		case OUTFLOW :
			for(int j=1; j<=jmax; j++){
				U[imax][j] = U[imax-1][j];
				V[imax+1][j] = V[imax][j];
			}
	}

	switch(wt){ //top wall indices: v(i,jmax), (u(i,jmax)+u(i,jmax+1))/2
		case NO_SLIP :
			for(int i=1; i<=imax; i++){
				V[i][jmax] = 0.0;
			        U[i][jmax+1] = -U[i][jmax];
			}
		case FREE_SLIP :
			for(int i=1; i<=imax; i++){
				V[i][jmax] = 0.0;
				U[i][jmax+1] = U[i][jmax];
			}
		case OUTFLOW :
			for(int i=1; i<=imax; i++){
				V[i][jmax] = V[i][jmax-1];
				U[i][jmax+1] = U[i][jmax];
			}
	}

	switch(wb){ //bottom wall indices: v(i,0), (u(i,0)+u(i,1))/2
		case NO_SLIP :
			for(int i=1; i<=imax; i++){
				V[i][0] = 0.0;
			        U[i][0] = -U[i][1];
			}
		case FREE_SLIP :
			for(int i=1; i<=imax; i++){
				V[i][0] = 0.0;
				U[i][0] = U[i][1];
			}
		case OUTFLOW :
			for(int i=1; i<=imax; i++){
				V[i][0] = V[i][1];
				U[i][0] = U[i][1];
			}
	}

	//Boundary Values for F, G, P - discrete Neumann cond.
/*	for (int j=1; j<=jmax; j++){
		P[0][j] = P[1][j];
		P[imax+1][j] = P[imax][j];

		F[0][j] = U[0][j];
		F[imax][j] = U[imax][j];
	}
	for (int i=1; i<=imax; i++){
		P[i][0] = P[i][1];
		P[i][jmax+1] = P[i][jmax];

		G[i][0] = V[i][0];
		G[i][jmax] = V[i][jmax];
	}*/


	//special boundaries
	spec_boundary_val(problem, imax, jmax, U, V, Flag, UI);

}

void spec_boundary_val(char *problem, int imax, int jmax, double **U, double **V, double **Flag, double UI){
	if (strcmp(problem,"KARMAN.pgm")!=0 || strcmp(problem, "SHEAR.pgm")!=0){
		for (int j=1; j<=jmax; j++)
			U[0][j] = UI;
        } else if (strcmp(problem, "STEP.pgm")!=0){
		for (int j=1; j<=jmax/2; j++){
			U[0][j] = UI;
			for (int i=1; i<=imax; i++){
				U[i][j] = 0.0;
			}
		}
		for (int j=jmax/2; j<=jmax; j++){
			U[0][j] = UI;

	}
	//TODO case of pressure??

	//Take care of arbitrary geometries
	//8 cases
	for (int i=1; i<=imax; i++){
		for (int j=1; j<=jmax; j++){
		        switch(Flag[i][j]){
                		case B_N:
					V[i][j] = 0;
					U[i-1][j] = -U[i-1][j+1];
					U[i][j] = -U[i][j+1];
					G[i][j] = V[i][j];
					P[i][j] = P[i][j+1];
                		case B_O:
                                        U[i][j] = 0;
                                        V[i][j-1] = -V[i+1][j-1];
                                        V[i][j] = -V[i+1][j];
                                        F[i][j] = U[i][j];
                                        P[i][j] = P[i+1][j];
	                	case B_S:
                                        V[i][j-1] = 0;
                                        U[i-1][j] = -U[i-1][j-1];
                                        U[i][j] = -U[i-1][j-1];
                                        G[i][j-1] = V[i][j-1];
                                        P[i][j] = P[i][j-1];
				case B_W:
                                        U[i-1][j] = 0;
                                        V[i][j-1] = -V[i-1][j];
                                        V[i][j] = -V[i-1][j];
                                        F[i-1][j] = U[i-1][j];
                                        P[i][j] = P[i-1][j];

                		case B_NO:
					V[i][j] = 0;
					U[i][j] = 0;
					U[i-1][j] = -U[i-1][j+1];
					V[i][j-1] = -V[i+1][j-1];
					F[i][j] = U[i][j];
					G[i][j] = V[i][j];
					P[i][j] = (P[i+1][j] + P[i][j+1])*0.5;
                		case B_NW:
                                        V[i][j] = 0;
                                        U[i-1][j] = 0;
                                        U[i][j] = -U[i][j+1];
                                        V[i][j-1] = -V[i-1][j-1];
                                        F[i-1][j] = U[i-1][j];
                                        G[i][j] = V[i][j];
                                        P[i][j] = (P[i-1][j] + P[i][j+1])*0.5;
		                case B_SO:
                                        V[i][j-1] = 0;
                                        U[i][j] = 0;
                                        U[i-1][j] = -U[i-1][j-1];
                                        V[i][j] = -V[i+1][j];
                                        F[i][j] = U[i][j];
                                        G[i][j-1] = V[i][j-1];
                                        P[i][j] = (P[i+1][j] + P[i][j-1])*0.5;
		                case B_SW:
                                        V[i][j-1] = 0;
                                        U[i-1][j] = 0;
                                        U[i][j] = -U[i][j-1];
                                        V[i][j] = -V[i-1][j];
                                        F[i-1][j] = U[i-1][j];
                                        G[i][j-1] = V[i][j-1];
                                        P[i][j] = (P[i-1][j] + P[i][j-1])*0.5;
			}
		}
	}
}
