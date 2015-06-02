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
		    char *problem,
		    int **Flag,
		    double vel
                    ){

	switch(wl){ //left wall indices:u(0,i), (v(0,i)+v(1,i))/2
		case NO_SLIP :
			for(int j=1; j<=jmax; j++){
				U[0][j] = 0.0;
				V[0][j] = -V[1][j];
			}
			break;
		case FREE_SLIP :
			for(int j=1; j<=jmax; j++){
				U[0][j] = 0.0;
				V[0][j] = V[1][j];
			}
			break;
		case OUTFLOW :
			for(int j=1; j<=jmax; j++){
				U[0][j] = U[1][j];
				V[0][j] = V[1][j];
			}
			break;
	}

	switch(wr){ //right wall indices: u(imax,i), (v(imax,i)+v(imax+1,i))/2
		case NO_SLIP :
			for(int j=1; j<=jmax; j++){
				U[imax][j] = 0.0;
			        V[imax+1][j] = -V[imax][j];
			}
			break;
		case FREE_SLIP :
			for(int j=1; j<=jmax; j++){
				U[imax][j] = 0.0;
				V[imax+1][j] = V[imax][j];
			}
			break;
		case OUTFLOW :
			for(int j=1; j<=jmax; j++){
				U[imax][j] = U[imax-1][j];
				V[imax+1][j] = V[imax][j];
			}
			break;
	}

	switch(wt){ //top wall indices: v(i,jmax), (u(i,jmax)+u(i,jmax+1))/2
		case NO_SLIP :
			for(int i=1; i<=imax; i++){
				V[i][jmax] = 0.0;
			        U[i][jmax+1] = -U[i][jmax];
			}
			break;
		case FREE_SLIP :
			for(int i=1; i<=imax; i++){
				V[i][jmax] = 0.0;
				U[i][jmax+1] = U[i][jmax];
			}
			break;
		case OUTFLOW :
			for(int i=1; i<=imax; i++){
				V[i][jmax] = V[i][jmax-1];
				U[i][jmax+1] = U[i][jmax];
			}
			break;
	}

	switch(wb){ //bottom wall indices: v(i,0), (u(i,0)+u(i,1))/2
		case NO_SLIP :
			for(int i=1; i<=imax; i++){
				V[i][0] = 0.0;
			        U[i][0] = -U[i][1];
			}
			break;
		case FREE_SLIP :
			for(int i=1; i<=imax; i++){
				V[i][0] = 0.0;
				U[i][0] = U[i][1];
			}
			break;
		case OUTFLOW :
			for(int i=1; i<=imax; i++){
				V[i][0] = V[i][1];
				U[i][0] = U[i][1];
			}
			break;
	}

	//special boundaries
	spec_boundary_val(problem, imax, jmax, U, V, Flag, vel);

}

void spec_boundary_val(char *problem, int imax, int jmax, double **U, double **V, int **Flag, double vel){
	//take care of inflow velocity in different scenarios
	if(Flag[0][jmax/2]!=C_P){ //using 1 cell in left boundary to check if P given
		if (strcmp(problem,"KARMAN.pgm")!=0 || strcmp(problem, "SHEAR.pgm")!=0){
			for (int j=1; j<=jmax; j++)
				U[0][j] = vel;
	        } else if (strcmp(problem, "STEP.pgm")!=0){
			for (int j=0; j<=jmax/2; j++){//j=1
				//U[0][j] = vel;
				for (int i=0; i<=imax; i++){ //behind the step velocity is zero
					U[i][j] = 0.0;
				}
			}
			for (int j=jmax/2; j<=jmax; j++){
				U[0][j] = vel;
			}
		} else if (strcmp(problem, "DRIVEN_CAVITY.pgm")!=0){
			for (int i=0; i<=imax; i++){
				U[i][jmax+1] = vel*2.0 - U[i][jmax];
			}
		}
	}

	//take care of arbitrary boundaries
	for (int i=1; i<=imax; i++){
		for (int j=1; j<=jmax; j++){
		        switch(Flag[i][j]){ //should there be a check for a case of given P? (see Eva's ws2 notes, pg13)
                		case B_N:
					V[i][j] = 0;
					U[i-1][j] = -U[i-1][j+1];
					U[i][j] = -U[i][j+1];
					break;
                		case B_O:
                                        U[i][j] = 0;
                                        V[i][j-1] = -V[i+1][j-1];
                                        V[i][j] = -V[i+1][j];
					break;
	                	case B_S:
                                        V[i][j-1] = 0;
                                        U[i-1][j] = -U[i-1][j-1];
                                        U[i][j] = -U[i-1][j-1];
					break;
				case B_W:
                                        U[i-1][j] = 0;
                                        V[i][j-1] = -V[i-1][j];
                                        V[i][j] = -V[i-1][j];
					break;
                		case B_NO:
					V[i][j] = 0;
					U[i][j] = 0;
					U[i-1][j] = -U[i-1][j+1];
					V[i][j-1] = -V[i+1][j-1];
					break;
                		case B_NW:
                                        V[i][j] = 0;
                                        U[i-1][j] = 0;
                                        U[i][j] = -U[i][j+1];
                                        V[i][j-1] = -V[i-1][j-1];
					break;
		                case B_SO:
                                        V[i][j-1] = 0;
                                        U[i][j] = 0;
                                        U[i-1][j] = -U[i-1][j-1];
                                        V[i][j] = -V[i+1][j];
					break;
		                case B_SW:
                                        V[i][j-1] = 0;
                                        U[i-1][j] = 0;
                                        U[i][j] = -U[i][j-1];
                                        V[i][j] = -V[i-1][j];
					break;
				case C_B://added so the step wont look red.... is it okay?
					V[i][j] = 0;
					U[i][j] = 0;
					break;
			}
		}
	}
}
