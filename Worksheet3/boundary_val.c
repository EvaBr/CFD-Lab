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
		    char *problem
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


	        if (strcmp(problem,"KARMAN")!=0){
                	        //TODO
        	}
	        else if (strcmp(problem,"SHEAR")!=0){
                	        //TODO
        	}
	        else if (strcmp(problem,"STEP")!=0){

			//ignore interior obstacles
			for(int i=1; i<jmax/2; i++){ //B_N
				V[i][jmax/2] = 0;
				U[i-1][jmax/2] = -U[i-1][jmax/2+1];
				U[i][jmax/2] = -U[i][jmax/2+1];
				G[i][jmax/2] = V[i][jmax/2];
				P[i][jmax/2] = P[i][jmax/2+1];
			}
			for(int j=1; j<jmax/2; j++){ //B_O
				U[jmax/2][j] = 0;
				V[jmax/2][j] = -V[jmax/2+1][j];
				V[jmax/2][j-1] = -V[jmax/2+1][j-1];
				F[jmax/2][j] = U[jmax/2][j];
				P[jmax/2][j] = P[jmax/2+1][j];
			}
			//B_NO
			U[jmax/2][jmax/2] = 0;
			V[jmax/2][jmax/2] = 0;
			U[jmax/2-1][jmax/2] = -U[jmax/2-1][jmax/2+1];
			V[jmax/2][jmax/2-1] = -V[jmax/2+1][jmax/2-1];
			F[jmax/2][jmax/2] = U[jmax/2][jmax/2];
			G[jmax/2][jmax/2] = V[jmax/2][jmax/2];
			P[jmax/2][jmax/2] = (P[jmax/2][jmax/2+1] + P[jmax/2+1][jmax])/2;
		}
}

void spec_boundary_var(char *problem, int imax, int jmax, double **U, double **V){
	
}
