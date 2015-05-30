#include "init.h"
#include "boundary_val.h"

/**
 * The boundary values of the 4 domain walls are set.
 */
void boundaryvalues(
                    int imax,
                    int jmax,
                    double **U,
                    double **V,
		    int wl,
		    int wr,
		    int wt,
		    int wb,
                    ){

	switch(wl){ //left wall indeces:u(0,i), (v(0,i)+v(1,i))/2
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

	switch(wr){ //right wall indeces: u(imax,i), (v(imax,i)+v(imax+1,i))/2
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

	switch(wt){ //top wall indeces: v(i,jmax), (u(i,jmax)+u(i,jmax+1))/2
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

	switch(wb){ //bottom wall indeces: v(i,0), (u(i,0)+u(i,1))/2
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

}

void spec_boundary_var(char *problem, int imax, int jmax, double **U, double **V){
	
}
