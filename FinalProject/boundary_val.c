#include "init.h"
#include "boundary_val.h"
#include "helper.h"

/**
 * No slip boundary values of the 6 domain walls are set.
 */

void boundaryvalues_no_slip(
                    int imax,
                    int jmax,
                    int kmax,
                    double ***U,
                    double ***V,
                    double ***W,
                    double ***P,
                    double ***F,
                    double ***G,
                    double ***H,
                    int ***Flag
                    ){
        // No-slip boundary conditions for U, V and W.
        for(int j=1; j<=jmax; j++){
                for(int k=1; k<=kmax; k++){
                        U[0][j][k] = 0.0;
                        U[imax][j][k] = 0.0;
                        V[0][j][k] = -V[1][j][k];
                        V[imax+1][j][k] = -V[imax][j][k];
                        W[0][j][k] = -W[1][j][k];
                        W[imax+1][j][k] = -W[imax][j][k];
                }
        }

        for(int i=1; i<=imax; i++){
                for(int k=1; k<=kmax; k++){
                        V[i][0][k] = 0.0;
                        V[i][jmax][k] = 0.0;
                        U[i][0][k] = -U[i][1][k];
                        U[i][jmax+1][k] = -U[i][jmax][k];
                        W[i][0][k] = -W[i][1][k];
                        W[i][jmax+1][k] = -W[i][jmax][k];
                }
        }

        for(int i=1; i<=imax; i++){
                for(int j=1; j<=jmax; j++){
                        W[i][j][0] = 0.0;
                        W[i][j][kmax] = 0.0;
                        U[i][j][0] = -U[i][j][1];
                        U[i][j][kmax+1] = -U[i][j][kmax];
                        V[i][j][0] = -V[i][j][1];
                        V[i][j][kmax+1] = -V[i][j][kmax];
                }
        }

}


/**
 * The boundary values of the 4 domain walls are set.
 */

void boundaryvalues(
        int imax,
        int jmax,
        int kmax,
        double ***U,
        double ***V,
        double ***W,
        double ***P,
		    int wl,
		    int wr,
        int wf,
        int wh,
		    int wt,
		    int wb,
		    double ***F,
		    double ***G,
        double ***H,
		    char *problem,
		    int ***Flag,
		    double velIN,
        double velMW
                    ){

	switch(wl){ //left wall indices:u(0,i), (v(0,i)+v(1,i))/2
		case NO_SLIP :
      for(int j=1; j<=jmax; j++){
        for(int k=1; k<=kmax; k++){
				    U[0][j][k] = 0.0;
				    V[0][j][k] = -V[1][j][k];
            W[0][j][k] = -W[1][j][k];
			}
			break;
		case FREE_SLIP :
			for(int j=1; j<=jmax; j++){
        for(int k=1; k<=kmax; k++){
				  U[0][j][k] = 0.0;
          V[0][j][k] = V[1][j][k];
          W[0][j][k] = W[1][j][k];
			}
			break;
    case OUTFLOW :
			for(int j=1; j<=jmax; j++){
        for(int k=1; k<=kmax; k++){
	          U[0][j][k] = U[1][j][k];
  			    V[0][j][k] = V[1][j][k];
  		      W[0][j][k] = W[1][j][k];
      }
  		break;
    case INFLOW : //currently, only supported inflow is one perpendicular to the wall of the inflow source
  		for(int j=1; j<=jmax; j++){
        for(int k=1; k<=kmax; k++){
  	        U[0][j][k] = velIN;
  			    V[0][j][k] = - V[1][j][k];
  		      W[0][j][k] = - W[1][j][k];
      }
  		break;
    case MOWING_WALL : //currently only supported moving wall direction is to the + direction of next (first nonfixed) coordinate (e.g. if x/y/z is fixed, were moving in +y/+z/+x)
  		for(int j=1; j<=jmax; j++){
        for(int k=1; k<=kmax; k++){
  	        U[0][j][k] = 0.0;
  			    V[0][j][k] = 2.0*velMW - V[1][j][k];
  		      W[0][j][k] = - W[1][j][k];
      }
  		break;
	}

	switch(wr){ //right wall indices: u(imax,i), (v(imax,i)+v(imax+1,i))/2
  case NO_SLIP :
    for(int j=1; j<=jmax; j++){
      for(int k=1; k<=kmax; k++){
          U[imax][j][k] = 0.0;
          V[imax+1][j][k] = -V[imax][j][k];
          W[imax+1][j][k] = -W[imax][j][k];
    }
    break;
  case FREE_SLIP :
    for(int j=1; j<=jmax; j++){
      for(int k=1; k<=kmax; k++){
        U[imax][j][k] = 0.0;
        V[imax+1][j][k] = V[imax][j][k];
        W[imax+1][j][k] = W[imax][j][k];
    }
    break;
  case OUTFLOW :
    for(int j=1; j<=jmax; j++){
      for(int k=1; k<=kmax; k++){
          U[imax][j][k] = U[imax-1][j][k];
          V[imax+1][j][k] = V[imax][j][k];
          W[imax+1][j][k] = W[imax][j][k];
    }
    break;
  case INFLOW : //currently, only supported inflow is one perpendicular to the wall of the inflow source
    for(int j=1; j<=jmax; j++){
      for(int k=1; k<=kmax; k++){
          U[imax][j][k] = -velIN; //TODO: is this right?
          V[imax+1][j][k] = - V[imax][j][k];
          W[imax+1][j][k] = - W[imax][j][k];
    }
    break;
  case MOWING_WALL : //currently only supported moving wall direction is to the + direction of next (first nonfixed) coordinate (e.g. if x/y/z is fixed, were moving in +y/+z/+x)
    for(int j=1; j<=jmax; j++){
      for(int k=1; k<=kmax; k++){
          U[imax][j][k] = 0.0;
          V[imax+1][j][k] = 2.0*velMW - V[imax][j][k];
          W[imax+1][j][k] = - W[imax][j][k];
    }
    break;
	}

	switch(wt){ //top wall indices: v(i,jmax), (u(i,jmax)+u(i,jmax+1))/2
  case NO_SLIP :
    for(int i=1; i<=imax; i++){
      for(int j=1; j<=jmax; j++){
          U[i][j][kmax+1] = -U[i][j][kmax];
          V[i][j][kmax+1] = -V[i][j][kmax];
          W[i][j][kmax] = 0.0;
    }
    break;
  case FREE_SLIP :
    for(int i=1; i<=imax; i++){
      for(int j=1; j<=jmax; j++){
        U[i][j][kmax+1] = U[i][j][kmax];
        V[i][j][kmax+1] = V[i][j][kmax];
        W[i][j][kmax] = 0.0;
    }
    break;
  case OUTFLOW :
    for(int i=1; i<=imax; i++){
      for(int j=1; j<=jmax; j++){
          U[i][j][kmax+1] = U[i][j][kmax];
          V[i][j][kmax+1] = V[i][j][kmax];
          W[i][j][kmax] = W[i][j][kmax-1];
    }
    break;
  case INFLOW : //currently, only supported inflow is one perpendicular to the wall of the inflow source
    for(int i=1; i<=imax; i++){
      for(int j=1; j<=jmax; j++){
          U[i][j][kmax+1] = - U[i][j][kmax];
          V[i][j][kmax+1] = - V[i][j][kmax];
          W[i][j][kmax] = - velIN;
    }
    break;
  case MOWING_WALL : //currently only supported moving wall direction is to the + direction of next (first nonfixed) coordinate (e.g. if x/y/z is fixed, were moving in +y/+z/+x)
    for(int i=1; i<=imax; i++){
      for(int j=1; j<=jmax; j++){
          U[i][j][kmax+1] = 2.0*velMW - U[i][j][kmax];
          V[i][j][kmax+1] = - V[i][j][kmax];
          W[i][j][kmax] = 0.0;
    }
    break;
  }

	switch(wb){ //bottom wall indices: v(i,0), (u(i,0)+u(i,1))/2
  case NO_SLIP :
    for(int i=1; i<=imax; i++){
      for(int j=1; j<=jmax; j++){
          U[i][j][0] = -U[i][j][1];
          V[i][j][0] = -V[i][j][1];
          W[i][j][0] = 0.0;
    }
    break;
  case FREE_SLIP :
    for(int i=1; i<=imax; i++){
      for(int j=1; j<=jmax; j++){
        U[i][j][0] = U[i][j][1];
        V[i][j][0] = V[i][j][1];
        W[i][j][0] = 0.0;
    }
    break;
  case OUTFLOW :
    for(int i=1; i<=imax; i++){
      for(int j=1; j<=jmax; j++){
          U[i][j][0] = U[i][j][1];
          V[i][j][0] = V[i][j][1];
          W[i][j][0] = W[i][j][1];
    }
    break;
  case INFLOW : //currently, only supported inflow is one perpendicular to the wall of the inflow source
    for(int i=1; i<=imax; i++){
      for(int j=1; j<=jmax; j++){
          U[i][j][0] = - U[i][j][1];
          V[i][j][0] = - V[i][j][1];
          W[i][j][0] = velIN;
    }
    break;
  case MOWING_WALL : //currently only supported moving wall direction is to the + direction of next (first nonfixed) coordinate (e.g. if x/y/z is fixed, were moving in +y/+z/+x)
    for(int i=1; i<=imax; i++){
      for(int j=1; j<=jmax; j++){
          U[i][j][0] = 2.0*velMW - U[i][j][1];
          V[i][j][0] = - V[i][j][1];
          W[i][j][0] = 0.0;
    }
    break;
	}


  switch(wf){
  case NO_SLIP :
    for(int k=1; k<=kmax; k++){
      for(int i=1; i<=imax; i++){
          U[i][0][k] = - U[i][1][k];
          V[i][0][k] = 0.0;
          W[i][0][k] = - W[i][1][k];
    }
    break;
  case FREE_SLIP :
    for(int k=1; k<=kmax; k++){
      for(int i=1; i<=imax; i++){
        U[i][0][k] = U[i][1][k];
        V[i][0][k] = 0.0;
        W[i][0][k] = W[i][1][k];
    }
    break;
  case OUTFLOW :
  for(int k=1; k<=kmax; k++){
    for(int i=1; i<=imax; i++){
      U[i][0][k] = U[i][1][k];
      V[i][0][k] = V[i][1][k];
      W[i][0][k] = W[i][1][k];
    }
    break;
  case INFLOW : //currently, only supported inflow is one perpendicular to the wall of the inflow source
  for(int k=1; k<=kmax; k++){
    for(int i=1; i<=imax; i++){
      U[i][0][k] = - U[i][1][k];
      V[i][0][k] = velIN;
      W[i][0][k] = - W[i][1][k];
    }
    break;
  case MOWING_WALL : //currently only supported moving wall direction is to the + direction of next (first nonfixed) coordinate (e.g. if x/y/z is fixed, were moving in +y/+z/+x)
  for(int k=1; k<=kmax; k++){
    for(int i=1; i<=imax; i++){
      U[i][0][k] = - U[i][1][k];
      V[i][0][k] = 0.0;
      W[i][0][k] = 2*velMW - W[i][1][k];
    }
    break;
  }

  switch(wh){
  case NO_SLIP :
    for(int k=1; k<=kmax; k++){
      for(int i=1; i<=imax; i++){
          U[i][jmax+1][k] = - U[i][jmax][k];
          V[i][jmax][k] = 0.0;
          W[i][jmax+1][k] = - W[i][jmax][k];
    }
    break;
  case FREE_SLIP :
    for(int k=1; k<=kmax; k++){
      for(int i=1; i<=imax; i++){
        U[i][jmax+1][k] = U[i][jmax][k];
        V[i][jmax][k] = 0.0;
        W[i][jmax+1][k] = W[i][jmax][k];
    }
    break;
  case OUTFLOW :
  for(int k=1; k<=kmax; k++){
    for(int i=1; i<=imax; i++){
      U[i][jmax+1][k] = U[i][jmax][k];
      V[i][jmax][k] = V[i][jmax-1][k];
      W[i][jmax+1][k] = W[i][jmax][k];
    }
    break;
  case INFLOW : //currently, only supported inflow is one perpendicular to the wall of the inflow source
  for(int k=1; k<=kmax; k++){
    for(int i=1; i<=imax; i++){
      U[i][jmax+1][k] = - U[i][jmax][k];
      V[i][jmax][k] = - velIN;
      W[i][jmax+1][k] = - W[i][jmax][k];
    }
    break;
  case MOWING_WALL : //currently only supported moving wall direction is to the + direction of next (first nonfixed) coordinate (e.g. if x/y/z is fixed, were moving in +y/+z/+x)
  for(int k=1; k<=kmax; k++){
    for(int i=1; i<=imax; i++){
      U[i][jmax+1][k] = - U[i][jmax][k];
      V[i][jmax][k] = 0.0;
      W[i][jmax+1][k] = 2*velMW - W[i][jmax][k];
    }
    break;
  }
	//special boundaries
	//spec_boundary_val(problem, imax, jmax, U, V, Flag, vel);

}

/*
void spec_boundary_val(char *problem, int imax, int jmax, double **U, double **V, int **Flag, double vel){
	//take care of inflow velocity in different scenarios
	if((Flag[0][jmax/2] & 32) == 0){ //using 1 cell in left boundary to check if P given
		if (strcmp(problem,"KARMAN.pgm")!=0 || strcmp(problem, "SHEAR.pgm")!=0){
			for (int j=1; j<=jmax; j++){
				U[0][j] = vel;
				V[0][j] = -V[1][j]; //V is set to 0 on the boundary
			}
	        } else if (strcmp(problem, "STEP.pgm")!=0){
			for (int j=jmax/2; j<=jmax; j++){
				U[0][j] = vel;
				V[0][j] = -V[1][j];
			}
		} else if (strcmp(problem, "DRIVEN_CAVITY.pgm")!=0){
			for (int i=1; i<=imax; i++){
				U[i][jmax+1] = vel*2.0 - U[i][jmax];
			}
		}
	}
*/
	//take care of arbitrary boundaries
	for (int i=1; i<=imax; i++){
		for (int j=1; j<=jmax; j++){
		        switch(Flag[i][j]){
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
                                        U[i][j] = -U[i][j-1];
					break;
				case B_W:
                                        U[i-1][j] = 0;
                                        V[i][j-1] = -V[i-1][j-1];
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
				case C_B://added so the insides of obstacles wont be red
					V[i][j] = 0;
					U[i][j] = 0;
					break;
			}
		}
	}
	//added another forloops for setting outside boundary C_B cells to U=0 and V=0
	for (int i=0; i<=imax+1; i++){
		if((Flag[i][0] & 31) == C_B){ //(Flag[i][0] & 31) gets rid of the C_P bit, for left boundary
			V[i][0] = 0;
			U[i][0] = 0;
		}
		if((Flag[i][jmax+1] & 31) == C_B){ //(Flag[i][0] & 31) gets rid of the C_P bit, for right boundary
			V[i][jmax+1] = 0;
			U[i][jmax+1] = 0;
		}
	}
	for (int j=0; j<=jmax+1; j++){
		if(Flag[0][j] == C_B){
			V[0][j] = 0;
			U[0][j] = 0;
		}
		if(Flag[imax+1][j] == C_B){
			V[imax+1][j] = 0;
			U[imax+1][j] = 0;
		}

	}
}
