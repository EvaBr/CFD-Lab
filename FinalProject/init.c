#include "helper.h"
#include "init.h"
#include <math.h>

int read_parameters( const char *szFileName,       /* name of the file */
                    double *Re,                /* reynolds number   */
                    double *UI,                /* velocity x-direction */
                    double *VI,                /* velocity y-direction */
                    double *WI,                /* velocity z-direction */
                    double *PI,                /* pressure */
                    double *GX,                /* gravitation x-direction */
                    double *GY,                /* gravitation y-direction */
                    double *GZ,                /* gravitation z-direction */
                    double *t_end,             /* end time */
                    double *xlength,           /* length of the domain x-dir.*/
                    double *ylength,           /* length of the domain y-dir.*/
                    double *zlength,           /* length of the domain z-dir.*/
                    double *dt,                /* time step */
                    double *dx,                /* length of a cell x-dir. */
                    double *dy,                /* length of a cell y-dir. */
                    double *dz,                /* length of a cell z-dir. */
                    int  *imax,                /* number of cells x-direction*/
                    int  *jmax,                /* number of cells y-direction*/
                    int  *kmax,                /* number of cells z-direction*/
                    double *alpha,             /* uppwind differencing factor*/
                    double *omg,               /* relaxation factor */
                    double *tau,               /* safety factor for time step*/
                    int  *itermax,             /* max. number of iterations for pressure per time step */
                    double *eps,               /* accuracy bound for pressure*/
		                double *dt_value,          /* time for output */
                    int *wl,              		 /*initial boundary for left wall*/
                    int *wr,              		 /*initial boundary for right wall*/
                    int *wf,		               /*initial boundary for front wall*/
                    int *wh,			             /*initial boundary for back wall*/
                    int *wt,			             /*initial boundary for top wall*/
                    int *wb,			             /*initial boundary for bottom wall*/
                    char *problem,		         /*problem to solve*/
//                  double *presLeft,		       /*pressure at the left wall*/
//                  double *presRight,		     /*pressure at the right wall*/
//                  double *presDelta,		     /*pressure difference across the domain*/
		                double *velIN,             /*velocity of inflow*/
                    double *velMW )		         /*velocity of wall (in U direction)*/
{
   READ_DOUBLE( szFileName, *xlength );
   READ_DOUBLE( szFileName, *ylength );
   READ_DOUBLE( szFileName, *zlength );

   READ_INT   ( szFileName, *imax );
   READ_INT   ( szFileName, *jmax );
   READ_INT   ( szFileName, *kmax );

   READ_DOUBLE( szFileName, *dt    );
   READ_DOUBLE( szFileName, *t_end );
   READ_DOUBLE( szFileName, *tau   );

   READ_DOUBLE( szFileName, *dt_value );
   READ_INT   ( szFileName, *itermax );

   READ_DOUBLE( szFileName, *eps   );
   READ_DOUBLE( szFileName, *omg   );
   READ_DOUBLE( szFileName, *alpha );

   READ_DOUBLE( szFileName, *Re    );

   READ_DOUBLE( szFileName, *GX );
   READ_DOUBLE( szFileName, *GY );
   READ_DOUBLE( szFileName, *GZ );

   READ_DOUBLE( szFileName, *PI );
   READ_DOUBLE( szFileName, *UI );
   READ_DOUBLE( szFileName, *VI );
   READ_DOUBLE( szFileName, *WI );

   READ_INT   ( szFileName, *wl );
   READ_INT   ( szFileName, *wr );
   READ_INT   ( szFileName, *wf );
   READ_INT   ( szFileName, *wh );
   READ_INT   ( szFileName, *wt );
   READ_INT   ( szFileName, *wb );

   READ_STRING( szFileName, problem );

/*   READ_DOUBLE( szFileName, *presLeft);
   READ_DOUBLE( szFileName, *presRight);
   READ_DOUBLE( szFileName, *presDelta); */

   READ_DOUBLE( szFileName, *velIN );
   READ_DOUBLE( szFileName, velMW[0] );
   READ_DOUBLE( szFileName, velMW[1] );
   READ_DOUBLE( szFileName, velMW[2] );

   //take care of (in)valid pressure input
  /*if (*presDelta<=0){
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
   if (*presDelta>0){ // if pressure given, left and right bound. set to outflow (they probably already are set in the input file, but just in case)
	*wl = 3;
	*wr = 3;
}*/
   *dx = *xlength / (double)(*imax);
   *dy = *ylength / (double)(*jmax);
   *dz = *zlength / (double)(*kmax);
   return 1;
}


/**
 * The arrays U,V,W and P are initialized to the constant values UI, VI, WI and PI on
 * the whole domain.
 */
void init_uvwp(
  double UI,
  double VI,
  double WI,
  double PI,
  int imax,
  int jmax,
  int kmax,
  double ***U,
  double ***V,
  double ***W,
  double ***P,
  char * problem
) {
  init_matrix2(U, 0, imax+1, 0, jmax+1, 0, kmax+1, UI);
  init_matrix2(V, 0, imax+1, 0, jmax+1, 0, kmax+1, VI);
	init_matrix2(W, 0, imax+1, 0, jmax+1, 0, kmax+1, WI);
	init_matrix2(P, 0, imax+1, 0, jmax+1, 0, kmax+1, PI);
}


/**
 * The integer array Flag is initialized to constants C_F for fluid cells and C_B
 * for obstacle cells as specified by the parameter problem.
 */
void init_flag(
  char *problem,
  int imax,
  int jmax,
  int kmax,
//  double presDelta,
  int ***Flag,
  int wl,
  int wr,
  int wf,
  int wh,
  int wt,
  int wb
) {
	int i,j,k;
	int temp, temp2;

  //temporary array of mapping, for easier computation of flags
  int tarr[] = {1, 2, 0, 3, 4, 7, 8};
	//read the geometry
  int **Pic = read_pgm(problem);
	for (k=1; k<kmax+1; k++) {
	  //initialisation to C_F and C_B  Flag[i][0][0]
	  for (int i=1; i<jmax+1; i++){ //i is the line, j the column. so left, right is j-1,j+1, and north, south is i-1, i+1. up/down (bigger/smaller z) is -/+ (ymax+2)*k
		    for (int j=1; j<imax+1; j++){ //in Pic, 1 is where it's fluid, and 0 where it's air, apart from that: different boundary cond.
			     temp = tarr[Pic[i][j]]*pow2(2, 12) + min(Pic[i][j+1]+1, 3)*pow2(2, 10) + min(Pic[i][j-1]+1, 3)*pow2(2, 8) + min(Pic[i+1][j]+1, 3)*pow2(2, 6) +  min(Pic[i-1][j]+1, 3)*pow2(2, 4) + min(Pic[i-k*(jmax+2)][j]+1, 3)*4 + min(Pic[i+k*(jmax+2)][j]+1, 3); //use min() bcs obstacle neighbours will have numbrs 3-6 in the picture, but they should be flagged with (11)_2 = 3 in the flag field.
           //check for forbidden cells:
			     //if ( ((temp > pow(2, 12)*3) || (temp < pow(2, 12))) /*so it is bound.*/ & (E,W water / N,S water / D,U water) ) { error }
           if (Pic[i][j]>1 && min(Pic[i][j+1]+Pic[i][j-1], Pic[i-1][j]+Pic[i+1][j])!=0  && !(Pic[i-k*(jmax+2)][j]+Pic[i+k*(jmax+2)]))  {
             ERROR("Invalid geometry! Forbidden boundary cell found.\n");
           }
           temp2 = getcelltype(temp)%16;
           if (( ((temp>>12)&15) == 8 ) && (temp2==5 || temp2==6 || temp2==9 || temp2==10)){ //tryin to set moving wall for a B_??? cell. not allowed!
             ERROR("Invalid geometry! Forbidden boundary cell found.\n");
           }
           Flag[j][jmax+1-i][k] = temp;
		    }
	  } //da bo to delal more bit slika tk narjena, da ma vsaka 2D podslika okoliinokoli ghost boundary. (s poljubno boundary cifro, tj med 2 in 6)
  }
  //set outer boundary vortices, so you wont use them uninitialised. used getbit() bcs theyre on the edges of domain, so have only 2 or 3 neighbs, all boundary
  for (k=0; k<kmax+2; k++) {
    Flag[0][0][k] = getbit(wf);            //xz0
    Flag[0][jmax+1][k] = getbit(wh);       //xz1
    Flag[imax+1][0][k] = getbit(wf);       //xz0
    Flag[imax+1][jmax+1][k] = getbit(wh);  //xz1
  }
  for (j=0; j<jmax+2; j++) {
    Flag[0][j][0] = getbit(wb);            //xy0
    Flag[0][j][kmax+1] = getbit(wt);       //xy1
    Flag[imax+1][j][0] = getbit(wb);       //xy0
    Flag[imax+1][j][kmax+1] = getbit(wt);  //xy1
  }
  for (i=0; i<imax+2; i++) {
    Flag[i][0][0] = getbit(wf);            //xy0
    Flag[i][jmax+1][0] = getbit(wh);       //xy0
    Flag[i][0][kmax+1] = getbit(wf);       //xy1
    Flag[i][jmax+1][kmax+1] = getbit(wh);  //xy1
	}

  //set outer boundary flags - top and bottom:
	for (i=1; i<=imax; i++){
    for (j=1; j<=jmax; j++){
		    if ((Flag[i][j][1]&pow2(2,12))==(Flag[i][j][1]&pow2(2,13))) {//if our only inner cell neighbour is a boundary(obstacle), we set this cell to the same kind of boundary.
			       Flag[i][j][0] = getbit(0)/*3*(pow(2,10)+pow(2,8)+pow(2,6)+pow(2,4)+4+1)*/ + (Flag[i][j][1]&(pow2(2,12)*15)/*info bout the boundary at the beginning*/);
		    } else { //Otherwise set it to what you read in the picture (<- remark: NO. we set it to what we read in the parameter file.). plus subtract/set the neighbouring cell, which seems to be water or air
			       Flag[i][j][0] = (getbit(wb)-3) + (Flag[i][j][1]&3);
		    }
        if ((Flag[i][j][kmax]&pow2(2,12))==(Flag[i][j][kmax]&pow2(2,13))) {//if our only inner cell neighbour is a boundary(obstacle), we set this cell to the same kind of bc.
             Flag[i][j][kmax+1] = getbit(0) + (Flag[i][j][kmax]&(pow2(2,12)*15)); /*info bout the boundary at the beginning*/
        } else { //Otherwise set it to what you read in the picture <- remark: NO. we set it to what we read in the parameter file.
             Flag[i][j][kmax+1] = (getbit(wt)-12) + (Flag[i][j][kmax]&12);
        }
    }
	}
	//set the outer boundary flags - right and left:
  for (j=1; j<=jmax; j++){
    for (k=1; k<=kmax; k++){
        if ((Flag[1][j][k]&pow2(2,12))==(Flag[1][j][k]&pow2(2,13))) {
             Flag[0][j][k] = getbit(0) + (Flag[1][j][k]&(pow2(2,12)*15));
        } else { //Otherwise set it to what you read in the picture <- remark: NO. we set it to what we read in the parameter file.
             Flag[0][j][k] = (getbit(wl)- 3*pow2(2,10)) + (Flag[1][j][k]& (3*pow2(2,10)));
        }
        if ((Flag[imax][j][k]&pow2(2,12))==(Flag[imax][j][k]&pow2(2,13))) {
             Flag[imax+1][j][k] = getbit(0) + (Flag[imax][j][k]&(pow2(2,12)*15));
        } else {
             Flag[imax+1][j][k] = (getbit(wr)-3*pow2(2,8)) + (Flag[imax][j][k]&(3*pow2(2,8)));
        }
    }
  }
  //set the outer boundary flags - front and back:
  for (i=1; i<=imax; i++){
    for (k=1; k<=kmax; k++){
        if ((Flag[i][1][k]&pow2(2,12))==(Flag[i][1][k]&pow2(2,13))) {
             Flag[i][0][k] = getbit(0) + (Flag[i][1][k]&(pow2(2,12)*15));
        } else {
             Flag[i][0][k] = (getbit(wf)-16*3) + (Flag[i][1][k]&(16*3));
        }
        if ((Flag[i][jmax][k]&pow2(2,12))==(Flag[i][jmax][k]&pow2(2,13))) {
             Flag[i][jmax+1][k] = getbit(0) + (Flag[i][jmax][k]&(pow2(2,12)*15));
        } else {
             Flag[i][jmax+1][k] = (getbit(wh)-64*3) + (Flag[i][jmax][k]&(64*3));
        }
    }
  }
	// take care of the case when pressure is given
/*	for (i=0; i<=imax+1; i++){
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
	}*/
}
