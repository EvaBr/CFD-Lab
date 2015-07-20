#ifndef __INIT_H_
#define __INIT_H_




#include "surface.h"

/**
 * This operation initializes all the local variables reading a configuration
 * file. For every variable a macro like READ_INT() is called passing it the
 * input filename and the variable to be written to. This macro calls
 * an operation read_int() augmenting the parameter set with the name of the
 * variable to be read. The read_int() operation parses the input file, extracts
 * the value of the variable, sets the variable and finally prints some debug
 * information. This is possible as the macro adds the name of the variable to
 * be set. All the helper operations can be found within helper.h and helper.c.
 *
 * @param szFileName char pointer to the filename
 * @param Re         Reynolds number
 * @param UI         initial velocity in  x-direction - used by init_uvp()
 * @param VI         initial velocity y-direction - used by init_upv()
 * @param PI         initial pressure - used by init_upv()
 * @param GX         gravitation x-direction
 * @param GY         gravitation y-direction
 * @param t_end      end time (not discrete in time steps)
 * @param xlength    domain length x-direction
 * @param ylength    domain lenght y-direction
 * @param dt         time step length: dividing t_end by dt gives the number of
 *                   time steps to perform. Actually dt is determined by a
 *                   function, so manipulating this value within the
 *                   configuration file should not affect the solution process
 *                   at all
 * @param dx         cell length x-direction
 * @param dy         cell length y-direction
 * @param imax       number of cells in x-direction
 * @param jmax       number of cells in Y-direction
 * @param alpha      uppwind-differencing-factor alpha
 * @param omg        relaxation factor omega
 * @param tau        safety parameter for time step calculation
 * @param itermax    max. number of pressure iterations
 * @param eps        tolerance limit for pressure calculation
 * @param dt_value   time steps for output (after how many time steps one should
 *                   write into the output file)
 * @param wl	     initial boundary for left wall
 * @param wr	     initial boundary for right wall
 * @param wt	     initial boundary for top wall
 * @param wb	     initial boundary for bottom wall
 * @param problem    problem to solve
 */
int read_parameters(
		const char *szFileName,
		double *Re,
		double *UI,
		double *VI,
		double *WI,
		double *PI,
		double *GX,
		double *GY,
		double *GZ,
		double *t_end,
		double *xlength,
		double *ylength,
		double *zlength,
		double *dt,
		double *dx,
		double *dy,
		double *dz,
		int  *imax,
		int  *jmax,
		int  *kmax,
		double *alpha,
		double *omg,
		double *tau,
		int  *itermax,
		double *eps,
		double *dt_value,
		int *wl,
		int *wr,
		int *wf,
		int *wh,
		int *wt,
		int *wb,
		char *problem,
		/*double *presLeft,
  double *presRight,
  double *presDelta,*/
		double *velIN,
		double *velMW,
		int *particles
);


void init_particles(
		int ***Flag,
		double dx,
		double dy,
		double dz,
		int imax,
		int jmax,
		int kmax,
		int ppc,
		struct particleline *Partlines);



/**
 * The arrays U,V,W and P are initialized to the constant values UI, VI and PI on
 * the whole domain.
 */
void init_uvwp(
		double UI,
		double VI,
		double WI,
		double PI,
		int ***Flag,
		int imax,
		int jmax,
		int kmax,
		double ***U,
		double ***V,
		double ***W,
		double ***P,
		char * problem
);

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
);

#endif