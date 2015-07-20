#ifndef __VISUAL_H__
#define __VISUAL_H__

#include "surface.h"
#include <stdio.h>

void write_vtkFile_debug(const char *szProblem,
		  int    timeStepNumber,
		  double xlength,
      double ylength,
			double zlength,
      int    imax,
      int    jmax,
			int    kmax,
		  double dx,
		  double dy,
			double dz,
      double ***U,
      double ***V,
			double ***W,
      double ***P,
      int ***Flag);

/**
 * Method for writing header information in vtk format.
 *
 * @param fp      File pointer for writing info.
 * @param imax    Maximum number of entries (minus 2) in x-direction
 * @param jmax    Maximum number of entries (minus 2) in y-direction
 * @param dx      mesh size dx
 * @param dy      mesh size dy
 *
 * @author Tobias Neckel
 */
void write_vtkHeader_debug( FILE *fp, int imax, int jmax, int kmax,
                      double dx, double dy, double dz);

/**
 * Method for writing grid coordinate information in vtk format.
 *
 * @param fp      File pointer for writing info.
 * @param imax    Maximum number of entries (minus 2) in x-direction
 * @param jmax    Maximum number of entries (minus 2) in y-direction
 * @param dx      mesh size dx
 * @param dy      mesh size dy
 *
 * @author Tobias Neckel
 */
void write_vtkPointCoordinates_debug( FILE *fp, int imax, int jmax, int kmax,
                                double dx, double dy, double dz);











void write_particles_debug(const char *szProblem, int    timeStepNumber, int N, struct particleline *Partlines);




/**
 * Method for writing header information in vtk format.
 *
 * The name of the file consists of the problem name (szProblem)
 * and of the current time step. It gets the suffix .vtk.
 *
 * @param szProblem      File pointer for writing info.
 * @param timeStepNumber Number of the current time step to be printed.
 * @param xlength Length in x-direction
 * @param ylength Length in y-direction
 * @param imax    Maximum number of entries (?) in x-direction
 * @param jmax    Maximum number of entries (?) in y-direction
 * @param dx      Mesh size in x-direction
 * @param dy      Mesh size in x-direction
 * @param U       Velocities in x-direction
 * @param V       Velocities in y-direction
 * @param P       Pressure data
 *
 * @author Tobias Neckel
 */
void write_vtkFile_debug(const char *szProblem,
		  int    timeStepNumber,
		  double xlength,
      double ylength,
			double zlength,
      int    imax,
      int    jmax,
			int    kmax,
		  double dx,
		  double dy,
			double dz,
      double ***U,
      double ***V,
			double ***W,
      double ***P,
      int ***Flag);

/**
 * Method for writing header information in vtk format.
 *
 * @param fp      File pointer for writing info.
 * @param imax    Maximum number of entries (minus 2) in x-direction
 * @param jmax    Maximum number of entries (minus 2) in y-direction
 * @param dx      mesh size dx
 * @param dy      mesh size dy
 *
 * @author Tobias Neckel
 */
void write_vtkHeader_debug( FILE *fp, int imax, int jmax, int kmax,
                      double dx, double dy, double dz);

/**
 * Method for writing grid coordinate information in vtk format.
 *
 * @param fp      File pointer for writing info.
 * @param imax    Maximum number of entries (minus 2) in x-direction
 * @param jmax    Maximum number of entries (minus 2) in y-direction
 * @param dx      mesh size dx
 * @param dy      mesh size dy
 *
 * @author Tobias Neckel
 */
void write_vtkPointCoordinates_debug( FILE *fp, int imax, int jmax, int kmax,
                                double dx, double dy, double dz);




void write_particles(const char *szProblem, int    timeStepNumber, int N, struct particleline *Partlines);




/**
 * Method for writing header information in vtk format.
 *
 * The name of the file consists of the problem name (szProblem)
 * and of the current time step. It gets the suffix .vtk.
 *
 * @param szProblem      File pointer for writing info.
 * @param timeStepNumber Number of the current time step to be printed.
 * @param xlength Length in x-direction
 * @param ylength Length in y-direction
 * @param imax    Maximum number of entries (?) in x-direction
 * @param jmax    Maximum number of entries (?) in y-direction
 * @param dx      Mesh size in x-direction
 * @param dy      Mesh size in x-direction
 * @param U       Velocities in x-direction
 * @param V       Velocities in y-direction
 * @param P       Pressure data
 *
 * @author Tobias Neckel
 */
void write_vtkFile(const char *szProblem,
		  int    timeStepNumber,
		  double xlength,
      double ylength,
			double zlength,
      int    imax,
      int    jmax,
			int    kmax,
		  double dx,
		  double dy,
			double dz,
      double ***U,
      double ***V,
			double ***W,
      double ***P,
      int ***Flag);

/**
 * Method for writing header information in vtk format.
 *
 * @param fp      File pointer for writing info.
 * @param imax    Maximum number of entries (minus 2) in x-direction
 * @param jmax    Maximum number of entries (minus 2) in y-direction
 * @param dx      mesh size dx
 * @param dy      mesh size dy
 *
 * @author Tobias Neckel
 */
void write_vtkHeader( FILE *fp, int imax, int jmax, int kmax,
                      double dx, double dy, double dz);

/**
 * Method for writing grid coordinate information in vtk format.
 *
 * @param fp      File pointer for writing info.
 * @param imax    Maximum number of entries (minus 2) in x-direction
 * @param jmax    Maximum number of entries (minus 2) in y-direction
 * @param dx      mesh size dx
 * @param dy      mesh size dy
 *
 * @author Tobias Neckel
 */
void write_vtkPointCoordinates( FILE *fp, int imax, int jmax, int kmax,
                                double dx, double dy, double dz);

#endif
