/*
 * surface.h
 *
 */

#ifndef SURFACE_H_
#define SURFACE_H_


struct particle{
     double x;
     double y;
     double z;
     double vel;
     struct particle *next;
};

struct particleline{
     int length;
     struct particle *Particles;
};

struct cell{
	double *p;
	double *u[2];
	double *v[2];
	double *w[2];
};

void getcell(struct cell* c,double ***U,double ***V,double ***W,double ***P,int i, int j, int k,int imax, int jmax, int kmax, double* zero);

void compute_surface_values_1(double ***U,double ***V,double ***W,double ***P,int ***Flag,int i,int j,int k,int imax, int jmax, int kmax,double dx,double dy,double dz,double dt,double GX,double GY,double GZ, double Re);
void compute_surface_values_2(double ***U,double ***V,double ***W,double ***P,int ***Flag,int i,int j,int k,int imax, int jmax, int kmax,double dx,double dy,double dz,double dt,double GX,double GY,double GZ, double Re);
void compute_surface_values_3(double ***U,double ***V,double ***W,double ***P,int ***Flag,int i,int j,int k,int imax, int jmax, int kmax,double dx,double dy,double dz,double dt,double GX,double GY,double GZ, double Re);


void compute_surface_pressure(double ***U,double ***V,double ***W,double ***P,int ***Flag,int i,int j,int k,int imax, int jmax, int kmax,double dx,double dy,double dz,double dt,double GX,double GY,double GZ, double Re);




double trilinearInterpolation(double ***m,int i, int j, int k, double dx, double dz, double dy,double x,double y, double z, double x1,  double x2,double y1,double y2,double z1,double z2);

void add_cell_particles(struct particleline *Partline, double i, double j, double k,double dx,double dy,double dz,int ppc);

struct particle *create_particle(double x, double y, double z);

void add_particle(struct particleline *Partline, double x,double y,double z);

void mark_cells(int ***Flag,double dx,double dy, double dz,int imax,int jmax,int kmax,int N,struct particleline *Partlines);

void set_uvwp_surface(double ***U,double ***V,double ***W,double ***P,int ***Flag,double dx,double dy,double dz, int imax,int jmax, int kmax,double GX,double GY,double GZ,double dt,double Re);

void advance_particles(double dx,double dy, double dz,int imax,int jmax,int kmax, double dt,double ***U,double ***V,double ***W,int N,struct particleline *Partlines);

void get_particle_speed(double dx,double dy, double dz,int imax,int jmax,int kmax, double dt,double ***U,double ***V,double ***W,int N,struct particleline *Partlines);


#endif /* SURFACE_H_ */
