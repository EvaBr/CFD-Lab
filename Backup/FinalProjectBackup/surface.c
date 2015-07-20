#include "surface.h"
#include "helper.h"

void add_cell_particles(struct particleline *Partline, double i, double j, double k,double dx,double dy,double dz,int ppc){
	int ic,jc,kc;
	double x,y,z;
	for (ic=1;ic<=ppc;ic++){
		x = (i-1)*dx+(ic-0.5)/((double)ppc)*dx;
		for (jc=1;jc<=ppc;jc++)
		{
			y = (j-1)*dy+(jc-0.5)/((double)ppc)*dy;
			for (kc=1;kc<=ppc;kc++)
			{
				z = (k-1)*dz+(kc-0.5)/((double)ppc)*dz;
				add_particle(Partline,x,y,z);
			}
		}
	}
}


struct particle *create_particle(double x, double y, double z){
	struct particle *p  = (struct particle *) malloc(sizeof(struct particle));
	if(p == 0)
	{
		ERROR("Could not create particle!");
	}
	p->next = 0;
	p->x   = x;
	p->y   = y;
	p->z   = z;
	p->vel = 0;
	return  p;
}


void add_particle(struct particleline *Partline, double x,double y,double z)
{
	struct particle *p;
	p =   create_particle(x, y, z);
	p->next = (*Partline).Particles->next;
	(*Partline).Particles->next = p;
	(*Partline).length++;
}





void mark_cells(int ***Flag,double dx,double dy, double dz,int imax,int jmax,int kmax,int N,struct particleline *Partlines){
	int i,j,k,n;
	struct particle *t,*p;
	double x,y,z;

	for (i=0;i<=imax+1;i++){
		for (j=0;j<=jmax+1;j++){
			for(k=0;k<=kmax+1;k++){
				if (!isboundary(Flag[i][j][k])){
					setcelltype(&Flag[i][j][k],AIR);
				}
			}
		}
	}
	for (n=0;n<N;n++){
		for(p=Partlines[n].Particles; p->next != 0; p=p->next)
		{
			x = p->next->x;
			y = p->next->y;
			z = p->next->z;
			i = (int)(x/dx)+1;
			j = (int)(y/dy)+1;
			k = (int)(z/dz)+1;
			if (isboundary(Flag[i][j][k]))
			{
				t = p->next->next;
				free(p->next);
				p->next = t;
				Partlines[n].length--;
			}
			else{
				setcelltype(&Flag[i][j][k],FLUID);
			}

		}
	}

	for (j=1;j<=jmax;j++){
		for (i=1;i<=imax;i++)
		{
			for (i=1;i<=imax;i++)
			{

			}
		}
	}


}


void set_uvwp_surface(double ***U,double ***V,double ***W,double ***P,int ***Flag,double dx,double dy,double dz, int imax,int jmax, int kmax,double GX,double GY,double GZ,double Re){

}

inline double trilinearInterpolation(double ***m,int i, int j, int k, double dx, double dz, double dy,double x,double y, double z, double x1,  double x2,double y1,double y2,double z1,double z2){
	double t1,t2;

	t1 = ((y2-y)*((x2-x)*m[i-1][j-1][k]  +(x-x1)*m[i][j-1][k])   +
			(y-y1)*((x2-x)*m[i-1][j][k]    +(x-x1)*m[i][j][k])
	)/dx/dy;

	t2 = ((y2-y)*((x2-x)*m[i-1][j-1][k-1] + (x-x1)*m[i][j-1][k-1])   +
				(y-y1)*((x2-x)*m[i-1][j][k-1]   + (x-x1)*m[i][j][k-1])
		)/dx/dy;

	return((z2-z)*t2+(z-z1)*t1)/dz;

}

void advance_particles(double dx,double dy, double dz,int imax,int jmax,int kmax, double dt,double ***U,double ***V,double ***W,int N,struct particleline *Partlines){

	int i, j, k,s;
	int i1, j1, k1;
	int i2, j2, k2;
	double x,y,z;
	double x1_1,x1_2,x2_1,x2_2;
	double y1_1,y1_2,y2_1,y2_2;
	double z1_1,z1_2,z2_1,z2_2;

	double u;
	double v;
	double w;

	struct particle *p;
	struct particle *t;


	for(s=0;s<N;s++){
		for(p=Partlines[s].Particles; p->next != NULL; p=p->next){
			x = p->next->x;
			y = p->next->y;
			z = p->next->z;

			i1 = (int)(x/dx)+1;
			j1 = (int)(y/dy)+1;
			k1 = (int)(z/dz)+1;

			i2 = (int)((x+0.5*dx)/dx)+1;
			j2 = (int)((y+0.5*dy)/dy)+1;
			k2 = (int)((z+0.5*dz)/dz)+1;

			if (i1<=0 || j1<=0 || k1<=0 || i2>imax || j2>jmax || k2>kmax ){

				t = p->next;
				p->next = p->next->next;
				free(t);
				Partlines[s].length--;
				if (p->next==NULL)
					break;
			}
			else
			{
				i = i1;
				j = j2;
				k = k2;

				x1_1 = (i-1)*dx;
				y1_1 = (j-1)*dy;
				z1_1 = (k-1)*dz;

				x1_2 = ((i-1)-0.5)*dx;
				y1_2 = ((j-1)-0.5)*dy;
				z1_2 = ((k-1)-0.5)*dz;

				x2_1 = i*dx;
				y2_1 = j*dy;
				z2_1 = k*dz;

				x2_2 = (i-0.5)*dx;
				y2_2 = (j-0.5)*dy;
				z2_2 = (k-0.5)*dz;

				u = trilinearInterpolation(U,i1, j2, k2,  dx, dz, dy,x,y,z,  x1_1, x2_1,y1_2,y2_2,z1_2,z2_2);

				i = i2;
				j = j1;
				k = k2;

				x1_1 = (i-1)*dx;
				y1_1 = (j-1)*dy;
				z1_1 = (k-1)*dz;

				x1_2 = ((i-1)-0.5)*dx;
				y1_2 = ((j-1)-0.5)*dy;
				z1_2 = ((k-1)-0.5)*dz;

				x2_1 = i*dx;
				y2_1 = j*dy;
				z2_1 = k*dz;

				x2_2 = (i-0.5)*dx;
				y2_2 = (j-0.5)*dy;
				z2_2 = (k-0.5)*dz;


				v = trilinearInterpolation(V,i2, j1, k2,  dx, dz, dy,x,y,z,  x1_2, x2_2,y1_1,y2_1,z1_2,z2_2);

				i = i2;
				j = j2;
				k = k1;

				x1_1 = (i-1)*dx;
				y1_1 = (j-1)*dy;
				z1_1 = (k-1)*dz;

				x1_2 = ((i-1)-0.5)*dx;
				y1_2 = ((j-1)-0.5)*dy;
				z1_2 = ((k-1)-0.5)*dz;

				x2_1 = i*dx;
				y2_1 = j*dy;
				z2_1 = k*dz;

				x2_2 = (i-0.5)*dx;
				y2_2 = (j-0.5)*dy;
				z2_2 = (k-0.5)*dz;


				w = trilinearInterpolation(W,i2, j2, k1,  dx, dz, dy,x,y,z,  x1_2, x2_2,y1_2,y2_2,z1_1,z2_1);


				x += dt*u;
				y += dt*v;
				z += dt*w;

				i = (int)(x/dx);
				j = (int)(y/dy);
				k = (int)(z/dz);

				if (i<0 || j<0 || k<0 || i>imax || j>jmax || k>kmax ){
					t = p->next;
					p->next = p->next->next;
					free(t);
					Partlines[s].length--;
					if (p->next==NULL)
						break;

				}
				else{

					/* todo boundary */

					p->next->x = x;
					p->next->y = y;
					p->next->z = z;
					p->next->vel = 2*sqrt(u*u+v*v*w*w);
				}
			}
		}
	}
}
