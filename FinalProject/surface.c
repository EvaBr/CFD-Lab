#include "surface.h"
#include "helper.h"
#include <omp.h>

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

void clean_cells(int ***Flag,int imax,int jmax, int kmax){
	int i,j,k;
#pragma omp parallel for private(j,k)
	for (i=0;i<=imax+1;i++){
		for (j=0;j<=jmax+1;j++){
			for(k=0;k<=kmax+1;k++){
				if (!isboundary(Flag[i][j][k])){
					setcelltype(&Flag[i][j][k],AIR);

				}
			}
		}
	}
}

void update_cell(int ***Flag,int i, int j, int k ){

	int *flag=0,*flagc=0;
	flagc = &Flag[i][j  ][k  ];
	//[east][west][north][south][bottom][top]
	if(!isboundary(*flagc) )
	{
		flag = &Flag[i+1][j  ][k  ];
		if(isempty(*flag)){
			changebit(flagc,10,0);//B_O
			changebit(flagc,11,1);

		}
		else if (isfluid(*flag)){
			changebit(flagc,10,1);
			changebit(flagc,11,0);

		}


		flag = &Flag[i-1][j  ][k  ];
		if(isempty(*flag)){
			changebit(flagc,8,0); //B_W
			changebit(flagc,9,1);
		}
		else if (isfluid(*flag)){
			changebit(flagc,8,1);
			changebit(flagc,9,0);

		}

		flag = &Flag[i][j+1][k  ];
		if(isempty(*flag)){
			//printf("markcell!\n");
			changebit(flagc,6,0); //B_N
			changebit(flagc,7,1);
		}
		else if(isfluid(*flag)){
			changebit(flagc,6,1);
			changebit(flagc,7,0);

		}

		flag = &Flag[i][j-1][k  ];
		if(isempty(*flag)){
			//printf("markcell!\n");
			changebit(flagc,4,0); //B_S
			changebit(flagc,5,1);
		}
		else if(isfluid(*flag)){
			changebit(flagc,4,1);
			changebit(flagc,5,0);

		}

		flag = &Flag[i][j  ][k-1];
		if(isempty(*flag)){

			changebit(flagc,2,0); //B_D
			changebit(flagc,3,1);
		}
		else if (isfluid(*flag)){

			changebit(flagc,2,1);
			changebit(flagc,3,0);
		}

		flag = &Flag[i][j  ][k+1];
		if (isfluid(*flag)){
			changebit(flagc,0,1);
			changebit(flagc,1,0);
		}
		else if(isempty(*flag)){
			changebit(flagc,0,0);//B_U
			changebit(flagc,1,1);
		}
	}
}


void mark_cells(int ***Flag,double dx,double dy, double dz,int imax,int jmax,int kmax,int N,struct particleline *Partlines,struct p_pointer *PP1,int*FluidCells){
	int i,j,k,n;
	struct particle *p;//*t,
	double x,y,z;
	int PP1counter = 0;
	clean_cells(Flag,imax,jmax, kmax);
	for (n=0;n<N;n++){
		for(p=Partlines[n].Particles; p->next != 0; p=p->next)
		{
			x = p->next->x;
			y = p->next->y;
			z = p->next->z;
			i = (int)(x/dx)+1;
			j = (int)(y/dy)+1;
			k = (int)(z/dz)+1;
			if (i<1 || j<1 || k<1 || i>imax || j>jmax || k>kmax){
				/*
				t = p->next;
				p->next = p->next->next;
				free(t);
				Partlines[n].length--;
				if (p->next==NULL)
					break;
				 */
			}
			else{
				if (!isboundary(Flag[i][j][k])){
					setcelltype(&Flag[i][j][k],FLUID);


				}
			}
		}
	}



#pragma omp parallel for private(j,k)
	for (i=1;i<=imax;i++)
	{
		for (j=1;j<=jmax;j++){
			for (k=1;k<=kmax;k++)
			{
				update_cell(Flag,i,j,k);
			}
		}

	}

	int found_p = 0;
	int found_u = 0;
	int found_v = 0;
	int found_w = 0;
	for (i=1;i<=imax;i++)
	{
		for (j=1;j<=jmax;j++){
			for (k=1;k<=kmax;k++)
			{
				found_p = 0;
				found_u = 0;
				found_v = 0;
				found_w = 0;
				if(isfluid(Flag[i][j][k])){

					if(!emptyneighbor(Flag[i][j][k])){
						found_p = 1;
					}
					if(isfluid(Flag[i+1][j][k])){
						found_u = 1;
					}
					if(isfluid(Flag[i][j+1][k])){
						found_v = 1;
					}
					if(isfluid(Flag[i][j][k+1])){
						found_w = 1;
					}
				}


				if(found_p||found_u||found_v||found_w){
					PP1[PP1counter].i = i;
					PP1[PP1counter].j = j;
					PP1[PP1counter].k = k;
					PP1[PP1counter].p = found_p;
					PP1[PP1counter].u = found_u;
					PP1[PP1counter].v = found_v;
					PP1[PP1counter].w = found_w;
					PP1counter++;
				}
			}

		}
	}

	*FluidCells = PP1counter;
	PP1[PP1counter].i = -1;

}

void clean_empty_space(int*** Flag,double ***U,double ***V,double ***W,double ***P,int imax,int jmax, int kmax){
	int i,j,k;
#pragma omp parallel for private(j,k)
	for (i=0;i<=imax-1;i++){
		for (j=0;j<=jmax;j++){
			for (k=0;k<=kmax;k++){
				if (!isfluid(Flag[i][j][k]) && !isfluid(Flag[i+1][j][k])){
					U[i][j][k] = 0;
				}
			}
		}
	}
#pragma omp parallel for private(j,k)
	for (i=0;i<=imax-1;i++){
		for (j=0;j<=jmax-1;j++){
			for (k=0;k<=kmax;k++){
				if (!isfluid(Flag[i][j][k]) && !isfluid(Flag[i][j+1][k])){
					V[i][j][k] = 0;
				}
			}
		}
	}
#pragma omp parallel for private(j,k)
	for (i=0;i<=imax;i++){
		for (j=0;j<=jmax;j++){
			for (k=0;k<=kmax-1;k++){
				if (!isfluid(Flag[i][j][k]) && !isfluid(Flag[i][j][k+1])){
					W[i][j][k] = 0;
				}
			}
		}
	}
#pragma omp parallel for private(j,k)
	for (i=1;i<=imax;i++){
		for (j=1;j<=jmax;j++){
			for (k=1;k<=kmax;k++){
				if (!isempty(Flag[i][j][k] )) {
					P[i][j][k] = 0;
				}
			}
		}
	}

}

inline void getcell(struct cell* c,double ***U,double ***V,double ***W,double ***P,int i, int j, int k,int imax, int jmax, int kmax, double* zero){
	c->p    = &P[i][j][k];
	if(i>0)      c->u[0] = &U[i-1][j][k];   else c->u[0] = zero;
	if(i<=imax+1)c->u[1] = &U[i][j][k];     else c->u[1] = zero;

	if(j>0)       c->v[0] = &V[i][j-1][k];  else c->v[0] = zero;
	if(j<=jmax+1) c->v[1] = &V[i][j][k];    else c->v[1] = zero;

	if(k>0)       c->w[0] = &W[i][j][k-1];  else c->w[0] = zero;
	if(k<=kmax+1) c->w[1] = &W[i][j][k];    else c->w[1] = zero;
}


inline void compute_surface_values_1(double ***U,double ***V,double ***W,double ***P,int ***Flag,int i,int j,int k,int imax, int jmax, int kmax,double dx,double dy,double dz,double dt,double GX,double GY,double GZ, double Re){
	//int type;
	double zero = 0;
	int nx=0,ny=0,nz=0;
	int xdb=0,ydb=0,zdb=0;
	int ui=0,vi=0,wi=0;
	int num = 0;

	getsurfacetype(Flag[i][j][k],&ui,&vi,&wi,&nx,&ny,&nz,&xdb,&ydb,&zdb,&num);
	//printf("type: %d %d,%d,%d\n",type,i,j,k);
	struct cell cells[3][3][3] ;
	struct cell* c = &cells[1][1][1];


	getcell(c,U,V,W,P,i,j,k,imax, jmax,kmax,&zero);


	/*double ddu=0,ddv=0,ddk=0;*/

	double du1=0,dv1=0,dk1=0;

	if(nx==0){
		du1 = ((*c->u[1]-*c->u[0])/dx);
	}
	if(ny==0){
		dv1 = ((*c->v[1]-*c->v[0])/dy);
	}
	if(nz==0){
		dk1 = ((*c->w[1]-*c->w[0])/dz);
	}

	if(xdb==1){
		*c->u[0] = *c->u[0] + dt*GX;
		*c->u[1] = *c->u[1] + dt*GX;
	}
	else if(nx!=0){
		*c->u[ui] = *c->u[1-ui]-((dv1+dk1)*dx)*nx;

		/*printf("surface_1 (%d,%d,%d) %f, %f, %f (%f)-> u = %f  \n",i,j,k,
		 *c->u[1-ui],
				dv1,
				dk1,dx,
		 *c->u[ui]);
		if(dk1!=0){
			printf("dk1: %f, %f, %f, %f\n",*c->w[1],*c->w[0],dz,1/dz);
		}
		 */

	}

	if(ydb==1){

		*c->v[0] = *c->v[0] + dt*GY;
		*c->v[1] = *c->v[1] + dt*GY;
	}
	else if(ny!=0){
		*c->v[vi] = *c->v[1-vi]-((du1+dk1)*dy)*ny;
	}

	if(zdb==1){
		*c->w[0] = *c->w[0]+dt*GZ;
		*c->w[1] = *c->w[1]+dt*GZ;

	}
	else if(nz!=0){
		*c->w[wi] = *c->w[1-wi]-((du1+dv1)*dz)*nz;

	}
}

inline void compute_surface_values_2(double ***U,double ***V,double ***W,double ***P,int ***Flag,int i,int j,int k,int imax, int jmax, int kmax,double dx,double dy,double dz,double dt,double GX,double GY,double GZ, double Re){
	//int type;
	double zero = 0;
	int nx=0,ny=0,nz=0;
	int xdb=0,ydb=0,zdb=0;
	int ui=0,vi=0,wi=0;
	int num = 0;
	int mode = 0;
	getsurfacetype(Flag[i][j][k],&ui,&vi,&wi,&nx,&ny,&nz,&xdb,&ydb,&zdb,&num);

	//printf("type: %d %d,%d,%d\n",type,i,j,k);
	struct cell cells[3][3][3] ;
	struct cell* c = &cells[1][1][1];

	getcell(c,U,V,W,P,i,j,k,imax, jmax,kmax,&zero);


	if(nx!=0 && xdb==0){


		getcell(&cells[1+nx][1   ][1   ],U,V,W,P,i+nx,j   ,k   ,imax, jmax,kmax,&zero);
		getcell(&cells[1   ][0   ][1   ],U,V,W,P,i   ,j-1 ,k   ,imax, jmax,kmax,&zero);
		getcell(&cells[1   ][1   ][0   ],U,V,W,P,i   ,j   ,k-1,imax, jmax,kmax,&zero);

		mode = 0;
		if(j>1){
			if(isboundary(Flag[i+nx][j   ][j-2])){
				mode = 1;
			}
		}
		if(!isempty(Flag[i+nx][j-1][k   ])){
			mode = 2;
		}
		if(mode ==0){
			*cells[1+nx][1   ][1   ].v[0] =*c->v[0]-(dx/dy)*(*c->u[ui]-*cells[1][0][1].u[ui]) *nx;
		}
		else if(mode==1){
			*cells[1+nx][1   ][1   ].v[0] =*c->v[0];
		}

		mode = 0;
		if(k>1){
			if(isboundary(Flag[i+nx][j   ][k-2])){
				mode = 1;
			}
		}
		if(!isempty(Flag[i+nx][j   ][k-1])){
			mode = 2;
		}
		if(mode==0){
			*cells[1+nx][1   ][1   ].w[0] =*c->w[0]-(dx/dz)*(*c->u[ui]-*cells[1][1][0].u[ui]) *nx;
			/*printf("(%d,%d,%d)   %f, %f, %f -> w = %f  \n\n",i,j,k,
			 *c->w[0],
			 *c->u[ui],
			 *cells[1][1][0].u[ui],
			 *cells[1+nx][1   ][1   ].w[0]);*/
		}
		else if(mode==1){
			*cells[1+nx][1   ][1   ].w[0] =*c->w[0];
		}



	}


	if(ny!=0&& ydb == 0){

		getcell(&cells[1   ][1+ny][1   ],U,V,W,P,i   ,j+ny,k   ,imax, jmax,kmax,&zero);
		getcell(&cells[0   ][1   ][1   ],U,V,W,P,i-1 ,j   ,k   ,imax, jmax,kmax,&zero);
		getcell(&cells[1   ][1   ][0   ],U,V,W,P,i   ,j   ,k-1 ,imax, jmax,kmax,&zero);


		mode = 0;
		if(i>1){
			if(isboundary(Flag[i-2][j+ny][k   ])){
				mode = 1;
			}
		}
		if(!isempty(Flag[i-1][j+ny][k   ])){
			mode = 2;
		}
		if(mode==0){
			*cells[1   ][1+ny][1   ].u[0] =*c->u[0]-(dy/dx)*(*c->v[vi]-*cells[0][1][1].v[vi])*ny;
		}
		else if(mode==1){
			*cells[1   ][1 +ny][1   ].u[0] =*c->u[0];
		}


		mode = 0;
		if(k>1){
			if(isboundary(Flag[i   ][j+ny][k-2])){
				mode = 1;
			}
		}
		if(!isempty(Flag[i   ][j+ny][k-1])){
			mode = 2;
		}
		if(mode==0){
			*cells[1   ][1+ny][1   ].w[0] =*c->w[0]-(dy/dz)*(*c->v[vi]-*cells[1][1][0].v[vi])*ny;
		}
		else if(mode==1){
			*cells[1   ][1+ny][1   ].w[0] =*c->w[0];
		}


	}

	if(nz!=0 && zdb==0){
		getcell(&cells[1   ][1   ][1+nz],U,V,W,P,i   ,j   ,k+nz,imax, jmax,kmax,&zero);
		getcell(&cells[0   ][1   ][1   ],U,V,W,P,i-1,j   ,k   ,imax, jmax,kmax,&zero);
		getcell(&cells[1   ][0   ][1   ],U,V,W,P,i   ,j-1,k   ,imax, jmax,kmax,&zero);

		mode = 0;
		if(i>1){
			if(isboundary(Flag[i-2  ][j   ][k+nz])){
				mode = 1;
			}
		}
		if(!isempty(Flag[i-1  ][j   ][k+nz])){
			mode = 2;
		}
		if(mode==0){
			*cells[1   ][1   ][1+nz].u[0] =*c->u[0]-(dz/dx)*(*c->w[wi]-*cells[0][1][1].w[wi])*nz;
		}
		else if(mode==1){
			*cells[1   ][1   ][1+nz].u[0] =*c->u[0];
		}

		mode = 0;
		if(j>1){
			if(isboundary(Flag[i   ][j-2][k+nz])){
				mode = 1;
			}
		}
		if(!isempty(Flag[i   ][j-1][k+nz])){
			mode = 2;
		}
		if(mode==0){

			*cells[1   ][1   ][1+nz].v[0] =*c->v[0]-(dz/dy)*(*c->w[wi]-*cells[1][0][1].w[wi])*nz;
		}
		else if(mode==1){
			*cells[1   ][1   ][1+nz].v[0] =*c->v[0];
		}
	}

}

inline void compute_surface_values_3(double ***U,double ***V,double ***W,double ***P,int ***Flag,int i,int j,int k,int imax, int jmax, int kmax,double dx,double dy,double dz,double dt,double GX,double GY,double GZ, double Re){
	//int type;
	double zero = 0;
	int nx=0,ny=0,nz=0;
	int xdb=0,ydb=0,zdb=0;
	int ui=0,vi=0,wi=0;
	int num = 0;

	getsurfacetype(Flag[i][j][k],&ui,&vi,&wi,&nx,&ny,&nz,&xdb,&ydb,&zdb,&num);

	if(num>1){
		//printf("%d,%d,%d,NUM > 2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",i,j,k);
		struct cell cells[3][3][3] ;
		struct cell* c = &cells[1][1][1];
		getcell(c,U,V,W,P,i,j,k,imax, jmax,kmax,&zero);

		if(nx!=0 && xdb == 0){
			*c->u[ui] =  *c->u[1-ui];
		}

		if(ny!=0 && ydb == 0){
			*c->v[vi] =  *c->v[1-vi];
		}

		if(nz!=0 && zdb == 0){
			*c->w[wi] =  *c->w[1-wi];
		}

		if(nx!=0 && ny!=0 && xdb == 0 && ydb ==0){
			getcell(&cells[1+nx][1   ][1   ],U,V,W,P,i+nx,j   ,k   ,imax, jmax,kmax,&zero);
			getcell(&cells[1-nx][1   ][1   ],U,V,W,P,i-nx,j   ,k   ,imax, jmax,kmax,&zero);
			getcell(&cells[1   ][1+ny][1   ],U,V,W,P,i   ,j+ny,k   ,imax, jmax,kmax,&zero);
			getcell(&cells[1   ][1-ny][1   ],U,V,W,P,i   ,j-ny,k   ,imax, jmax,kmax,&zero);

			if(isempty(Flag[i-nx][j+ny][k   ])){
				*cells[1   ][1+ny][1   ].u[1-ui] =*c->u[1-ui]-(dy/dx)*(*c->v[vi]-*cells[1-nx][1][1].v[vi]);
			}
			if(isempty(Flag[i+nx][j-ny][k   ])){
				*cells[1+nx][1   ][1   ].v[1-vi] =*c->v[1-vi]-(dx/dy)*(*c->u[ui]-*cells[1][1-ny][1].u[ui]);
			}

			*cells[1   ][1+ny][1   ].u[ui] =*c->u[ui];
			*cells[1+nx][1   ][1   ].v[vi] =*c->v[vi];

		}
		if(ny!=0 && nz!=0 && ydb == 0 && zdb ==0 ){
			getcell(&cells[1   ][1+ny][1   ],U,V,W,P,i   ,j+ny,k   ,imax, jmax,kmax,&zero);
			getcell(&cells[1   ][1-ny][1   ],U,V,W,P,i   ,j-ny,k   ,imax, jmax,kmax,&zero);
			getcell(&cells[1   ][1   ][1+nz],U,V,W,P,i   ,j   ,k+nz,imax, jmax,kmax,&zero);
			getcell(&cells[1   ][1   ][1-nz],U,V,W,P,i   ,j   ,k-nz,imax, jmax,kmax,&zero);

			if(isempty(Flag[i  ][j+ny][k-nz])){
				*cells[1   ][1+ny][1   ].w[1-wi] =*c->w[1-wi]-(dy/dz)*(*c->v[vi]-*cells[1][1][1-nz].v[vi]);
			}
			if(isempty(Flag[i   ][j-ny][k+nz])){
				*cells[1   ][1   ][1+nz  ].v[1-vi] =*c->v[1-vi]-(dz/dy)*(*c->w[wi]-*cells[1][1-ny][1].w[wi]);
			}

			*cells[1   ][1   ][1+nz].v[vi]   =*c->v[vi];
			*cells[1   ][1+ny][1   ].w[wi]   =*c->w[wi];



		}
		if(nx!=0 && nz!=0 && xdb == 0 && zdb ==0){
			getcell(&cells[1+nx][1   ][1   ],U,V,W,P,i+nx,j   ,k   ,imax, jmax,kmax,&zero);
			getcell(&cells[1-nx][1   ][1   ],U,V,W,P,i-nx,j   ,k   ,imax, jmax,kmax,&zero);
			getcell(&cells[1   ][1   ][1+nz],U,V,W,P,i   ,j   ,k+nz,imax, jmax,kmax,&zero);
			getcell(&cells[1   ][1   ][1-nz],U,V,W,P,i   ,j   ,k-nz,imax, jmax,kmax,&zero);


			if(isempty(Flag[i-nx][j   ][k+nz])){
				*cells[1   ][1   ][1+nz].u[1-ui] =*c->u[1-ui]-(dz/dx)*(*c->w[wi]-*cells[1-nx][1][1].w[vi]);
			}

			if(isempty(Flag[i+nx][j   ][k-nz])){
				*cells[1+nx][1][1  ].w[1-wi] =*c->w[1-wi]-(dx/dz)*(*c->u[ui]-*cells[1][1][1-nz].u[ui]);
			}

			*cells[1   ][1   ][1+nz].u[ui]   =*c->u[ui];
			*cells[1+nx][1   ][1   ].w[wi]   =*c->w[wi];


		}
	}
}


inline void compute_surface_pressure(double ***U,double ***V,double ***W,double ***P,int ***Flag,int i,int j,int k,int imax, int jmax, int kmax,double dx,double dy,double dz,double dt,double GX,double GY,double GZ, double Re){
	//int type;
	double zero = 0;
	int nx=0,ny=0,nz=0;
	int xdb=0,ydb=0,zdb=0;
	int ui=0,vi=0,wi=0;
	int num = 0;

	getsurfacetype(Flag[i][j][k],&ui,&vi,&wi,&nx,&ny,&nz,&xdb,&ydb,&zdb,&num);
	//if(i== 1 && j == 12 && k == 2){printf("type: %d %d,%d,%d\n",type,i,j,k);};

	struct cell cells[3][3][3] ;
	struct cell* c = &cells[1][1][1];
	getcell(c,U,V,W,P,i,j,k,imax,jmax,kmax,&zero);
	double du=0,dv=0,dk=0;
	double du2=0,dv2=0,dk2=0;
	//if(i== 1 && j == 12 && k == 2){printf("set du, dv, dk\n");};
	du = ((*c->u[1]-*c->u[0])/dx);

	dv = ((*c->v[1]-*c->v[0])/dy);

	dk = ((*c->w[1]-*c->w[0])/dz);

	//if(i== 1 && j == 12 && k == 2){printf("type: %d %d,%d,%d\n",type,i,j,k);}

	double np = 0;
	//if(i== 1 && j == 12 && k == 2){printf("nx!= ? 0\n");};
	if(nx!=0){
		du2 = du;

	}
	//if(i== 1 && j == 12 && k == 2){printf("nx!= ? 0\n");};
	if(ny!=0){
		dv2 = dv;

	}
	//if(i== 1 && j == 12 && k == 2){printf("nz!= ? 0\n");};
	if(nz!=0){
		dk2 = dk;
	}
	//if(i== 1 && j == 12 && k == 2){printf("calc p\n");};
	if (xdb==1||ydb==1||zdb==1){
		/*printf("-> cell <-\n");*/
		*c->p = 0;
	}
	else{
		if(num==1){
			//if(i== 1 && j == 12 && k == 2){printf("num=1");};
			*c->p = (2.0/Re)*(du2+dv2+dk2);
		}
		else {
			np = 0;
			if(nz!=0 && ny!=0){

				getcell(&cells[1][1-ny][1],U,V,W,P,i,j-ny,k,imax,jmax,kmax,&zero);
				getcell(&cells[1][1][1-nz],U,V,W,P,i,j,k-nz,imax,jmax,kmax,&zero);
				//if(i== 1 && j == 12 && k == 2){printf("ny!=0 && nz!=0");};
				np = -(1.0/(2.0*Re))*
						(
								((*cells[1][1][1-nz].v[1] + *cells[1][1][1-nz].v[0] - *cells[1][1   ][1].v[1]-*cells[1][1   ][1].v[0])/dz)+
								((*cells[1][1][1   ].w[1] + *cells[1][1][1   ].w[0] - *cells[1][1-ny][1].w[1]-*cells[1][1-ny][1].w[0])/dy));
			}
			if(nx!=0 && nz!=0){
				getcell(&cells[1-nx][1][1],U,V,W,P,i-nx,j,k,imax,jmax,kmax,&zero);
				getcell(&cells[1][1][1-nz],U,V,W,P,i,j,k-nz,imax,jmax,kmax,&zero);
				//if(i== 1 && j == 12 && k == 2){printf("nx!=0 && nz!=0");};
				np = np-(1.0/(2.0*Re))*(
						((*cells[1][1][1-nz].u[1] + *cells[1][1][1-nz].u[0] - *cells[1   ][1][1].u[1] - *cells[1   ][1][1].u[0])/dz) +
						((*cells[1][1][1   ].w[1] + *cells[1][1][1   ].w[0] - *cells[1-nx][1][1].w[1] - *cells[1-nx][1][1].w[0])/dx));
			}
			if(nx!=0 && ny!=0){
				getcell(&cells[1-nx][1][1],U,V,W,P,i-nx,j,k,imax,jmax,kmax,&zero);
				getcell(&cells[1][1-ny][1],U,V,W,P,i,j-ny,k,imax,jmax,kmax,&zero);
				//if(i== 1 && j == 12 && k == 2){printf("nx!=0 && ny!=0");};
				np = np-(1.0/(2.0*Re))*(
						((*cells[1][1-ny][1].u[1] + *cells[1][1-ny][1].u[0] - *cells[1   ][1][1].u[1] -*cells[1   ][1][1].u[0])/dy) +
						((*cells[1][1   ][1].v[1] + *cells[1][1   ][1].v[0] - *cells[1-nx][1][1].v[1] -*cells[1-nx][1][1].v[0])/dx));
			}
			//if(i== 1 && j == 12 && k == 2){printf("*c->p = np;");};
			*c->p = np;

		}
	}
	//	if(i== 1 && j == 12 && k == 2){printf("done!");}
}


void set_uvwp_surface(double ***U,double ***V,double ***W,double ***P,int ***Flag,double dx,double dy,double dz, int imax,int jmax, int kmax,double GX,double GY,double GZ,double dt,double Re){

	printf("clean_empty_space\n");
	clean_empty_space(Flag,U,V,W,P,imax,jmax, kmax);

	int i,j,k;
	//int x=0,y=0,z=0;


	printf("compute_surface_values_1 (free surface cells)\n");

	for (i=1;i<=imax;i++){
		for (j=1;j<=jmax;j++){
			for (k=1;k<=kmax;k++){
				if (issurface(Flag[i][j][k])){
					compute_surface_values_1(U,V,W,P,Flag,i, j, k,imax,jmax,kmax,dx,dy,dz,dt,GX,GY,GZ,Re);


				}
			}
		}
	}

	printf("compute_surface_values_2 (cells next to free surface, planes )\n");
	for (i=1;i<=imax;i++){
		for (j=1;j<=jmax;j++){
			for (k=1;k<=kmax;k++){
				if (issurface(Flag[i][j][k])){
					compute_surface_values_2(U,V,W,P,Flag,i, j, k,imax,jmax,kmax,dx,dy,dz,dt,GX,GY,GZ,Re);
				}
			}
		}
	}
	printf("compute_surface_values_3 (cells next to free surface, edges and corners )\n");
	for (i=1;i<=imax;i++){
		for (j=1;j<=jmax;j++){
			for (k=1;k<=kmax;k++){
				if (issurface(Flag[i][j][k])){
					compute_surface_values_3(U,V,W,P,Flag,i, j, k,imax,jmax,kmax,dx,dy,dz,dt,GX,GY,GZ,Re);
				}
			}
		}
	}

	printf("compute_surface_pressure\n");
#pragma omp parallel for private(j,k)
	for (i=1;i<=imax;i++){
		for (j=1;j<=jmax;j++){
			for (k=1;k<=kmax;k++){
				if (issurface(Flag[i][j][k])){
					compute_surface_pressure(U,V,W,P,Flag,i, j, k,imax,jmax,kmax,dx,dy,dz,dt,GX,GY,GZ,Re);
				}
			}
		}
	}
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

	//struct particle *t;


	for(s=0;s<N;s++){
		for(p=Partlines[s].Particles; p->next != NULL; p=p->next){
			x = p->next->x;
			i1 = (int)(x/dx)+1;
			i2 = (int)((x+0.5*dx)/dx)+1;

			y = p->next->y;
			j1 = (int)(y/dy)+1;
			j2 = (int)((y+0.5*dy)/dy)+1;

			z = p->next->z;
			k1 = (int)(z/dz)+1;
			k2 = (int)((z+0.5*dz)/dz)+1;

			if(i1<1)
				i1 = 1;
			if(i2>imax)
				i2 = imax;

			if(j1<1)
				j1 = 1;
			if(j2>jmax)
				j2 = jmax ;

			if(k1<1)
				k1 = 1;
			if(k2>kmax)
				k2 = kmax ;

			if (i1<0 || j1<0 || k1<0 || i2>imax+1 || j2>jmax+1 || k2>kmax+1 ){
				/*
				t = p->next;
				p->next = p->next->next;
				free(t);
				Partlines[s].length--;
				if (p->next==NULL)
					break;
				 */
			}
			else
			{

#pragma omp parallel sections private(i, j, k,s,x1_1,x1_2,x2_1,x2_2, y1_1,y1_2,y2_1,y2_2,z1_1,z1_2,z2_1,z2_2)
				{
#pragma omp section
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


						x += dt*u;
						i = (int)(x/dx);
						p->next->x = x;

					}
#pragma omp section
					{
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

						y += dt*v;
						j = (int)(y/dy);
						p->next->y = y;

					}

#pragma omp section
					{
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
						z += dt*w;
						k = (int)(z/dz);
						p->next->z = z;
					}


				}


				/*
				if (i<0 || j<0 || k<0 || i>imax || j>jmax || k>kmax ){

					t = p->next;
					p->next = p->next->next;
					free(t);
					Partlines[s].length--;
					if (p->next==NULL)
						break;

				}
				else{
				 */
				/* todo boundary */


				p->next->vel = sqrt(u*u+v*v+w*w);
				/*
				}
				 */
			}
		}
	}
}






void get_particle_speed(double dx,double dy, double dz,int imax,int jmax,int kmax, double dt,double ***U,double ***V,double ***W,int N,struct particleline *Partlines){

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
	//struct particle *t;


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

			if (!(i1<=0 || j1<=0 || k1<=0 || i2>imax || j2>jmax || k2>kmax)){
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

				p->next->vel = sqrt(u*u+v*v+w*w);

			}
		}
	}
}





