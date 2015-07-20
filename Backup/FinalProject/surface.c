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

void clean_cells(int ***Flag,int imax,int jmax, int kmax){
	int i,j,k;
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


void mark_cells(int ***Flag,double dx,double dy, double dz,int imax,int jmax,int kmax,int N,struct particleline *Partlines){
	int i,j,k,n;
	struct particle *p;//*t,
	double x,y,z;
	int *flag=0,*flagc=0;
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
			if (i<0 || j<0 || k<0 || i>imax || j>jmax || k>kmax){
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
	printf("markcell!\n");
	for (i=1;i<=imax;i++)
	{
		for (j=1;j<=jmax;j++){
			for (k=1;k<=kmax;k++)
			{
				flagc = &Flag[i][j  ][k  ];
				if( isfluid(*flagc) )
				{
					flag = &Flag[i+1][j  ][k  ];
					if(isempty(*flag)){  //B_O
						//Flag[i][j  ][k  ] = 0;
						//	printf("before: %d , %d\n ",Flag[i][j][k],issurface(Flag[i][j][k]));
						changebit(flagc,0,0);
						changebit(flagc,1,1);
						//	printf("is surface: %d, %d\n",Flag[i][j][k],issurface(Flag[i][j][k]));
						//printf("flag: %d\n",Flag[i][j  ][k  ]);
						//printf("is surface: %d\n",issurface(Flag[i][j][k]));
					}
					else if (isfluid(*flag)){
						changebit(flagc,0,1);
						changebit(flagc,1,0);
					}
					flag = &Flag[i-1][j  ][k  ];
					if(isempty(*flag)){
						//printf("markcell!\n");
						changebit(flagc,2,0); //B_W
						changebit(flagc,3,1);
					}
					else if (isfluid(*flag)){
						changebit(flagc,2,1);
						changebit(flagc,3,0);
					}

					flag = &Flag[i][j+1][k  ];
					if(isempty(*flag)){
						//printf("markcell!\n");
						changebit(flagc,4,0); //B_N
						changebit(flagc,5,1);
					}
					else if(isfluid(*flag)){
						changebit(flagc,4,1);
						changebit(flagc,5,0);
					}

					flag = &Flag[i][j-1][k  ];
					if(isempty(*flag)){
						//printf("markcell!\n");
						changebit(flagc,6,0); //B_S
						changebit(flagc,7,1);
					}
					else if(isfluid(*flag)){
						changebit(flagc,6,1);
						changebit(flagc,7,0);
					}

					flag = &Flag[i][j  ][k-1];
					if(isempty(*flag)){
						//printf("markcell!\n");
						changebit(flagc,8,0); //B_D
						changebit(flagc,9,1);
					}
					else if (isfluid(*flag)){
						changebit(flagc,8,1);
						changebit(flagc,9,0);
					}

					flag = &Flag[i][j  ][k+1];
					if(isempty(*flag)){
						//printf("markcell!\n");
						changebit(flagc,10,0);//B_U
						changebit(flagc,11,1);
					}
					else if (isfluid(*flag)){
						changebit(flag,10,1);
						changebit(flag,11,0);
					}
				}
			}
		}
	}
}

void clean_empty_space(int*** Flag,double ***U,double ***V,double ***W,double ***P,int imax,int jmax, int kmax){
	int i,j,k;
	for (i=0;i<=imax-1;i++){
		for (j=0;j<=jmax;j++){
			for (k=0;k<=kmax;k++){
				if (!isfluid(Flag[i][j][k]) && !isfluid(Flag[i+1][j][k])){
					U[i][j][k] = 0;
				}
			}
		}
	}
	for (i=0;i<=imax-1;i++){
		for (j=0;j<=jmax-1;j++){
			for (k=0;k<=kmax;k++){
				if (!isfluid(Flag[i][j][k]) && !isfluid(Flag[i][j+1][k])){
					V[i][j][k] = 0;
				}
			}
		}
	}
	for (i=0;i<=imax;i++){
		for (j=0;j<=jmax;j++){
			for (k=0;k<=kmax-1;k++){
				if (!isfluid(Flag[i][j][k]) && !isfluid(Flag[i][j][k+1])){
					W[i][j][k] = 0;
				}
			}
		}
	}
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




struct cell{
	double *p;
	double *ue;
	double *uw;
	double *vn;
	double *wt;
	double *wb;
};

void getcell(struct cell* c,double ***U,double ***V,double ***W,double ***P,int i, int j, int k){
	*c->p = P[i][j][k];
}


void boundaryvalues_surface(
			int i,
			int j,
			int k,
			double ***U,
			double ***V,
			double ***W,
			int ***Flag,
			int imax,int jmax,int kmax
			){
        // No-slip boundary conditions for U, V and W.
	// 26 (6 + 12 + 8) cases in total.

	int nx,ny,nz;
	int mx,my,mz;
	int num;

	int type =  getsurfacetype(Flag[i][j][k],&nx,&ny,&nz,&mx,&my,&mz,&num);

	for(int s=0;s<=1;s++){
		if(binMatch(type,B_O) ){
			i++;
		}
		if(binMatch(type,B_W) ){
			i--;
		}
		if(binMatch(type,B_N)){
			j++;
		}
		if(binMatch(type,B_S)){
			j--;
		}
		if(binMatch(type,B_U)){
			k++ ;
		}
		if(binMatch(type,B_D)){
			k-- ;
		}
		if(i>0&&j>0&&k>0&&i<imax&&j<jmax&&k<kmax){
			switch(type){
				case B_W:

					U[i  ][j  ][k  ]   = U[i+1][j  ][k  ];

					V[i  ][j-1][k  ] = V[i+1][j-1][k  ];
					V[i  ][j  ][k  ]   = V[i+1][j  ][k  ];

					W[i  ][j  ][k-1] = W[i+1][j  ][k-1];
					W[i  ][j  ][k  ]   = W[i+1][j  ][k  ];
					break;
				case B_O: //step outflow
					U[i-1][j  ][k  ] = U[i-2][j  ][k  ];

					V[i  ][j-1][k  ] = V[i-1][j-1][k  ];
					V[i  ][j  ][k  ] = V[i-1][j  ][k  ];

					W[i  ][j  ][k-1] = W[i-1][j  ][k-1];
					W[i  ][j  ][k  ] = W[i-1][j  ][k  ];
					break;
				case B_S:
					printf("S\n");
					V[i  ][j  ][k  ] = V[i  ][j+1][k  ];

					U[i-1][j  ][k  ] = U[i-1][j+1][k  ];
					U[i  ][j  ][k  ] = U[i  ][j+1][k  ];

					W[i  ][j  ][k-1] = W[i  ][j+1][k-1];
					W[i  ][j  ][k  ] = W[i  ][j+1][k  ];
					break;
				case B_N:
					V[i  ][j-1][k  ] = V[i  ][j-2][k  ];
					U[i-1][j  ][k  ] = U[i-1][j-1][k  ];
					U[i  ][j  ][k  ] = U[i  ][j-1][k  ];
					W[i  ][j  ][k-1] = W[i  ][j-1][k-1];
					W[i  ][j  ][k  ] = W[i  ][j-1][k  ];

					break;
				case B_D:
					W[i  ][j  ][k  ] = W[i  ][j  ][k+1];
					U[i-1][j  ][k  ] = U[i-1][j  ][k+1];
					U[i  ][j  ][k  ] = U[i  ][j  ][k+1];
					V[i  ][j-1][k  ] = V[i  ][j-1][k+1];
					V[i  ][j  ][k  ] = V[i  ][j  ][k+1];
					break;
				case B_U:
					W[i  ][j  ][k-1] = W[i  ][j  ][k-2];
					U[i-1][j  ][k  ] = U[i-1][j  ][k-1];
					U[i  ][j  ][k  ] = U[i  ][j  ][k-1];
					V[i  ][j-1][k  ] = V[i  ][j-1][k-1];
					V[i  ][j  ][k  ] = V[i  ][j  ][k-1];
					break;
				case B_SW:
					U[i  ][j  ][k  ] = U[i+1][j  ][k  ];
					U[i-1][j  ][k  ] = U[i-1][j+1][k  ];
					V[i  ][j  ][k  ] = V[i  ][j+1][k  ];
					V[i  ][j-1][k  ] = V[i+1][j-1][k  ];
					W[i  ][j  ][k  ] = (W[i  ][j+1][k  ] + W[i+1][j  ][k  ]) * 0.5;
					W[i  ][j  ][k-1] = (W[i  ][j+1][k-1] + W[i+1][j  ][k-1]) * 0.5;
					break;
				case B_SO:
					U[i-1][j  ][k  ] = U[i-2][j  ][k  ];
					U[i  ][j  ][k  ] = U[i  ][j+1][k  ];
					V[i  ][j  ][k  ] = V[i  ][j+1][k  ];
					V[i  ][j-1][k  ] = V[i-1][j-1][k  ];
					W[i  ][j  ][k  ] = (W[i  ][j+1][k  ] + W[i-1][j  ][k  ]) * 0.5;
					W[i  ][j  ][k-1] = (W[i  ][j+1][k-1] + W[i-1][j  ][k-1]) * 0.5;
					break;
				case B_SD:
					V[i  ][j  ][k  ] = V[i  ][j+1][k  ];
					V[i  ][j-1][k  ] = V[i  ][j-1][k+1];
					W[i  ][j  ][k  ] = W[i  ][j  ][k+1];
					W[i  ][j  ][k-1] = W[i  ][j+1][k-1];
					U[i  ][j  ][k  ] = (U[i  ][j+1][k  ] + U[i  ][j  ][k+1]) * 0.5;
					U[i-1][j  ][k  ] = (U[i-1][j+1][k  ] + U[i-1][j  ][k+1]) * 0.5;
					break;
				case B_SU:
					V[i  ][j  ][k  ] = V[i  ][j+1][k  ];
					V[i  ][j-1][k  ] = V[i  ][j-1][k-1];
					W[i  ][j  ][k-1] = W[i  ][j  ][k-2];
					W[i  ][j  ][k  ] = W[i  ][j+1][k  ];
					U[i  ][j  ][k  ] = (U[i  ][j+1][k  ] + U[i  ][j  ][k-1]) * 0.5;
					U[i-1][j  ][k  ] = (U[i-1][j+1][k  ] + U[i-1][j  ][k-1]) * 0.5;
					break;
				case B_NW:
					U[i  ][j  ][k  ] = U[i+1][j  ][k  ];
					U[i-1][j  ][k  ] = U[i-1][j-1][k  ];
					V[i  ][j-1][k  ] = V[i  ][j-2][k  ];
					V[i  ][j  ][k  ] = V[i+1][j  ][k  ];
					W[i  ][j  ][k  ] = (W[i  ][j-1][k  ] + W[i+1][j  ][k  ]) * 0.5;
					W[i  ][j  ][k-1] = (W[i  ][j-1][k-1] + W[i+1][j  ][k-1]) * 0.5;
					break;
				case B_NO:
					U[i-1][j  ][k  ] = U[i-2][j  ][k  ];
					U[i  ][j  ][k  ] = U[i  ][j-1][k  ];
					V[i  ][j-1][k  ] = V[i  ][j-2][k  ];
					V[i  ][j  ][k  ] = V[i-1][j  ][k  ];
					W[i  ][j  ][k  ] = (W[i  ][j-1][k  ] + W[i-1][j  ][k  ]) * 0.5;
					W[i  ][j  ][k-1] = (W[i  ][j-1][k-1] + W[i-1][j  ][k-1]) * 0.5;
					break;
				case B_ND:
					V[i  ][j-1][k  ] = V[i  ][j-2][k  ];
					V[i  ][j  ][k  ] = V[i  ][j  ][k+1];
					W[i  ][j  ][k  ] = W[i  ][j  ][k+1];
					W[i  ][j  ][k-1] = W[i  ][j-1][k-1];
					U[i  ][j  ][k  ] = (U[i  ][j-1][k  ] + U[i  ][j  ][k+1]) * 0.5;
					U[i-1][j  ][k  ] = (U[i-1][j-1][k  ] + U[i-1][j  ][k+1]) * 0.5;
					break;
				case B_NU:
					V[i  ][j-1][k  ] = V[i  ][j-2][k  ];
					V[i  ][j  ][k  ] = V[i  ][j  ][k-1];
					W[i  ][j  ][k-1] = W[i  ][j  ][k-2];
					W[i  ][j  ][k  ] = W[i  ][j-1][k  ];
					U[i  ][j  ][k  ] = (U[i  ][j-1][k  ] + U[i  ][j  ][k-1]) * 0.5;
					U[i-1][j  ][k  ] = (U[i-1][j-1][k  ] + U[i-1][j  ][k-1]) * 0.5;
					break;
				case B_WD:
					U[i  ][j  ][k  ] = U[i+1][j  ][k  ];
					U[i-1][j  ][k  ] = U[i-1][j  ][k+1];
					W[i  ][j  ][k  ] = W[i  ][j  ][k+1];
					W[i  ][j  ][k-1] = W[i+1][j  ][k-1];
					V[i  ][j  ][k  ] = (V[i+1][j  ][k  ] + V[i  ][j  ][k+1]) * 0.5;
					V[i  ][j-1][k  ] = (V[i+1][j-1][k  ] + V[i  ][j-1][k+1]) * 0.5;
					break;
				case B_OD:
					U[i-1][j  ][k  ] = U[i-2][j  ][k  ];
					U[i  ][j  ][k  ] = U[i  ][j  ][k+1];
					W[i  ][j  ][k  ] = W[i  ][j  ][k+1];
					W[i  ][j  ][k-1] = W[i-1][j  ][k-1];
					V[i  ][j  ][k  ] = (V[i-1][j  ][k  ] + V[i  ][j  ][k+1]) * 0.5;
					V[i  ][j-1][k  ] = (V[i-1][j-1][k  ] + V[i  ][j-1][k+1]) * 0.5;
					break;
				case B_WU:
					U[i  ][j  ][k  ] = U[i+1][j  ][k  ];
					U[i-1][j  ][k  ] = U[i-1][j  ][k-1];
					W[i  ][j  ][k-1] = W[i  ][j  ][k-2];
					W[i  ][j  ][k  ] = W[i+1][j  ][k  ];
					V[i  ][j  ][k  ] = (V[i+1][j  ][k  ] + V[i  ][j  ][k-1]) * 0.5;
					V[i  ][j-1][k  ] = (V[i+1][j-1][k  ] + V[i  ][j-1][k-1]) * 0.5;
					break;
				case B_OU:
					U[i-1][j  ][k  ] = U[i-2][j  ][k  ];
					U[i  ][j  ][k  ] = U[i  ][j  ][k-1];
					W[i  ][j  ][k-1] = W[i  ][j  ][k-2];
					W[i  ][j  ][k  ] = W[i-1][j  ][k  ];
					V[i  ][j  ][k  ] = (V[i-1][j  ][k  ] + V[i  ][j  ][k-1]) * 0.5;
					V[i  ][j-1][k  ] = (V[i-1][j-1][k  ] + V[i  ][j-1][k-1]) * 0.5;
					break;

				case B_SWD:
					U[i  ][j  ][k  ] = U[i+1][j  ][k  ];
					U[i-1][j  ][k  ] = (U[i-1][j+1][k  ] + U[i-1][j  ][k+1]) * 0.5;
					V[i  ][j  ][k  ] = V[i  ][j+1][k  ];
					V[i  ][j-1][k  ] = (V[i+1][j-1][k  ] + V[i  ][j-1][k+1]) * 0.5;
					W[i  ][j  ][k  ] = W[i  ][j  ][k+1];
					W[i  ][j  ][k-1] = (W[i+1][j  ][k-1] + W[i  ][j+1][k-1]) * 0.5;
					break;
				case B_SOD:
					U[i-1][j  ][k  ] = U[i-2][j  ][k  ];
					U[i  ][j  ][k  ] = (U[i  ][j+1][k  ] + U[i  ][j  ][k+1]) * 0.5;
					V[i  ][j  ][k  ] = V[i  ][j+1][k  ];
					V[i  ][j-1][k  ] = (V[i-1][j-1][k  ] + V[i  ][j-1][k+1]) * 0.5;
					W[i  ][j  ][k  ] = W[i  ][j  ][k+1];
					W[i  ][j  ][k-1] = (W[i-1][j  ][k-1] + W[i  ][j+1][k-1]) * 0.5;
					break;
				case B_SWU:
					U[i  ][j  ][k  ] = U[i+1][j  ][k  ];
					U[i-1][j  ][k  ] = (U[i-1][j+1][k  ] + U[i-1][j  ][k-1]) * 0.5;
					V[i  ][j  ][k  ] = V[i  ][j+1][k  ];
					V[i  ][j-1][k  ] = (V[i+1][j-1][k  ] + V[i  ][j-1][k-1]) * 0.5;
					W[i  ][j  ][k-1] = W[i  ][j  ][k-2];
					W[i  ][j  ][k  ] = (W[i+1][j  ][k  ] + W[i  ][j+1][k  ]) * 0.5;
					break;
				case B_SOU:
					U[i-1][j  ][k  ] = U[i-2][j  ][k  ];
					U[i  ][j  ][k  ] = (U[i  ][j+1][k  ] + U[i  ][j  ][k-1]) * 0.5;
					V[i  ][j  ][k  ] = V[i  ][j+1][k  ];
					V[i  ][j-1][k  ] = (V[i-1][j-1][k  ] + V[i  ][j-1][k-1]) * 0.5;
					W[i  ][j  ][k-1] = W[i  ][j  ][k-2];
					W[i  ][j  ][k  ] = (W[i-1][j  ][k  ] + W[i  ][j+1][k  ]) * 0.5;
					break;
				case B_NWD:
					U[i  ][j  ][k  ] = U[i+1][j  ][k  ];
					U[i-1][j  ][k  ] = (U[i-1][j-1][k  ] + U[i-1][j  ][k+1]) * 0.5;
					V[i  ][j-1][k  ] = V[i  ][j-2][k  ];
					V[i  ][j  ][k  ] = (V[i+1][j  ][k  ] + V[i  ][j  ][k+1]) * 0.5;
					W[i  ][j  ][k  ] = W[i  ][j  ][k+1];
					W[i  ][j  ][k-1] = (W[i+1][j  ][k-1] + W[i  ][j-1][k-1]) * 0.5;
					break;
				case B_NOD:
					U[i-1][j  ][k  ] = U[i-2][j  ][k  ];
					U[i  ][j  ][k  ] = (U[i  ][j-1][k  ] + U[i  ][j  ][k+1]) * 0.5;
					V[i  ][j-1][k  ] = V[i  ][j-2][k  ];
					V[i  ][j  ][k  ] = (V[i-1][j  ][k  ] + V[i  ][j  ][k+1]) * 0.5;
					W[i  ][j  ][k  ] = W[i  ][j  ][k+1];
					W[i  ][j  ][k-1] = (W[i-1][j  ][k-1] + W[i  ][j-1][k-1]) * 0.5;
					break;
				case B_NWU:
					U[i  ][j  ][k  ] = U[i+1][j  ][k  ];
					U[i-1][j  ][k  ] = (U[i-1][j-1][k  ] + U[i-1][j  ][k-1]) * 0.5;
					V[i  ][j-1][k  ] = V[i  ][j-2][k  ];
					V[i  ][j  ][k  ] = (V[i+1][j  ][k  ] + V[i  ][j  ][k-1]) * 0.5;
					W[i  ][j  ][k-1] = W[i  ][j  ][k-2];
					W[i  ][j  ][k  ] = (W[i+1][j  ][k  ] + W[i  ][j-1][k  ]) * 0.5;
					break;
				case B_NOU:
					U[i-1][j  ][k  ] = U[i-2][j  ][k  ];
					U[i  ][j  ][k  ] = (U[i  ][j-1][k  ] + U[i  ][j  ][k-1]) * 0.5;
					V[i  ][j-1][k  ] = V[i  ][j-2][k  ];
					V[i  ][j  ][k  ] = (V[i-1][j  ][k  ] + V[i  ][j  ][k-1]) * 0.5;
					W[i  ][j  ][k-1] = W[i  ][j  ][k-2];
					W[i  ][j  ][k  ] = (W[i-1][j  ][k  ] + W[i  ][j-1][k  ]) * 0.5;
					break;

				default: //case 0
					//printf("case 0\n");
					break;
			}
		}
	}

}



void compute_surface_values(double ***U,double ***V,double ***W,double ***P,int ***Flag,int i,int j,int k,double dt,double GX,double GY,double GZ,int imax, int jmax, int kmax){
	/*int type;*/
	struct cell c;
	int nx,ny,nz;
	int mx,my,mz;
	int num;
	getsurfacetype(Flag[i][j][k],&nx,&ny,&nz,&mx,&my,&mz,&num);
	getcell(&c,U,V,W,P,i,j,k);
	boundaryvalues_surface(i, j,k,U,V,W,Flag,imax,jmax,kmax);


}



void set_uvwp_surface(double ***U,double ***V,double ***W,double ***P,int ***Flag,double dx,double dy,double dz, int imax,int jmax, int kmax,double GX,double GY,double GZ,double dt,double Re){

	clean_empty_space(Flag,U,V,W,P,imax,jmax, kmax);

	int i,j,k;
	//int x=0,y=0,z=0;

	int type = 0;

	int nx,ny,nz;
	int mx,my,mz;
	int num;


	for (j=0;j<=jmax+1;j++){
		for (i=0;i<=imax+1;i++){
			for (k=0;k<=kmax+1;k++){

				if (issurface(Flag[i][j][k])){
					type = getsurfacetype(Flag[i][j][k],&nx,&ny,&nz,&mx,&my,&mz,&num);

					switch(type){

					case B_O   : P[i][j][k] =  P[i-1][j][k]; break;//2.0/Re/dx*(U[i][j][k]-U[i-1][j][k]); break;
					case B_W   : P[i][j][k] =  P[i][j+1][k]; break;//2.0/Re/dx*(U[i][j][k]-U[i-1][j][k]); break;
					case B_N   : P[i][j][k] =  P[i][j-1][k]; break;//2.0/Re/dy*(V[i][j][k]-V[i][j-1][k]); break;
					case B_S   : P[i][j][k] =  P[i][j+1][k]; break;//2.0/Re/dy*(V[i][j][k]-V[i][j-1][k]); break;

					case B_U   : P[i][j][k] =  P[i][j][k-1]; break;//2.0/Re/dz*(W[i][j][k]-W[i][j][k-1]); break;
					case B_D   : P[i][j][k] =  P[i][j][k+1]; break;//2.0/Re/dz*(W[i][j][k]-W[i][j][k-1]); break;

					default   :
						P[i][j][k] = 0; break;
						P[i][j][k] = 0; break;
						P[i][j][k] = 0; break;
						P[i][j][k] = 0; break;
						P[i][j][k] = 0; break;
						P[i][j][k] = 0; break;
						P[i][j][k] = 0; break;

						break;
					}

				}
			}
		}
	}

	for (j=2;j<jmax;j++){
		for (i=2;i<imax;i++){
			for (k=2;k<kmax;k++){

				if (issurface(Flag[i][j][k])){

					compute_surface_values(U,V,W,P,Flag,i, j, k,dt,GX,GY,GZ,imax,jmax,kmax);



				}
			}
		}
	}
	printf("%d",type);
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
			y = p->next->y;
			z = p->next->z;

			i1 = (int)(x/dx)+1;
			j1 = (int)(y/dy)+1;
			k1 = (int)(z/dz)+1;

			i2 = (int)((x+0.5*dx)/dx)+1;
			j2 = (int)((y+0.5*dy)/dy)+1;
			k2 = (int)((z+0.5*dz)/dz)+1;

			if (i1<=0 || j1<=0 || k1<=0 || i2>imax || j2>jmax || k2>kmax ){
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

				p->next->x = x;
				p->next->y = y;
				p->next->z = z;
				p->next->vel = sqrt(u*u+v*v+w*w);
				/*
				}
				 */
			}
		}
	}
}
