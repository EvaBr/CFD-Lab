#include "computeCellValues.h"
#include "LBDefinitions.h"

void computeDensity(const double *const currentCell, double *density){
  /* TODO */
	*density = 0; //debug
	for (int i=0; i<Q; i++)
		*density += *(currentCell+i);
}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity){
  /* TODO */
	double sum_of_product[3] = { 0 }; // all elements 0
	for (int j=0; j<D; j++){
		for (int i=0; i<Q; i++){
			*(sum_of_product+j) += *(currentCell+i) * LATTICEVELOCITIES[i][j];
		}
	}

	for (int j=0; j<D; j++)
		*(velocity+j) = *(sum_of_product+j) / (*density);

}

void computeFeq(const double * const density, const double * const velocity, double *feq){
 /* TODO */

	double c_u_inner_product = 0;
        double u_u_inner_product = 0;
        for (int i=0; i<Q; i++){
		for (int j=0; j<D; j++){
			c_u_inner_product += *(velocity+j) * LATTICEVELOCITIES[i][j];
			u_u_inner_product += *(velocity+j) * (*(velocity+j));
		}
               	*(feq+i) = LATTICEWEIGHTS[i] * (*density) * (1 + c_u_inner_product/(C_S*C_S) +
		c_u_inner_product*c_u_inner_product/(2*C_S*C_S*C_S*C_S) - u_u_inner_product/(2*C_S*C_S));
        	c_u_inner_product = 0; //debug
        	u_u_inner_product = 0; //debug
	}
}
