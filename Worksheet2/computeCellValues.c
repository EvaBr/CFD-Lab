#include "computeCellValues.h"
#include "LBDefinitions.h"

void computeDensity(const double *const currentCell, double *density){
	*density = 0;
	for (int i=0; i<Q; i++)
		*density += *(currentCell+i);
}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity){
	double sum_of_product[3] = { 0 }; // initialize all elements to 0

	// compute the sum of c_i*f_i for velocity
	for (int j=0; j<D; j++){
		for (int i=0; i<Q; i++){
			*(sum_of_product+j) += *(currentCell+i) * LATTICEVELOCITIES[i][j];
		}
	}
	// set velocity components
	for (int j=0; j<D; j++)
		*(velocity+j) = *(sum_of_product+j) / (*density);

}

void computeFeq(const double * const density, const double * const velocity, double *feq){
	double c_u_inner_product = 0;
        double u_u_inner_product = 0;

        for (int i=0; i<Q; i++){
		for (int j=0; j<D; j++){
			// compute inner product of c_i and u
			c_u_inner_product += *(velocity+j) * LATTICEVELOCITIES[i][j];
			// compute scalar product of u
			u_u_inner_product += *(velocity+j) * (*(velocity+j));
		}
		// compute equilibrium distribution
               	*(feq+i) = LATTICEWEIGHTS[i] * (*density) * (1 + c_u_inner_product/(C_S*C_S) +
		c_u_inner_product*c_u_inner_product/(2*C_S*C_S*C_S*C_S) - u_u_inner_product/(2*C_S*C_S));

		// initialize inner_prod variables back to 0 for further use
        	c_u_inner_product = 0;
        	u_u_inner_product = 0;
	}
}
