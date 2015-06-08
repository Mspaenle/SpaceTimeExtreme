#include <stdio.h>
#include <stdlib.h> /* for malloc & free */
#include <math.h> 
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>	/* NA handling */
#include <Rmath.h>	
#include <R_ext/Random.h>	/* ..RNGstate */
#include <R_ext/Applic.h>	/* NA handling */


double L(double *x, double *n, double *theta) {
	double val, coef;
	int i;
	
	val = 0;
	coef =  (1/theta[2]) + 1 ;
	for (i = 0; i < (*n); i++) { 
		if (x[i] > (theta[0]-theta[1]/theta[2]) ) {
			val += log(1 + theta[2] * ((x[i]-theta[0])/theta[1]));
    	}
    }
	return (coef*val);
}

void spLikelihood(double *theta, double *x, double *n, double *val) {
	// *val = ( -(*n) * log(theta[1]) ) - L(x,n,theta); // Maximizing the log likelihood
	*val = ( (*n) * log(theta[1]) ) + L(x,n,theta); // Minimizing the minus log likelihood
	printf("likelihood: %f,mu: %f , sigma: %f, xi: %f \n",*val,theta[0],theta[1],theta[2]);
}

void printTheta(double *theta) {
	printf("mu: %f , sigma: %f, xi: %f\n",theta[0],theta[1],theta[2]);
}
