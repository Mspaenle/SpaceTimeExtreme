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
	coef = ( (1/theta[2]) + 1 );
	for (i = 0; i < (*n); i++) { 
		val += log(1 + theta[2] * ((x[i]-theta[0])/theta[1]));
	}
	
	return (coef*val);
}

void spLikelihood(double *x, double *n, double *theta, double *val) {
	*val = ( -(*n) * log(theta[1]) ) - L(x,n,theta);
}


