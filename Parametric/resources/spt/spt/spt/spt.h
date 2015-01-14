//
//  spt.h
//  spt
//
//  Created by rchailan on 24/11/14.
//
//

#ifndef __spt__spt__
#define __spt__spt__

#include <stdio.h>

#endif /* defined(__spt__spt__) */


// Test function for openmp implementation of pwl
void testpwlOpenmp(double *x, double *y, double *nsites, double *ntimes, double *val);

// Parallel implementation of testpwlOpenmp(...)
void partestpwlOpenmp(double *x, double *y, double *nsites, double *ntimes, double *val);