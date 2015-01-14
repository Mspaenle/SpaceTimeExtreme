//
//  timer.c
//  spt
//
#include "timer.h"
#include <sys/time.h>

double mysecond() {
    struct timeval tp;
    struct timezone tzp;
    gettimeofday(&tp,&tzp);
    return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}


void mysecondR(double *top) {
    *top = mysecond();
}