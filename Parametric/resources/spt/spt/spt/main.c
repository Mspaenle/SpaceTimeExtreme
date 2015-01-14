//
//  main.c
//  spt
//
//  Created by rchailan on 24/11/14.
//  Copyright (c) 2014 rchailan. All rights reserved.
//

#include <stdio.h>
#include "spt.h"
#include "timer.h"

int main(int argc, const char * argv[]) {
    printf("Hello!\n");
    
    double val;
    double x = 10;
    double y = 10;
    double nsites = 5000;
    double ntimes = 1000;
    
    
    // Sequential implementation
    double start = mysecond();
    testpwlOpenmp(&x, &y, &nsites, &ntimes, &val);
    double stop = mysecond();
    
    printf("Sequential mode:\t");
    printf("The value summed is: %lf ; ", val);
    printf("Job took %lf seconds\n", stop-start);
    
    // Parallel implementation
    start = mysecond();
    partestpwlOpenmp(&x, &y, &nsites, &ntimes, &val);
    stop = mysecond();
    
    printf("Parallel mode:\t");
    printf("The value summed is: %lf ; ", val);
    printf("Job took %lf seconds\n", stop-start);
    
    return 0;
}