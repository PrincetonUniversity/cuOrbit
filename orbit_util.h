#ifndef SET_ORBIT_UTIL_H_
#define SET_ORBIT_UTIL_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#ifdef __NVCC__
#include <cuda.h>
#include <cuda_runtime.h>
#endif


double darray_max(double*, size_t);

double darray_min(double*, size_t);

#ifdef __NVCC__
__host__ __device__
#endif
int imax(int a, int b);

#ifdef __NVCC__
__host__ __device__
#endif
int imin(int a, int b);

bool atob(const char* a);

double rand_double();

double rand_double_range(const double low, const double high);

#endif
