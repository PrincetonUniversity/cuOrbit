/*  Copyright 2019, 2020 Garrett Wright, Princeton Plasma Physic Lab

    This file is part of CuOrbit.

    CuOrbit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CuOrbit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CuOrbit.  If not, see <https://www.gnu.org/licenses/>.
*/
                                  
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
