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

#include <float.h>
#include "orbit_util.h"

double darray_max(double* A, size_t nA){
  double val = DBL_MIN;
  size_t n;
  for(n=0; n<nA; n++){
    if(A[n] > val){
      val =A[n];
    }
  }
  return val;
}

double darray_min(double* A, size_t nA){
  double val = DBL_MAX;
  size_t n;
  for(n=0; n<nA; n++){
    if(A[n] < val){
      val =A[n];
    }
  }
  return val;
}

#ifdef __NVCC__
__host__ __device__
#endif
int imax(int a, int b){
  return a > b ? a : b;
}

#ifdef __NVCC__
__host__ __device__
#endif
int imin(int a, int b){
  return a < b ? a : b;
}

bool atob(const char* a){
  if (strncmp(a, "true", 4) == 0){
    return true;
  } else if (strncmp(a, "false", 5) == 0){
    return false;
  } else {
    fprintf(stderr, "atob expects \"true\"/\"false\""
            "received invalid string: %s\n", a);
    exit(1);
  }
}

double rand_double(){
  /* stupid simple c rng between 0 1  */
  return rand()/((double)RAND_MAX + 1);
}

double rand_double_range(const double low, const double high){
    /* stupid simple c rng between low and high  */
   return (high - low) * rand_double() + low;
}
