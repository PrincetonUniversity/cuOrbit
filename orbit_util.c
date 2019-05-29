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
