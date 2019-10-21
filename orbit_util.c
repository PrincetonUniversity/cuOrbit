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

int imax(int a, int b){
  return a > b ? a : b;
}
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
  return rand()/(double)RAND_MAX;
}

double rand_double_range(const double low, const double high){
    /* stupid simple c rng between low and high  */
   return (high - low) * rand_double() + low;
}
