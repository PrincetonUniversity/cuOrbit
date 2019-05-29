#ifndef SET_ORBIT_UTIL_H_
#define SET_ORBIT_UTIL_H_

#include <stdlib.h>

double darray_max(double*, size_t);
double darray_min(double*, size_t);

inline int imax(int a, int b){
  return a > b ? a : b;
};
inline int imin(int a, int b){
  return a < b ? a : b;
};

#endif
