#ifndef SET_ORBIT_UTIL_H_
#define SET_ORBIT_UTIL_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

double darray_max(double*, size_t);
double darray_min(double*, size_t);

inline int imax(int a, int b){
  return a > b ? a : b;
}
inline int imin(int a, int b){
  return a < b ? a : b;
}

inline bool atob(const char* a){
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
#endif
