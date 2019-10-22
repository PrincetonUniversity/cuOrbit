#ifndef SET_ORBIT_MAIN_H_
#define SET_ORBIT_MAIN_H_

#include "orbit_config_api.h"

#ifdef __NVCC__
#include <cuda.h>
#include <cuda_runtime.h>
#endif

void main_loop(Config_t* cfg_ptr);

#endif
