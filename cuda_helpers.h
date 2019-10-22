#ifndef _CUDA_HELPERS_H_
#define _CUDA_HELPERS_H_

#include <stdlib.h>
#include <stdio.h>

#ifdef __NVCC__
#define HANDLE_ERROR(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}
#endif


void* umacalloc(size_t num, size_t size);

#endif
