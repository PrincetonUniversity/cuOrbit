#include <stdlib.h>
#include <string.h>

#ifdef __NVCC__
#include <cuda.h>
#include <cuda_runtime.h>
#endif

#include "cuda_helpers.h"

/* Make an equiv to calloc for CUDA UMA*/
void* umacalloc(size_t num, size_t size){
  void* ptr;
#ifdef __NVCC__
  /* XXXX add cuda check on both of these returns later */
  HANDLE_ERROR(cudaMallocManaged(&ptr, num*size, cudaMemAttachGlobal));
  /* check if host memset faster, our arrays are generally small... */
  /* cudaMemset(ptr, 0, num*size); */
  memset(ptr, 0, num*size);
#else
  ptr = calloc(num, size);
#endif
  return ptr;
}
