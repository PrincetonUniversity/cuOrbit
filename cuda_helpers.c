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

void umafree(void** ptr){
#ifdef __NVCC__
  cudaFree(*ptr);
#else
  free(*ptr);
#endif
  *ptr=NULL;
}
