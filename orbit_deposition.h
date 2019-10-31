#ifndef SET_ORBIT_DEPOSITION_H_
#define SET_ORBIT_DEPOSITION_H_

#include "orbit_config.h"

#ifdef __NVCC__
#include <cuda.h>
#include <cuda_runtime.h>
#endif

typedef struct Deposition Deposition_t;
void initialize_Deposition(Deposition_t*, Config_t*);
Deposition_t* Deposition_ctor();
void initialize_pdedp(Deposition_t*);

#ifdef __NVCC__
__host__ __device__
#endif
bool compute_pdedp(Deposition_t*);

bool initial_update_pdedp(Deposition_t*);
bool pdedp_optimize(Deposition_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double get_pdedp_dtsamp(Deposition_t*);

#ifdef __NVCC__
__host__ __device__
#endif
int get_pdedp_tskip(Deposition_t*);

#ifdef __NVCC__
__host__ __device__
#endif
void set_pdedp_tskip(Deposition_t*, double);

#ifdef __NVCC__
__host__ __device__
#endif
bool get_initial_update_pdedp(Deposition_t*);

#ifdef __NVCC__
__host__ __device__
#endif
void set_initial_update_pdedp(Deposition_t*, bool);

#ifdef __NVCC__
__host__ __device__
#endif
bool get_pdedp_focusdep(Deposition_t*);

#ifdef __NVCC__
__host__ __device__
#endif
void set_pdedp_focusdep(Deposition_t*, bool);

#ifdef __NVCC__
__host__ __device__
#endif
size_t sizeof_pdedp(Deposition_t*);

void pdedp_read(Deposition_t*, Config_t*);
void pdedp_init(Deposition_t*);
void fulldepmp(Config_t*, Deposition_t*);
void fulldepmp_co(Config_t*, Deposition_t*);
void fullredepmp(Config_t* , Deposition_t*);
void class_kdomain(Config_t*, int);
void class_domain(Config_t*);
void pdedp_checkbdry(Config_t*, Deposition_t*);
void pdedp_finalize(Deposition_t*);
void pdedp_out(Deposition_t*);
void pdedp_rcrd_resid(Config_t*, Deposition_t*);
void rcrd_bfield(Config_t*, Deposition_t*);
void check_res_ptc(Config_t*, int);

#ifdef __NVCC__
__host__ __device__
#endif
void rcrd_vararr(Config_t* cfg_ptr, int k, int step);


#endif
