#ifndef SET_ORBIT_EQUILIBRIUM_H_
#define SET_ORBIT_EQUILIBRIUM_H_

#include "orbit_config.h"

#ifdef __NVCC__
#include <cuda.h>
#include <cuda_runtime.h>
#endif

/* this is generally speaking the values taken in from "spdata" */
typedef struct Equilib Equilib_t;

void initialize_Equilib(Equilib_t*, Config_t*);
Equilib_t* Equilib_ctor();

#ifdef __NVCC__
__host__ __device__
#endif
double** get_B(Equilib_t*);
#ifdef __NVCC__
__host__ __device__
#endif
double** get_G(Equilib_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double** get_R(Equilib_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double** get_X(Equilib_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double** get_Z(Equilib_t*);


#ifdef __NVCC__
__host__ __device__
#endif
double** get_QD(Equilib_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double** get_RD(Equilib_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double** get_GD(Equilib_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double** get_PD(Equilib_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double** get_PS(Equilib_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double** get_RP(Equilib_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double** get_VD(Equilib_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double get_pw(Equilib_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double get_rmaj(Equilib_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double get_ped(Equilib_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double get_lsp(Equilib_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double get_lst(Equilib_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double get_nrip(Equilib_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double get_krip(Equilib_t*);


/* utils */
#ifdef __NVCC__
__host__ __device__
#endif
int compute_jd(Equilib_t*, double);

double gfun(Equilib_t*, double);
double qfun(Equilib_t*, double);
double rifun(Equilib_t*, double);
double bfield(Equilib_t*, double, double);
double giac(Equilib_t*, double, double);
void vspline(Equilib_t*);

double xproj(Equilib_t*, double, double);
double zproj(Equilib_t*, double, double);

#ifdef __NVCC__
__host__ __device__
#endif
double rpol(Equilib_t*, double);

double polr_mp(Equilib_t*, double, double);
double compute_eps(Equilib_t*);

#endif
