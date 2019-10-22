#ifndef SET_ORBIT_PARTICLES_H_
#define SET_ORBIT_PARTICLES_H_

#include <stdlib.h>

#ifdef __NVCC__
#include <cuda.h>
#include <cuda_runtime.h>
#endif

//not sure if this is better way yet
/* typedef struct Particle { */
/*   double chrg;  /\* 1 for ion, -1 for electron *\/ */
/*   double zprt; */
/*   double prot;  /\* proton mass in proton units *\/ */

/* } Particle_t; */

typedef struct Particles Particles_t;

void initialize_Particles(Particles_t*, Config_t*);
Particles_t* Particles_ctor();

#ifdef __NVCC__
__host__ __device__
#endif
double* get_b(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_g(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_q(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_en(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_pol(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_rho(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_rmu(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
int* get_otp(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_ptch(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_thet(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_pot(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_time(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_dt(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_tim1(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_wt(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_zet(Particles_t*);


#ifdef __NVCC__
__host__ __device__
#endif
double get_zprt(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double get_prot(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double get_prot(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double get_ekev(Particles_t*);


#ifdef __NVCC__
__host__ __device__
#endif
double* get_dadp(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_dadt(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_dadz(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_padt(Particles_t*);


void field(Config_t*, int);

#ifdef __NVCC__
__host__ __device__
#endif
void kfield(Config_t*, int);


#ifdef __NVCC__
__host__ __device__
#endif
void ptrbak(Config_t*, int);

#ifdef __NVCC__
__host__ __device__
#endif
void ptrb2k(Config_t*, int);

/* launcher... */
#ifdef __NVCC__
__global__
#endif
void do_particles(Config_t* );

#ifdef __NVCC__
__host__ __device__
#endif
void konestep(Config_t*, int);

#ifdef __NVCC__
__host__ __device__
#endif
void kupdate(Config_t* cfg_ptr, int k);

#endif
