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
double* get_ri(Particles_t* ptcl_ptr);

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

/* generic high level dispatch*/
void do_particles(Config_t*);

/* launches for host/device */
#ifdef __NVCC__
__global__
void do_particles_dev(Config_t*);
#endif
void do_particles_host(Config_t*);

/* per particle kernel */
#ifdef __NVCC__
__host__ __device__
#endif
void do_particle_kernel(Config_t*, int);

#ifdef __NVCC__
__host__ __device__
#endif
void konestep(Config_t*, int);

#ifdef __NVCC__
__host__ __device__
#endif
void kupdate(Config_t* cfg_ptr, int k);

#ifdef __NVCC__
__host__ __device__
#endif
int get_idm(Particles_t*);

#ifdef __NVCC__
__host__ __device__
#endif
void modestep(Config_t*);

void sampledep(Config_t*);

#endif
