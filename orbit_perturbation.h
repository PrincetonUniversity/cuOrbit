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

#ifndef SET_ORBIT_PERTURBATION_H_
#define SET_ORBIT_PERTURBATION_H_

#include "orbit_config.h"
#include "orbit_particles.h"
#include "orbit_equilibrium.h"

#ifdef __NVCC__
#include <cuda.h>
#include <cuda_runtime.h>
#endif

/* we can pack this in the struct if we really need it... */
extern const int NAMP_;

typedef struct Perturb Perturb_t;

void initialize_Perturb(Perturb_t*, Config_t*, Equilib_t*, Particles_t*);
Perturb_t* Perturb_ctor();

void splna(Perturb_t* , Equilib_t*, Particles_t*);
void splnx(Perturb_t* , Equilib_t*, Particles_t*);
void set_omeg0(Perturb_t*, double);

#ifdef __NVCC__
__host__ __device__
#endif
double get_omeg0(Perturb_t*);

#ifdef __NVCC__
__host__ __device__
#endif
int get_nflr(Perturb_t*);

#ifdef __NVCC__
__host__ __device__
#endif
int get_lpt(Perturb_t*);

#ifdef __NVCC__
__host__ __device__
#endif
int get_md1(Perturb_t*);

#ifdef __NVCC__
__host__ __device__
#endif
int get_md2(Perturb_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_phaz(Perturb_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_omegv(Perturb_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_amp(Perturb_t*);

#ifdef __NVCC__
__host__ __device__
#endif
int* get_mmod(Perturb_t*);

#ifdef __NVCC__
__host__ __device__
#endif
int* get_nmod(Perturb_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_xi1(Perturb_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_xi2(Perturb_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_xi3(Perturb_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_a1(Perturb_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_a2(Perturb_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double* get_a3(Perturb_t*);

double pol2pot(Config_t*, double);

#ifdef __NVCC__
__host__ __device__
#endif
void modestep(Config_t* cfg_ptr);

#endif
