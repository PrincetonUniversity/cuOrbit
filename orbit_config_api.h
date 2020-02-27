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

#ifndef SET_ORBIT_CONFIG_API_H_
#define SET_ORBIT_CONFIG_API_H_

#include <stdbool.h>
#include "orbit_config.h"
#include "orbit_perturbation.h"
#include "orbit_equilibrium.h"
#include "orbit_particles.h"
#include "orbit_deposition.h"

#ifdef __NVCC__
#include <cuda.h>
#include <cuda_runtime.h>
#endif

/* this is set in the .c file for now */
/* xxx, get rid of this static IDP business */
extern const int IDP;

/* alias to reduce likelyhood of name collision */
typedef Config_t orbit_Config_t;

struct Config {
  /* meta */
  char* name;

  /* model */
  int npert;  /* selects mode of pertubation codes */
  int nmds;  /* probably dont need this */
  int nrmds;  /* dont remember this */
  int seed;   /* used by the RNG */
  int nstep_all;
  double bsum; //xxx i suspect this will go to another struct...
  double dsum; //xxx i suspect this will go to another struct...
  double esum; //xxx i suspect this will go to another struct...
  double pvol;  /* stash plasma volume, another struct? */
  double engn;
  double bkg;
  double bmin;
  double bmax;
  double rmaj;
  double trun;
  double tran;  /* transit time for particle at the mag axi with pitch=1 */
  double dt0;
  double xc;
  double eps;
  double bax;
  double *dwal;  /* wall definition */
  double pamp;
  double rprof;

  /* inputs */
  char* spdata_file;
  char* alphas_file;
  char* displ_file;
  char* fbmdata_file;
  /* outputs */
  char* pdedp_file;
  bool output_sparse;
  char* pdedp_sparse_file;
  char* bfield_file;

  Equilib_t* eqlb_ptr;

  /* Pertubation */
  Perturb_t* ptrb_ptr;
  double falf;
  double ascale;
  double alimit;
  double global_scaling_factor;
  double freq_scaling_factor;
  double sng;
  bool do_modestep;

  /* particle */
  Particles_t* ptcl_ptr;
  int nprt;
  double zprt;
  double chrg;
  double prot;
  double ekev;
  double polo_scale;
  double p1_scale;
  double p2_scale;
  double pchi;

  /* Deposition*/
  Deposition_t* depo_ptr;
  /* xxx stochastic, does this belong here? */
  double mubk_scale;
  int emink;
  int emaxk;
  int dmubk;
  int nstoche;
  int nstochp;
  /* pdedp */
  int nruns;
  int pdedp_nbinE;
  int pdedp_nbinPz;
  int pdedp_nbinmu;
  int pdedp_nbinDE;
  int pdedp_nbinDPz;
  double pdedp_Emin;
  double pdedp_Emax;
  bool compute_pdedp;
  bool initial_update_pdedp_from_file;
  double deposit_on_bins_after_fraction;
  double pdedp_dtrun;
  double pdedp_dtsamp;
  double pdedp_dtav;
  int pdedp_tskip;
  double pdedp_otpup;
  bool pdedp_focusdep;
  bool pdedp_optimize;

};


orbit_Config_t* orbit_Config_ctor();
void orbit_initialize_Config(orbit_Config_t*, char*);

#ifdef __NVCC__
__host__ __device__
#endif
void set_xc(Config_t*, double);

#ifdef __NVCC__
__host__ __device__
#endif
double get_xc(Config_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double get_bkg(Config_t*);

#ifdef __NVCC__
__host__ __device__
#endif
void set_eps(Config_t*, double);

#ifdef __NVCC__
__host__ __device__
#endif

#ifdef __NVCC__
__host__ __device__
#endif
double get_eps(Config_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double get_engn(Config_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double get_bax(Config_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double get_bmax(Config_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double get_bmin(Config_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double get_pamp(Config_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double get_rprof(Config_t*);

#ifdef __NVCC__
__host__ __device__
#endif
int get_nstep_all(Config_t*);

#ifdef __NVCC__
__host__ __device__
#endif
double get_dt0(Config_t*);

void orbit_main_loop(orbit_Config_t* cfg_ptr);

#ifdef __NVCC__
__host__ __device__
#endif
int get_npert(Config_t*);

#endif
