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
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>

#include "orbit_config_api.h"
#include "orbit_deposition.h"
#include "orbit_equilibrium.h"  /* bfield */
#include "orbit_particles.h"  /* ekev */
#include "orbit_perturbation.h"
#include "orbit_util.h"
#include "cuda_helpers.h"

const size_t MAXFNAME=255;

struct Deposition {
  /* pdedp */
  char* pdedp_file;
  bool output_sparse;
  char* pdedp_sparse_file;
  char* bfield_file;
  int nruns;
  bool compute_pdedp;
  bool initial_update_pdedp_from_file;
  double deposit_on_bins_after_fraction;
  double pdedp_dtrun;
  double pdedp_dtsamp;
  double pdedp_dtav;
  int pdedp_tskip;
  double pdedp_otpup;
  bool pdedp_focusdep;
  bool pdedp_do_focusdep;
  bool pdedp_optimize;

  bool pdedp_initialized;
  /* the 5d dist */
  double* pdedp_pdedp;
  double* res_id_arr;
  /* internal vars */
  int pdedp_nbinDE;
  int pdedp_nbinDPz;
  int pdedp_nbinE;
  int pdedp_nbinPz;
  int pdedp_nbinmu;
  double* pdedp_varDE;
  double* pdedp_varDPz;
  double* pdedp_varE;
  double* pdedp_varPz;
  double* pdedp_varmu;
  double pdedp_Emin;
  double pdedp_Emax;
  double pdedp_Pzmin;
  double pdedp_Pzmax;
  double pdedp_mumin;
  double pdedp_mumax;
  double pdedp_DEmin;
  double pdedp_DEmax;
  double pdedp_DPzmin;
  double pdedp_DPzmax;
  double pdedp_maxDE;
  double pdedp_maxDPz;
  int res_id_arr_j;
  int res_id_arr_i;

  /* xxx stochastic, does this belong here? */
  double mubk;
  int emink;
  int emaxk;
  int dmubk;
  int nstoche;
  int nstochp;

  /* less annoying prints */
  bool recursing;

};

Deposition_t* Deposition_ctor(){
  return (Deposition_t*)umacalloc(1, sizeof(Deposition_t));
}

void initialize_Deposition(Deposition_t* Depo_ptr, Config_t* cfg_ptr){
  size_t sparse_filename_sz;
  Depo_ptr->pdedp_file = strndup(cfg_ptr->pdedp_file, MAXFNAME);
  Depo_ptr->output_sparse = cfg_ptr->output_sparse;
  if(Depo_ptr->output_sparse){
    if(! cfg_ptr->pdedp_sparse_file){
      printf("Error: output_sparse is `true` but pdedp_sparse_file is unset." \
             " Please set a pdedp_sparse_file name. Aborting\n");
      exit(1);
    } else {
      sparse_filename_sz = strnlen(cfg_ptr->pdedp_sparse_file, MAXFNAME);
      assert(sparse_filename_sz > 2);
      assert(sparse_filename_sz != MAXFNAME);  /* actually null term'd */
    }
    Depo_ptr->pdedp_sparse_file = strndup(cfg_ptr->pdedp_sparse_file, MAXFNAME);
  }
  Depo_ptr->bfield_file = strndup(cfg_ptr->bfield_file, MAXFNAME);
  Depo_ptr->nruns = cfg_ptr->nruns;
  Depo_ptr->compute_pdedp = cfg_ptr->compute_pdedp;
  Depo_ptr->initial_update_pdedp_from_file = cfg_ptr->initial_update_pdedp_from_file;
  Depo_ptr->deposit_on_bins_after_fraction = cfg_ptr->deposit_on_bins_after_fraction;
  Depo_ptr->pdedp_dtrun = cfg_ptr->pdedp_dtrun;
  Depo_ptr->pdedp_dtsamp = cfg_ptr->pdedp_dtsamp;
  Depo_ptr->pdedp_dtav = cfg_ptr->pdedp_dtav;
  /*Depo_ptr->pdedp_tskip = cfg_ptr->pdedp_tskip;*/
  Depo_ptr->pdedp_otpup = cfg_ptr->pdedp_otpup;
  Depo_ptr->pdedp_focusdep = cfg_ptr->pdedp_focusdep;
  Depo_ptr->pdedp_do_focusdep = false;
  Depo_ptr->pdedp_optimize = cfg_ptr->pdedp_optimize;
  if(Depo_ptr->initial_update_pdedp_from_file){
    printf("Warning, found initial pdedp_update TRUE, setting optimize FALSE\n");
    Depo_ptr->pdedp_optimize = 0;
    abort();  /* for now, abort */
  }
  
  /* Energy range and grid points */
  Depo_ptr->pdedp_Emin = cfg_ptr->pdedp_Emin;
  Depo_ptr->pdedp_Emax = cfg_ptr->pdedp_Emax;
  Depo_ptr->pdedp_nbinE = cfg_ptr->pdedp_nbinE;
  Depo_ptr->pdedp_nbinPz = cfg_ptr->pdedp_nbinPz;
  Depo_ptr->pdedp_nbinmu = cfg_ptr->pdedp_nbinmu;
  Depo_ptr->pdedp_nbinDE = cfg_ptr->pdedp_nbinDE;
  Depo_ptr->pdedp_nbinDPz = cfg_ptr->pdedp_nbinDPz;
  
  
  /* Allocate array with particles' variables vs time */
  
  /* Depo_ptr->res_id_arr_j = 10000;  xxx should be computed not sure we know value yet in current implientation */
  /* const int nstep_all = round(get_pdedp_dtrun(depo_ptr) / dum); */
  double dum = 1E3 * cfg_ptr->dt0 / get_omeg0(cfg_ptr->ptrb_ptr);
  const int nstep_all = round(get_pdedp_dtrun(Depo_ptr) / dum);
  cfg_ptr->nstep_all = nstep_all;
  set_pdedp_tskip(Depo_ptr,imax(ceil(nstep_all / 2E4),1));    /* stay in array bounds */
  cfg_ptr->pdedp_tskip = Depo_ptr->pdedp_tskip;
  Depo_ptr->res_id_arr_j = ceil(cfg_ptr->nstep_all / cfg_ptr->pdedp_tskip); /* MP: get size based on run time */
  Depo_ptr->res_id_arr_i = 4;
  Depo_ptr->res_id_arr = (double*)umacalloc((unsigned)(cfg_ptr->nprt *
                                                       Depo_ptr->res_id_arr_j *
                                                       Depo_ptr->res_id_arr_i),
                                            sizeof(double));
  
  
  /* this allocs at runtime, so has own init func */
  Depo_ptr->pdedp_initialized = false;

  /* xxx stochastic,  does this actually belong here*/
  Depo_ptr->mubk = cfg_ptr->mubk_scale *  get_ekev(cfg_ptr->ptcl_ptr);
  Depo_ptr->emink = cfg_ptr->emink;
  Depo_ptr->emaxk = cfg_ptr->emaxk;
  Depo_ptr->dmubk = cfg_ptr->dmubk;
  Depo_ptr->nstoche = cfg_ptr->nstoche;
  Depo_ptr->nstochp = cfg_ptr->nstochp;

  /* less annoying recursive printing */
  Depo_ptr->recursing = false;

  return;
}

void initialize_pdedp(Deposition_t* Depo_ptr){
  const size_t sz = sizeof(double);
  Depo_ptr->pdedp_varDE = (double*)umacalloc((unsigned)Depo_ptr->pdedp_nbinDE, sz);
  Depo_ptr->pdedp_varDPz = (double*)umacalloc((unsigned)Depo_ptr->pdedp_nbinDPz, sz);
  Depo_ptr->pdedp_varE = (double*)umacalloc((unsigned)Depo_ptr->pdedp_nbinE, sz);
  Depo_ptr->pdedp_varPz = (double*)umacalloc((unsigned)Depo_ptr->pdedp_nbinPz, sz);
  Depo_ptr->pdedp_varmu = (double*)umacalloc((unsigned)Depo_ptr->pdedp_nbinmu, sz);
  Depo_ptr->pdedp_pdedp = (double*)umacalloc(sizeof_pdedp(Depo_ptr), sz);
  Depo_ptr->pdedp_initialized = true;
}


#ifdef __NVCC__
__host__ __device__
#endif
bool compute_pdedp(Deposition_t* Depo_ptr){
  return Depo_ptr->compute_pdedp;
}

/* bool initial_update_pdedp_from_file(Deposition_t* Depo_ptr){ */
/*   return Depo_ptr->initial_update_pdedp_from_file; */
/* } */

bool pdedp_optimize(Deposition_t* Depo_ptr){
  return Depo_ptr->pdedp_optimize;
}

#ifdef __NVCC__
__host__ __device__
#endif
double get_pdedp_dtrun(Deposition_t* Depo_ptr){
  return Depo_ptr->pdedp_dtrun;
}

#ifdef __NVCC__
__host__ __device__
#endif
double get_pdedp_dtsamp(Deposition_t* Depo_ptr){
  return Depo_ptr->pdedp_dtsamp;
}

#ifdef __NVCC__
__host__ __device__
#endif
int get_pdedp_tskip(Deposition_t* Depo_ptr){
  return Depo_ptr->pdedp_tskip;
}

#ifdef __NVCC__
__host__ __device__
#endif
void set_pdedp_tskip(Deposition_t* Depo_ptr, int pdedp_tskip){
  Depo_ptr->pdedp_tskip = pdedp_tskip;
  return;
}

#ifdef __NVCC__
__host__ __device__
#endif
bool get_initial_update_pdedp_from_file(Deposition_t* Depo_ptr){
  return Depo_ptr->initial_update_pdedp_from_file;
}

#ifdef __NVCC__
__host__ __device__
#endif
void set_initial_update_pdedp_from_file(Deposition_t* Depo_ptr, bool val){
  Depo_ptr->initial_update_pdedp_from_file = val;
  return;
}

#ifdef __NVCC__
__host__ __device__
#endif
bool get_pdedp_focusdep(Deposition_t* Depo_ptr){
  return Depo_ptr->pdedp_focusdep;
}

#ifdef __NVCC__
__host__ __device__
#endif
void set_pdedp_focusdep(Deposition_t* Depo_ptr, bool val){
  Depo_ptr->pdedp_focusdep = val;
  return;
}

#ifdef __NVCC__
__host__ __device__
#endif
bool get_pdedp_do_focusdep(Deposition_t* Depo_ptr){
  return Depo_ptr->pdedp_do_focusdep;
}

#ifdef __NVCC__
__host__ __device__
#endif
void set_pdedp_do_focusdep(Deposition_t* Depo_ptr, bool val){
  Depo_ptr->pdedp_do_focusdep = val;
  return;
}

#ifdef __NVCC__
__host__ __device__
#endif
size_t sizeof_pdedp(Deposition_t* Depo_ptr){
  return (size_t) (Depo_ptr->pdedp_nbinE * Depo_ptr->pdedp_nbinPz * Depo_ptr->pdedp_nbinmu *
                   Depo_ptr->pdedp_nbinDE * Depo_ptr->pdedp_nbinDPz);
}

#ifdef __NVCC__
__host__ __device__
#endif
static inline int get_pdedp_ind(Deposition_t* Depo_ptr, int iE, int iPz, int imu, int iDE, int iDPz){
  const int ind = Depo_ptr->pdedp_nbinPz * Depo_ptr->pdedp_nbinmu * Depo_ptr->pdedp_nbinDE * Depo_ptr->pdedp_nbinDPz * iE +
      Depo_ptr->pdedp_nbinmu * Depo_ptr->pdedp_nbinDE * Depo_ptr->pdedp_nbinDPz * iPz +
      Depo_ptr->pdedp_nbinDE * Depo_ptr->pdedp_nbinDPz * imu +
      Depo_ptr->pdedp_nbinDPz * iDE + iDPz;
  return ind;
}

#ifdef __NVCC__
__host__ __device__
#endif
static inline int get_epmu_ind(Deposition_t* Depo_ptr, int iE, int iPz, int imu){
  const int ind = Depo_ptr->pdedp_nbinPz * Depo_ptr->pdedp_nbinmu * iE +
      Depo_ptr->pdedp_nbinmu * iPz + imu;
  return ind;
}

#ifdef __NVCC__
__host__ __device__
#endif
static inline int get_bin(double val, double* arr, const int dim){
  int k;
  int indx = -1;
  double epsF = 1E12;
  double darr = 0.5 * (arr[1] - arr[0]);
  double tmp;
  darr *= 1.0001;  /* something about round off */
  for(k=0; k<dim; k++){
    tmp = fabs(arr[k] - val);
    if(tmp <= epsF && tmp <= darr){
      indx = k;
      epsF = tmp;
    }
    if(epsF<darr) break;  /* min found */
  }
  return indx;
}

static void write_grid(Deposition_t* Depo_ptr, FILE* fh){
  int k;
  /* variables DEbins,DPbins */
  for(k=0; k < Depo_ptr->pdedp_nbinDE; k++){
    fprintf(fh, "%f ", Depo_ptr->pdedp_varDE[k]);
  }
  fprintf(fh,"\n");

  for(k=0; k < Depo_ptr->pdedp_nbinDPz; k++){
    fprintf(fh, "%f ", Depo_ptr->pdedp_varDPz[k]);
  }
  fprintf(fh,"\n");

  /* variables Ebins,Pbins,Mubins */
  for(k=0; k < Depo_ptr->pdedp_nbinE; k++){
    fprintf(fh, "%f ", Depo_ptr->pdedp_varE[k]);
  }
  fprintf(fh,"\n");

  for(k=0; k < Depo_ptr->pdedp_nbinPz; k++){
    fprintf(fh, "%f ", Depo_ptr->pdedp_varPz[k]);
  }
  fprintf(fh,"\n");

  for(k=0; k < Depo_ptr->pdedp_nbinmu; k++){
    fprintf(fh, "%f ", Depo_ptr->pdedp_varmu[k]);
  }
  fprintf(fh,"\n");
  return;
}

void pdedp_read(Deposition_t* Depo_ptr, Config_t* cfg_ptr){
  /* this reads the probability distribution data
     for P(DE,DP| E,P,mu) from a text file (UFILE)*/
  int i;
  int je, jp, jmu, jde, jdp;
  FILE *ifp;
  const char *mode = "r";

  class_domain(cfg_ptr);

  ifp = fopen(Depo_ptr->pdedp_file, mode);  /* xxx add f status */
  if (ifp == NULL) {
    fprintf(stderr, "\nCan't open input file %s!\n", Depo_ptr->pdedp_file);
    exit(1);
  }
  printf("\nParsing P(DE,DP) file %s\n",  Depo_ptr->pdedp_file);

  /* parse the file */
  for (i=0; i<9; i++){
    fscanf(ifp, "%*[^\n]\n");  /* skip line */
  }

  /* read in the time step */
  fscanf(ifp, "%lf %*[^\n]\n", &Depo_ptr->pdedp_dtsamp);
  Depo_ptr->pdedp_dtrun = Depo_ptr->pdedp_dtsamp*10.; /* default run duration when reading from pDEDP file */
  Depo_ptr->pdedp_dtav = Depo_ptr->pdedp_dtsamp/50.;

  /* continue skipping header */
  fscanf(ifp, "%*[^\n]\n");  /* skip line */
  fscanf(ifp, "%*[^\n]\n");  /* skip line */

  fscanf(ifp, "%d %*[^\n]\n", &Depo_ptr->pdedp_nbinDE);
  fscanf(ifp, "%d %*[^\n]\n", &Depo_ptr->pdedp_nbinDPz);
  fscanf(ifp, "%d %*[^\n]\n", &Depo_ptr->pdedp_nbinE);
  fscanf(ifp, "%d %*[^\n]\n", &Depo_ptr->pdedp_nbinPz);
  fscanf(ifp, "%d %*[^\n]\n", &Depo_ptr->pdedp_nbinmu);

  printf("\nNumber of points for p(DE,DP|E,P,mu) function:\n");
  printf("\t E\t P\t  mu\t DE\t DP\n");
  printf("\t%d\t%d\t%d\t%d\t%d\n",
         Depo_ptr->pdedp_nbinE, Depo_ptr->pdedp_nbinPz, Depo_ptr->pdedp_nbinmu,
         Depo_ptr->pdedp_nbinDE, Depo_ptr->pdedp_nbinDPz);

  /* malloc */
  if(! Depo_ptr->pdedp_initialized){
    initialize_pdedp(Depo_ptr);
  }

  for(i=0; i<Depo_ptr->pdedp_nbinDE; i++){
    fscanf(ifp, "%lf ", &Depo_ptr->pdedp_varDE[i]);
  }
  for(i=0; i<Depo_ptr->pdedp_nbinDPz; i++){
    fscanf(ifp, "%lf ", &Depo_ptr->pdedp_varDPz[i]);
  }
  for(i=0; i<Depo_ptr->pdedp_nbinE; i++){
    fscanf(ifp, "%lf ", &Depo_ptr->pdedp_varE[i]);
  }
  for(i=0; i<Depo_ptr->pdedp_nbinPz; i++){
    fscanf(ifp, "%lf ", &Depo_ptr->pdedp_varPz[i]);
  }
  for(i=0; i<Depo_ptr->pdedp_nbinmu; i++){
    fscanf(ifp, "%lf ", &Depo_ptr->pdedp_varmu[i]);
  }

  //     Define boundaries of computing grid.
  Depo_ptr->pdedp_Emin = Depo_ptr->pdedp_varE[0] - .5 * (
      Depo_ptr->pdedp_varE[Depo_ptr->pdedp_nbinE-1] - Depo_ptr->pdedp_varE[0]) /(
          Depo_ptr->pdedp_nbinE - 1.);

  Depo_ptr->pdedp_Emax = Depo_ptr->pdedp_varE[Depo_ptr->pdedp_nbinE - 1] + .5*(
      Depo_ptr->pdedp_varE[Depo_ptr->pdedp_nbinE] - Depo_ptr->pdedp_varE[0])/(
          Depo_ptr->pdedp_nbinE - 1.);

  Depo_ptr->pdedp_Pzmin = Depo_ptr->pdedp_varPz[0] - .5*(
      Depo_ptr->pdedp_varPz[Depo_ptr->pdedp_nbinPz-1] - Depo_ptr->pdedp_varPz[0])/(
          Depo_ptr->pdedp_nbinPz - 1.);

  Depo_ptr->pdedp_Pzmax = Depo_ptr->pdedp_varPz[Depo_ptr->pdedp_nbinPz-1] + .5*(
      Depo_ptr->pdedp_varPz[Depo_ptr->pdedp_nbinPz] - Depo_ptr->pdedp_varPz[0])/(
          Depo_ptr->pdedp_nbinPz - 1.);

  Depo_ptr->pdedp_mumin = Depo_ptr->pdedp_varmu[0] - .5*(
      Depo_ptr->pdedp_varmu[Depo_ptr->pdedp_nbinmu - 1] - Depo_ptr->pdedp_varmu[0])/(
          Depo_ptr->pdedp_nbinmu - 1.);

  Depo_ptr->pdedp_mumax = Depo_ptr->pdedp_varmu[Depo_ptr->pdedp_nbinmu-1] + .5*(
      Depo_ptr->pdedp_varmu[Depo_ptr->pdedp_nbinmu - 1] - Depo_ptr->pdedp_varmu[0])/(
          Depo_ptr->pdedp_nbinmu - 1.);

  Depo_ptr->pdedp_DEmin = Depo_ptr->pdedp_varDE[0] - .5*(
      Depo_ptr->pdedp_varDE[Depo_ptr->pdedp_nbinDE-1] - Depo_ptr->pdedp_varDE[0])/(
          Depo_ptr->pdedp_nbinDE - 1.);

  Depo_ptr->pdedp_DEmax = Depo_ptr->pdedp_varDE[Depo_ptr->pdedp_nbinDE - 1] + .5*(
      Depo_ptr->pdedp_varDE[Depo_ptr->pdedp_nbinDE - 1] - Depo_ptr->pdedp_varDE[0])/(
          Depo_ptr->pdedp_nbinDE - 1.);

  Depo_ptr->pdedp_DPzmin = Depo_ptr->pdedp_varDPz[0] - .5*(
      Depo_ptr->pdedp_varDPz[Depo_ptr->pdedp_nbinDPz - 1] - Depo_ptr->pdedp_varDPz[0]) / (
          Depo_ptr->pdedp_nbinDPz - 1. );

  Depo_ptr->pdedp_DPzmax = Depo_ptr->pdedp_varDPz[Depo_ptr->pdedp_nbinDPz - 1] + .5*(
      Depo_ptr->pdedp_varDPz[Depo_ptr->pdedp_nbinDPz - 1] -
      Depo_ptr->pdedp_varDPz[0]) / (Depo_ptr->pdedp_nbinDPz - 1.);

  /* read the 5d distribution in */
  /*  careful, the loop order is trecherous */
  for(je=0; je < Depo_ptr->pdedp_nbinE; je++){
    for(jp=0; jp < Depo_ptr->pdedp_nbinPz; jp++){
      for(jmu=0; jmu < Depo_ptr->pdedp_nbinmu; jmu++){
        for(jde=0; jde < Depo_ptr->pdedp_nbinDE; jde++){
          for(jdp=0; jdp < Depo_ptr->pdedp_nbinDPz; jdp++){
            i = get_pdedp_ind(Depo_ptr, je, jp, jmu, jde, jdp);
            fscanf(ifp, "%lf ", &Depo_ptr->pdedp_pdedp[i]);
          }
        }
      }
    }
  }

  fclose(ifp);

  return;
}

void pdedp_init(Deposition_t* Depo_ptr){
  int k;
  double stp;

  /* xxx I suspect we will want to migrate towards config */
  if(Depo_ptr->pdedp_Emin == 0){
     Depo_ptr->pdedp_Emin = 1.;
  }
  if(Depo_ptr->pdedp_Emax == 0){
     Depo_ptr->pdedp_Emax = 110.;
  }
  Depo_ptr->pdedp_Pzmin = -1.2;
  Depo_ptr->pdedp_Pzmax =  0.7;
  Depo_ptr->pdedp_mumin = 0.;  /* muB0/E */
  Depo_ptr->pdedp_mumax = 1.4;
  Depo_ptr->pdedp_DEmin = -0.01;
  Depo_ptr->pdedp_DEmax =  0.01;
  Depo_ptr->pdedp_DPzmin = -0.0001;
  Depo_ptr->pdedp_DPzmax =  0.0001;
  if(Depo_ptr->pdedp_nbinE == 0){
     Depo_ptr->pdedp_nbinE = 12;  /* number of bins */
  }
  if(Depo_ptr->pdedp_nbinPz == 0){
     Depo_ptr->pdedp_nbinPz = 40;
  }
  if(Depo_ptr->pdedp_nbinmu == 0){
     Depo_ptr->pdedp_nbinmu = 16;
  }
  if(Depo_ptr->pdedp_nbinDE == 0){
     Depo_ptr->pdedp_nbinDE = 29;   /* this must be an ODD number */
     assert(Depo_ptr->pdedp_nbinDE % 2 == 1);
  }
  if(Depo_ptr->pdedp_nbinDPz == 0){
     Depo_ptr->pdedp_nbinDPz = 29;  /* this must be an ODD number; */
     assert(Depo_ptr->pdedp_nbinDE % 2 == 1);
  }
  
  /* malloc */
  if(! Depo_ptr->pdedp_initialized){
    initialize_pdedp(Depo_ptr);
  }

  /* fill in grid vectors */
  /* Energy  */
  stp = (Depo_ptr->pdedp_Emax - Depo_ptr->pdedp_Emin) / Depo_ptr->pdedp_nbinE;
  for(k=0; k < Depo_ptr->pdedp_nbinE; k++){
    Depo_ptr->pdedp_varE[k] = Depo_ptr->pdedp_Emin + k * stp + stp/2.;
  }

  /* Pz */
  stp = (Depo_ptr->pdedp_Pzmax - Depo_ptr->pdedp_Pzmin) / Depo_ptr->pdedp_nbinPz;
  for(k=0; k < Depo_ptr->pdedp_nbinPz; k++){
    Depo_ptr->pdedp_varPz[k] = Depo_ptr->pdedp_Pzmin + k * stp + stp/2.;
  }

  /* mu Bo/E */
  stp = (Depo_ptr->pdedp_mumax - Depo_ptr->pdedp_mumin) / Depo_ptr->pdedp_nbinmu;
  for(int k=0; k < Depo_ptr->pdedp_nbinmu; k++){
    Depo_ptr->pdedp_varmu[k] = Depo_ptr->pdedp_mumin + k * stp + stp/2.;
  }

  /* Delta-Energy */
  stp = (Depo_ptr->pdedp_DEmax - Depo_ptr->pdedp_DEmin) / Depo_ptr->pdedp_nbinDE;
  for(k=0; k < Depo_ptr->pdedp_nbinDE; k++){
    Depo_ptr->pdedp_varDE[k] = Depo_ptr->pdedp_DEmin + k * stp + stp/2.;
  }

  /*Delta-Pz */
  stp = (Depo_ptr->pdedp_DPzmax - Depo_ptr->pdedp_DPzmin) / Depo_ptr->pdedp_nbinDPz;
  for(k=0; k < Depo_ptr->pdedp_nbinDPz; k++){
    Depo_ptr->pdedp_varDPz[k] = Depo_ptr->pdedp_DPzmin + k * stp + stp/2.;
  }

  /* reset to 0 , should have been alloc'd in initialize_pdedp */
  memset(Depo_ptr->pdedp_pdedp, 0, sizeof_pdedp(Depo_ptr)*sizeof(double));

  /*      -------------------------------------------------------
          Initialize variables to monitor the maximum
          kicks in energy and Pz. These are used to optimize
          the (DE,DPz) range on-the-fly if the flag
          pdedp_optimize is set true*/
  Depo_ptr->pdedp_maxDE = 0.;
  Depo_ptr->pdedp_maxDPz = 0.;

  return;
}

void class_kdomain(Config_t* cfg_ptr, int k){


  /*    -  Orbit classification
        otp=1, co-passing confined, otp = 2, co-passing lost
        otp=3, ctr-passing confined, otp = 4, ctr-passing lost
        otp=5, trapped confined, otp = 6, trapped lost
        otp=7, stagnation, otp = 8, conf potato, otp = 9, trap potato
  */

  double pdum,edum,dum;
  double pzdum,elt,ert,eax,etp1,etp2,mu;
  double rhod;
  double E_ax,E_min_lcfs,E_max_lcfs;
  double E_th0,E_thpi;
  const double ekev=get_ekev(cfg_ptr->ptcl_ptr);
  const double engn=get_engn(cfg_ptr);
  const double pw=get_pw(cfg_ptr->eqlb_ptr);
  const double *pol=get_pol(cfg_ptr->ptcl_ptr);
  const double *en=get_en(cfg_ptr->ptcl_ptr);
  const double *g=get_g(cfg_ptr->ptcl_ptr);
  const double *rho=get_rho(cfg_ptr->ptcl_ptr);
  const double *rmu=get_rmu(cfg_ptr->ptcl_ptr);
  int * const otp=get_otp(cfg_ptr->ptcl_ptr);  /* writes */
  const double *ptch=get_ptch(cfg_ptr->ptcl_ptr);
  const double bax = get_bax(cfg_ptr);
  const double bmax = get_bmax(cfg_ptr);
  const double bmin = get_bmin(cfg_ptr);

  mu = rmu[k]*ekev/engn;
  pdum = (g[k]*rho[k] - pol[k])/pw;     /* pz0[k] */
  edum = en[k]*ekev/engn;
  pzdum = pdum*pw;      /* pz0[k]*pw   Pz */
  rhod = (pzdum + pol[k])/g[k] ; /*  particle rho */
  dum = pow(pzdum+pw,2) * ekev / (
      engn * 2 * pow(gfun(cfg_ptr->eqlb_ptr, pw),2));

  /* get potential energies */
  E_ax = pol2pot(cfg_ptr, 0.);
  E_min_lcfs = pol2pot(cfg_ptr, pw);
  E_max_lcfs = pol2pot(cfg_ptr, -pzdum);
  E_th0 = pol2pot(cfg_ptr, -pzdum);
  E_thpi = pol2pot(cfg_ptr, -pzdum);

  ert =  dum*pow(bmin, 2) + mu*bmin + E_min_lcfs;
  elt =  dum*pow(bmax, 2) + mu*bmax + E_max_lcfs;
  dum = pow(rhod, 2) * ekev/(2*engn);
  eax = ekev*pow(pzdum,2) * pow(bax,2)/(
      engn * 2 * pow(gfun(cfg_ptr->eqlb_ptr, 0.),2)) + mu*bax + E_ax;
  etp1 = mu*bfield(cfg_ptr->eqlb_ptr, -pzdum, 0.) + E_th0;
  etp2 = mu*bfield(cfg_ptr->eqlb_ptr, -pzdum, M_PI) + E_thpi;

  /* Right */
  if(pdum > 0.){
    if(edum > eax) otp[k] = 1;
    if(edum < eax) otp[k] = 7;
    return;
  }
  /* Left */
  if(pdum < -1.){
    if(edum < elt) otp[k] = 4;
    if(edum >= elt) otp[k] = 3;
    return;
  }

  /* Middle */
  if(pdum > -1 && pdum < 0.){
    if(edum >= elt) otp[k] = 3;
    if(edum > eax && edum < etp2) otp[k] = 8;
    if(edum > eax && edum > ert) otp[k] = 2;
    if(edum < eax && edum < etp2) otp[k] = 5;
    if(edum > ert && edum < etp2) otp[k] = 6;
    if(edum > etp2){
      if(ptch[k] > 0.){
        if(edum < elt && edum < eax) otp[k] = 2;
        if(edum < ert && edum < eax) otp[k] = 1;
      }
      if(ptch[k] < 0.){
        if(edum < elt && edum < eax) otp[k] = 3;
      }
      if(edum < ert && edum > eax) otp[k] = 1;
    }
    if(edum < etp2){
      if(edum < ert && edum < eax) otp[k] = 5;
    }
    if(edum > elt) otp[k] = 3;
  }

  /* -  Stagnation */
  if(edum < etp1 && edum < eax) otp[k] = 7;
  /* - Confined Potato */
  if(edum > eax && edum < etp2) otp[k] = 8;
  /* - Lost Potato */
  if(edum > ert){
    if(edum > eax && edum < etp2) otp[k] = 9;
  }

  return;
}

void class_domain(Config_t* cfg_ptr){
  int k;
  for(k=0; k < cfg_ptr->nprt; k++){
    class_kdomain(cfg_ptr, k);
  }
  return;
}

void pdedp_finalize(Deposition_t* Depo_ptr){
  /* gets average numbers of counts/bin from non empty bins,
     fill in empty bins.

     then nomalizes.

     if pdedp_pdedp lives on card, this is trivial there*/

  printf("-> Finalize pDEDP computation ...\n");
  /* Get indexes of (DE,DPz)=(0,0) bin */
  const int iDE0 = get_bin(0., Depo_ptr->pdedp_varDE, Depo_ptr->pdedp_nbinDE);
  const int iDPz0 = get_bin(0., Depo_ptr->pdedp_varDPz, Depo_ptr->pdedp_nbinDPz);

  /* make sure center bin has exactly DE=0,DPz=0 */
  /* set to zero */
  Depo_ptr->pdedp_varDE[iDE0]=0.;
  Depo_ptr->pdedp_varDPz[iDPz0]=0.;

  const size_t num_nbins = (size_t)Depo_ptr->pdedp_nbinE
      * (size_t)Depo_ptr->pdedp_nbinPz
      * (size_t)Depo_ptr->pdedp_nbinmu;
  int* agg_nbins = (int*)umacalloc(num_nbins, sizeof(int));
  double* agg_cnt = (double*)umacalloc(num_nbins, sizeof(double));

  /*  */
#ifdef __NVCC__
  /* cuda stuff */
  /* in this particular case we know the dimensions we care about */
  dim3 dimBlock(1,
                (unsigned)Depo_ptr->pdedp_nbinPz,
                (unsigned)Depo_ptr->pdedp_nbinE);

  dim3 dimGrid((unsigned)Depo_ptr->pdedp_nbinmu,1);

  pdedp_finalize_agg_dev<<<dimGrid, dimBlock>>>(Depo_ptr, iDE0, iDPz0, agg_nbins, agg_cnt);
  HANDLE_ERROR(cudaPeekAtLastError());
  /* you might need this for UVM... */
  HANDLE_ERROR(cudaDeviceSynchronize());

#else
  pdedp_finalize_agg_host(Depo_ptr, iDE0, iDPz0, agg_nbins, agg_cnt);
#endif

  /* SYNC here across all blocks, before normalize */

  int total_nbins = 0;
  double total_cnt = 0;
  for(size_t ind=0; ind < num_nbins; ind++){
    total_nbins += agg_nbins[ind];
    total_cnt += agg_cnt[ind];
  }

  /* launch our norm kern */
#ifdef __NVCC__
  pdedp_finalize_norm_dev<<<dimGrid, dimBlock>>>(Depo_ptr, total_nbins, total_cnt, agg_cnt);
  HANDLE_ERROR(cudaPeekAtLastError());
  /* you might need this for UVM... */
  HANDLE_ERROR(cudaDeviceSynchronize());
#else
  pdedp_finalize_norm_host(Depo_ptr, total_nbins, total_cnt, agg_cnt);
#endif

}

void pdedp_finalize_agg_host(Deposition_t* Depo_ptr,
                         int iDE0, int iDPz0,
                         int* agg_nbins, double* agg_cnt){

  for(int iE=0; iE < Depo_ptr->pdedp_nbinE; iE++){
    for(int iPz=0; iPz < Depo_ptr->pdedp_nbinPz; iPz++){
      for(int imu=0; imu < Depo_ptr->pdedp_nbinmu; imu++){
        pdedp_finalize_kernel_agg(Depo_ptr, iE, iPz, imu, iDE0, iDPz0,
                                  agg_nbins, agg_cnt);
      }
    }
  }
}

#ifdef __NVCC__
__global__
void pdedp_finalize_agg_dev(Deposition_t* Depo_ptr,
                        int iDE0, int iDPz0,
                        int* agg_nbins, double* agg_cnt){
  /* one element in grid,  3D blocks */
  const int iE = blockIdx.z * blockDim.z + threadIdx.z;
  if(iE >= Depo_ptr->pdedp_nbinE) return;  /* check bounds */
  const int iPz = blockIdx.y * blockDim.y + threadIdx.y;
  if(iPz >= Depo_ptr->pdedp_nbinPz) return; /* check bounds */
  const int imu = blockIdx.x * blockDim.x + threadIdx.x;
  if(imu >= Depo_ptr->pdedp_nbinmu) return; /* check bounds */

  pdedp_finalize_kernel_agg(Depo_ptr, iE, iPz, imu, iDE0, iDPz0,
      agg_nbins, agg_cnt);


}
#endif

#ifdef __NVCC__
__host__ __device__
#endif
void pdedp_finalize_kernel_agg(Deposition_t* Depo_ptr, int iE, int iPz, int imu,
                               int iDE0, int iDPz0,
                               int* agg_nbins, double* agg_cnt){
  /* gets average numbers of counts/bin from non empty bins,
     fill in empty bins.

     then nomalizes.

     if pdedp_pdedp lives on card, this is trivial there*/
  int ind;
  double* const pdedp_pdedp = Depo_ptr->pdedp_pdedp;

  /*       Get average number of counts/bin from non-empty bins
           and fill in empty bins */
  int cnt_;
  int nbins=0;

  cnt_=0.;
  for(int iDE=0; iDE < Depo_ptr->pdedp_nbinDE; iDE++){
    for(int iDPz=0; iDPz < Depo_ptr->pdedp_nbinDPz; iDPz++){
      ind = get_pdedp_ind(Depo_ptr, iE, iPz, imu, iDE, iDPz);
      cnt_ += pdedp_pdedp[ind];
    }
  }


  if(cnt_ > 0) {
    /* update*/
    nbins += 1;
  } else {
    ind = get_pdedp_ind(Depo_ptr, iE, iPz, imu, iDE0, iDPz0);
    pdedp_pdedp[ind] = 1.;
  }


  /* scatter results*/
  agg_nbins[get_epmu_ind(Depo_ptr, iE, iPz, imu)] = nbins;
  agg_cnt[get_epmu_ind(Depo_ptr, iE, iPz, imu)] = cnt_;

  return;
}

void pdedp_finalize_norm_host(Deposition_t* Depo_ptr,
                              int total_nbins, double total_cnt, double* agg_cnt){
  for(int iE=0; iE < Depo_ptr->pdedp_nbinE; iE++){
    for(int iPz=0; iPz < Depo_ptr->pdedp_nbinPz; iPz++){
      for(int imu=0; imu < Depo_ptr->pdedp_nbinmu; imu++){
        pdedp_finalize_kernel_norm(Depo_ptr, iE, iPz, imu,
                                   total_nbins, total_cnt, agg_cnt);
      }
    }
  }
}

#ifdef __NVCC__
__global__
void pdedp_finalize_norm_dev(Deposition_t* Depo_ptr,
                             int total_nbins, double total_cnt, double* agg_cnt){

  /* one element in grid,  3D blocks */
  const int iE = blockIdx.z * blockDim.z + threadIdx.z;
  if(iE >= Depo_ptr->pdedp_nbinE) return;  /* check bounds */
  const int iPz = blockIdx.y * blockDim.y + threadIdx.y;
  if(iPz >= Depo_ptr->pdedp_nbinPz) return; /* check bounds */
  const int imu = blockIdx.x * blockDim.x + threadIdx.x;
  if(imu >= Depo_ptr->pdedp_nbinmu) return; /* check bounds */

  pdedp_finalize_kernel_norm(Depo_ptr, iE, iPz, imu,
                             total_nbins, total_cnt, agg_cnt);
}
#endif

#ifdef __NVCC__
__host__ __device__
#endif
void pdedp_finalize_kernel_norm(Deposition_t* Depo_ptr, int iE, int iPz, int imu,
                                int total_nbins, double total_cnt, double* agg_cnt){
  /* then nomalizes. */

  int ind;
  double* const pdedp_pdedp = Depo_ptr->pdedp_pdedp;
  double sum_p;

  if(total_nbins > 0) total_cnt /= total_nbins;

  /* Normalize */
  sum_p = fmax(1., agg_cnt[get_epmu_ind(Depo_ptr, iE, iPz, imu)]);
  for(int iDE=0; iDE < Depo_ptr->pdedp_nbinDE; iDE++){
    for(int iDPz=0; iDPz < Depo_ptr->pdedp_nbinDPz; iDPz++){
      ind = get_pdedp_ind(Depo_ptr, iE, iPz, imu, iDE, iDPz);
      pdedp_pdedp[ind] *= total_cnt/sum_p;
    }
  }

  return;
}



void pdedp_out(Deposition_t* Depo_ptr){
  /* writes out the dist file */
  int ind;
  int item;
  double val;
  FILE *ofp=NULL;
  FILE *ofs=NULL;

  printf("\nOutputting P(DE,DP) file %s\n",  Depo_ptr->pdedp_file);
  ofp = fopen(Depo_ptr->pdedp_file, "w");  /* add f status */
  if (ofp == NULL) {
    fprintf(stderr, "\nCan't open output file %s!\n", Depo_ptr->pdedp_file);
    exit(1);
  }

  if(Depo_ptr->output_sparse){
    printf("\nOutputting P(DE,DP) file in sparse format, %s\n",  Depo_ptr->pdedp_sparse_file);
    ofs = fopen(Depo_ptr->pdedp_sparse_file, "w");  /* add f status */
    if (ofs == NULL) {
      fprintf(stderr, "\nCan't open output file %s!\n", Depo_ptr->pdedp_sparse_file);
      exit(1);
    }
    fprintf(ofs, "# iDE iDPz iE iPz imu val\n");
  }


  /* place holders, apparently this will come from transp sometime
     right now, just mimic old behavior. */
  const int lshot = 123456;
  const char *date = "Apr. 2018";
  const int nd = 5;     /* 5-D data */
  const int nq = 0;     /* unknown parameter */
  const int nr = 6;     /* #of decimal places f13.6 */
  const int np = 0;     /* process code */
  const int ns = 1;    /*  # of scalars */
  const char *dev = "DEV";
  const char *labelx = "DEstep"; /* *20 */
  const char *unitsx = "kev";             /* *10 */
  const char *labely = "DPsteps"; /* *20 */
  const char *unitsy = "";             /* *10 */
  const char *labelu = "Evar"; /* *20 */
  const char *unitsu = "keV";             /* *10 */
  const char *labelv = "Pvar"; /* *20 */
  const char *unitsv = "";             /* *10 */
  const char *labelw = "Muvar";  /* *20 */
  const char *unitsw = "";             /* *10 */
  const char *com = ";----END-OF-DATA-----------------COMMENTS:-----------";
  const char *com2 = "UFILE WRITTEN BY ORBIT, see PDEDP_OUT";
  const char *com3 = "SMOOTHING FACTORS, DELAY FACTORS:";
  const char *com4 = "       NONE";
  const char *com5 = "USER COMMENTS:";
  const char *com6 = "       TEST FILE";

  fprintf(ofp, " %-6d %-4s %2d %2d %2d ;-SHOT #- F(X) DATA -PDEDP_OUT-%s\n",
          lshot, dev, nd, nq, nr, date);

  fprintf(ofp, " %10s          ;-SHOT DATE-  UFILES ASCII FILE SYSTEM\n", date);

  fprintf(ofp, " %3d                    ;-NUMBER OF ASSOCIATED SCALAR QUANTITIES-\n", ns);

  fprintf(ofp,"1.0000E+00                   ;-SCALAR, LABEL FOLLOWS:\n");

  fprintf(ofp, " %-20s %-10s ;-INDEPENDENT VARIABLE LABEL: X-\n", labelx, unitsx);

  fprintf(ofp, " %-20s %-10s ;-INDEPENDENT VARIABLE LABEL: Y-\n", labely, unitsy);

  fprintf(ofp, " %-20s %-10s ;-INDEPENDENT VARIABLE LABEL: U-\n", labelu, unitsu);

  fprintf(ofp, " %-20s %-10s ;-INDEPENDENT VARIABLE LABEL: V-\n", labelv, unitsv);

  fprintf(ofp, " %-20s %-10s ;-INDEPENDENT VARIABLE LABEL: W-\n", labelw, unitsw);

  fprintf(ofp, " %13.6e          ; TSTEPSIM  - TIME STEP USED IN SIMULATION [ms]\n",
          Depo_ptr->pdedp_dtsamp);

  fprintf(ofp, " PROBABILITY DATA              ;-DEPENDENT VARIABLE LABEL-\n");

  fprintf(ofp, " %-10d                   ;-PROC CODE- 0:RAW 1:AVG 2:SM 3:AVG+SM\n", np);

  fprintf(ofp, " %-10d          ;-# OF X PTS-\n", Depo_ptr->pdedp_nbinDE);

  fprintf(ofp, " %-10d          ;-# OF Y PTS-\n", Depo_ptr->pdedp_nbinDPz);

  fprintf(ofp, " %-10d          ;-# OF U PTS-\n", Depo_ptr->pdedp_nbinE);

  fprintf(ofp, " %-10d          ;-# OF V PTS-\n", Depo_ptr->pdedp_nbinPz);

  fprintf(ofp, " %-10d          ;-# OF W PTS- X,Y,U,V,W,F(X,Y,U,V,W) DATA FOLLOW:\n",
          Depo_ptr->pdedp_nbinmu);

  /*       -------------------------------------------------------
           Make sure center bin has exactly DE=0,DPz=0

           Get indexes of (DE,DPz)=(0,0) bin */
  double valdum=0.;
  int iDE0 =  get_bin(valdum, Depo_ptr->pdedp_varDE, Depo_ptr->pdedp_nbinDE);
  int iDPz0 = get_bin(valdum, Depo_ptr->pdedp_varDPz, Depo_ptr->pdedp_nbinDPz);

  /* set to zero */
  Depo_ptr->pdedp_varDE[iDE0]=0.;
  Depo_ptr->pdedp_varDPz[iDPz0]=0.;

  /* Write grid vectors */
  write_grid(Depo_ptr, ofp);

  /* if outputting sparse, write grid there too */
  /* variables DEbins,DPbins */
  if(Depo_ptr->output_sparse){
    write_grid(Depo_ptr, ofs);
  }

  /*       -------------------------------------------------------
           Write  p(DE,DP|E,Pz,mu)=0 for each bin.
           This is a big chunk of data with many
           nested, horrible, unelegant and
           inefficient 'for' loops.

           NOTE: at this point, the p(DE,DP)'s are NOT
           normalized. This is to make it easier
           to update the distributions with multiple
           (serial) calls with random distributions
           that simply add together.  */

  for(int iE=0; iE < Depo_ptr->pdedp_nbinE; iE++){
    for(int iPz=0; iPz < Depo_ptr->pdedp_nbinPz; iPz++){
      for(int imu=0; imu < Depo_ptr->pdedp_nbinmu; imu++){

        /*
          compute normalization factor for this bin
          !          pnorm=0.0
          !          do iDE=1,pdedp_nbinDE
          !            do iDPz=1,pdedp_nbinDPz
          !              pnorm=pnorm+pdedp_Pdedp(iDE,iDPz,iE,iPz,imu)
          !            enddo
          !          enddo
          !
          !          ! write normalized p(DE,DPz|E,Pz,mu)
        */

        for(int iDE=0; iDE < Depo_ptr->pdedp_nbinDE; iDE++){
          /*
            ! normalize total probability to 1
            !            if(pnorm.gt.0.) then
            !              do iDPz=1,pdedp_nbinDPz
            !                pdedp_Pdedp(iDE,iDPz,iE,iPz,imu)=
            !     >             pdedp_Pdedp(iDE,iDPz,iE,iPz,imu)/pnorm
            !               enddo
            !            endif
          */

          item = 0;
          for(int iDPz=0; iDPz < Depo_ptr->pdedp_nbinDPz; iDPz++){
            ind = get_pdedp_ind(Depo_ptr, iE, iPz, imu, iDE, iDPz);
            val = Depo_ptr->pdedp_pdedp[ind];
            /* adhere to file format circa 1950s */
            if(item>0 && (item % 6 == 0)) fprintf(ofp,"\n");
            fprintf(ofp, "%14.6e", val);

            /*optionally write out the sparse record  */
            if(Depo_ptr->output_sparse && val != 0.){
              fprintf(ofs, "%d %d %d %d %d %lf.17\n", iE, iPz, imu, iDE, iDPz, val);
            }
            item++;
          }
          fprintf(ofp,"\n");
        }
      }
    }
  }

  fprintf(ofp, "%s\n", com);
  fprintf(ofp, "%s\n", com2);
  fprintf(ofp, "%s\n", com3);
  fprintf(ofp, "%s\n", com4);
  fprintf(ofp, "%s\n", com5);
  fprintf(ofp, "%s\n", com6);

  fclose(ofp);
  if(Depo_ptr->output_sparse){
    fclose(ofs);
  }

  return;
}

#ifdef __NVCC__
__host__ __device__
#endif
static inline int get_res_id_ind(Config_t* cfg_ptr, int kptcl, int time, int i){
  /*  fname indices, sorry */
  const int nsamples = cfg_ptr->depo_ptr->res_id_arr_j;
  const int ni = cfg_ptr->depo_ptr->res_id_arr_i;
  /* we can/should play around with the array order here */
  //const int nprt = cfg_ptr->nprt;

  //return nprt*ni*time + ni*kptcl +i;
  return  kptcl * (nsamples*ni) + time*ni + i;
}

void pdedp_rcrd_resid(Config_t* cfg_ptr, Deposition_t* Depo_ptr){
  /*      pdedp_tskip : reduce loops, skip time steps */
  int j, j2, j3, k, ind, step;
  double dtdum=Depo_ptr->pdedp_tskip *
      1.0E3 * cfg_ptr->dt0 / get_omeg0(cfg_ptr->ptrb_ptr); /* [ms] */
  int jstart = round(0.2 * Depo_ptr->pdedp_dtsamp / dtdum);  /* one indexed! */
  int nsamples= cfg_ptr->nstep_all / Depo_ptr->pdedp_tskip;
  double newE, newPz;
  double Eav, Pzav, Muav;
  double dedum, dpzdum;
  double thisE;
  int iE, iDE, iPz, iDPz, imu;
  double * const res_id_arr = Depo_ptr->res_id_arr;


  printf("   ... pdedp_rcrd_resid computing p(DE,DPz) matrix ...\n");
  int nintv = round( Depo_ptr->pdedp_dtsamp / dtdum);
  int Nav = imax(1, round( Depo_ptr->pdedp_dtav / dtdum));  /* 17, matches */

  for(k=0; k < cfg_ptr->nprt; k++){
    for(j=jstart-1; j < nsamples; j++){
      Eav = 0.;
      Pzav = 0.;
      Muav = 0.;
      for(j3=0; j3 <= (Nav-1); j3++){
        step = j - j3;
        thisE = res_id_arr[get_res_id_ind(cfg_ptr, k, step, 0)];
        Eav  += thisE;  /* E */
        if (thisE <0) {
          printf("neg Energy! %.18lg\n", thisE);
          abort();
        }
        Pzav += res_id_arr[get_res_id_ind(cfg_ptr, k, step, 1)];  /* Pz */
        Muav += res_id_arr[get_res_id_ind(cfg_ptr, k, step, 2)];  /* mu */
      }

      Eav /= ((double)Nav);
      Pzav /= ((double)Nav);
      Muav /= ((double)Nav);

      iE =  get_bin(Eav, Depo_ptr->pdedp_varE, Depo_ptr->pdedp_nbinE);
      iPz =  get_bin(Pzav, Depo_ptr->pdedp_varPz, Depo_ptr->pdedp_nbinPz);
      imu =  get_bin(Muav, Depo_ptr->pdedp_varmu, Depo_ptr->pdedp_nbinmu);

      j2 = nintv;
      ind = get_res_id_ind(cfg_ptr, k, j, 3);  /* xxx check inds */
      if (j + j2 < nsamples &&
          res_id_arr[ind] > 0 &&
          res_id_arr[ind] < 1)
      {
        newE = 0.;
        newPz = 0.;
        for(j3=0; j3 <= (Nav-1); j3++){
          step = j+j2-j3;
          newE  += (res_id_arr[get_res_id_ind(cfg_ptr, k, step, 0)] / ((double)Nav));
          newPz += (res_id_arr[get_res_id_ind(cfg_ptr, k, step, 1)] / ((double)Nav));
        }
        dedum = newE - Eav;
        dpzdum = newPz - Pzav;

        iDE =  get_bin(dedum, Depo_ptr->pdedp_varDE, Depo_ptr->pdedp_nbinDE);
        iDPz =  get_bin(dpzdum, Depo_ptr->pdedp_varDPz, Depo_ptr->pdedp_nbinDPz);
        /* if all bins in range */
        if (newE > 0 &&
            iE >= 0 && iE < Depo_ptr->pdedp_nbinE  &&
            iPz >= 0 && iPz < Depo_ptr->pdedp_nbinPz  &&
            imu >= 0 && imu < Depo_ptr->pdedp_nbinmu  &&
            iDE >= 0 && iDE < Depo_ptr->pdedp_nbinDE  &&
            iDPz >= 0 && iDPz < Depo_ptr->pdedp_nbinDPz
            )
        {
          ind = get_pdedp_ind(Depo_ptr, iE, iPz, imu, iDE, iDPz);
          Depo_ptr->pdedp_pdedp[ind] += 1;
        }

        /* if we want to "optimize" and iE iPz imu bins in range
           we will keep track of max values */
        if (Depo_ptr->pdedp_optimize == 1 &&
            newE > 0 &&
            iE >= 0 && iE < Depo_ptr->pdedp_nbinE  &&
            iPz >= 0 && iPz < Depo_ptr->pdedp_nbinPz &&
            imu >= 0 && imu < Depo_ptr->pdedp_nbinmu)
        {
          Depo_ptr->pdedp_maxDE = fmax(fabs(1.05 * dedum),
                                     Depo_ptr->pdedp_maxDE);
          Depo_ptr->pdedp_maxDPz = fmax(fabs(1.05 * dpzdum),
                                      Depo_ptr->pdedp_maxDPz);
        }
      }  /* res_id 0,1 */

    }  /* j */
  }    /* k */

  /* these are only non zero if the fmax lines above ran for pdedp_optimize */
  if(Depo_ptr->pdedp_optimize == 1){
    printf("pdedp_rcrd_resid::pdedp_maxDE %g\n", Depo_ptr->pdedp_maxDE);
    printf("pdedp_rcrd_resid::pdedp_maxDPz %g\n", Depo_ptr->pdedp_maxDPz);
  }

  return;
}

void rcrd_bfield(Config_t* cfg_ptr, Deposition_t* Depo_ptr){
  /*  record P, Bfield(P,0), Bfield(P,pi)  */
  double pdum;
  FILE *ofp = fopen(Depo_ptr->bfield_file, "w");  /* add f status */
  if (ofp == NULL) {
    fprintf(stderr, "\nCan't open output file %s!\n", Depo_ptr->bfield_file);
    exit(1);
  }
  printf("\nOutputting bfield file %s\n",  Depo_ptr->bfield_file);

  const double pw = get_pw(cfg_ptr->eqlb_ptr);
  Equilib_t* eqlb_ptr = cfg_ptr->eqlb_ptr;
  fprintf(ofp, "P, Bfield(P,0), Bfield(P,pi)\n");
  for(int j=1; j<=200; j++){
    pdum = (-1. + .02*j)*pw;
    fprintf(ofp, "%f %f %f\n",
            pdum,
            bfield(eqlb_ptr, pdum, 0.),
            bfield(eqlb_ptr, pdum, M_PI));
  }

  return;
}

void pdedp_checkbdry(Config_t* cfg_ptr, Deposition_t* depo_ptr){
  /* Check the range used for compute p(DE,DP | E, Pz,mu)
     and adjust the (DE,DPz) range on-the-fly to optimize sampling.
  */

  /*  */
  Equilib_t* eqlb_ptr = cfg_ptr->eqlb_ptr;
  //  double** const GD = get_GD(eqlb_ptr);
  //  double** const B = get_B(eqlb_ptr);
  const double pw = get_pw(eqlb_ptr);
  const double engn = get_engn(cfg_ptr);

  /*don't touch grid if rel. difference is smaller than this */
  const double dthres = 1.;
  int recompute;
  int k;
  double mumax, Pzmin, Pzmax;
  double fctE, fctPz, fctMu;
  double stp;

  printf("Checking pDEDP boundaries ...\n");

  /* First, find boundary for mu and Pz based on
     energy range used in pDEDP calculation */

  /* maximum energy, normalized units */

  /* XXX ask MP which engn he cares about here...  */
  /* double pdedp_engn = depo_ptr->pdedp_Emax * 10.533 * */
  /*     get_prot(cfg_ptr->ptcl_ptr) * pow(GD[0][0],2) / */
  /*     (get_rmaj(eqlb_ptr) * get_zprt(cfg_ptr->ptcl_ptr) * get_bkg(cfg_ptr) * pow(B[0][0],2)); */
  /* /\* cf. engn definition in initial.f *\/ */

  /* min/max B field */ /*checked okay*/
  /*const double Bmn = bfield(eqlb_ptr, pw, 0.);
  const double Bmx = bfield(eqlb_ptr, pw, M_PI); */
  /* MP: get Bmin, Bmax from equilibrium initialization */
  const double Bmn = cfg_ptr->bmin;
  const double Bmx = cfg_ptr->bmax; 
  
  /* redefine range of mu
     add buffer to the actual range */
  mumax = 1./Bmn * (depo_ptr->pdedp_nbinmu + 1.) / depo_ptr->pdedp_nbinmu;
  mumax = 1.05*mumax;  /* add buffer region*/
  
  /* upper limit for Pz, add buffer */
  Pzmax = gfun(eqlb_ptr, 0)/pw*sqrt(2. * engn) *
      (depo_ptr->pdedp_nbinPz + 1.) / depo_ptr->pdedp_nbinPz;
  Pzmax = 1.05*Pzmax;  /* add buffer region*/

  /* lower limit for Pz, add buffer */
  Pzmin = -1. - gfun(eqlb_ptr, pw)/pw*sqrt(2. * engn) / Bmx *
      (depo_ptr->pdedp_nbinPz + 1.) /depo_ptr->pdedp_nbinPz;
  Pzmin=1.05*Pzmin; /* add buffer region*/
  

  /* Check wheter the Pz,mu range needs to be adjusted. */

  fctMu = (depo_ptr->pdedp_mumax - mumax ) / depo_ptr->pdedp_mumax;
  fctPz = fmax(fabs((depo_ptr->pdedp_Pzmax - Pzmax) / depo_ptr->pdedp_Pzmax),
               fabs((depo_ptr->pdedp_Pzmax - Pzmin) / depo_ptr->pdedp_Pzmin));

  if(fctMu > dthres || fctPz > dthres ||
     Pzmin < depo_ptr->pdedp_Pzmin || Pzmax > depo_ptr->pdedp_Pzmax ||
     mumax > depo_ptr->pdedp_mumax) {
    /* update */

    /* display info with updated grid */
    printf("  -> New Pz,mu grid computed:\n");
    printf("  original: Pz1= %f,    Pz2= %f,    mu=%f\n",
           depo_ptr->pdedp_Pzmin,
           depo_ptr->pdedp_Pzmax,
           depo_ptr->pdedp_mumax);
    printf("  updated: Pz1= %f,    Pz2= %f,    mu=%f\n",
           Pzmin, Pzmax, mumax);

    depo_ptr->pdedp_Pzmax = Pzmax;
    depo_ptr->pdedp_Pzmin = Pzmin;
    depo_ptr->pdedp_mumin = 0.;
    depo_ptr->pdedp_mumax = mumax;


    /* Pz */
    stp=(depo_ptr->pdedp_Pzmax - depo_ptr->pdedp_Pzmin)/depo_ptr->pdedp_nbinPz;
    for(k=0; k < depo_ptr->pdedp_nbinPz; k++){
      depo_ptr->pdedp_varPz[k] = depo_ptr->pdedp_Pzmin + ((double)k)* stp + stp/2.;
    }

    /* mu Bo/E */
    stp=(depo_ptr->pdedp_mumax - depo_ptr->pdedp_mumin) / depo_ptr->pdedp_nbinmu;
    for(k=0; k < depo_ptr->pdedp_nbinmu; k++){
      depo_ptr->pdedp_varmu[k] = depo_ptr->pdedp_mumin +((double)k) * stp + stp/2.;
    }

    recompute = 1;
  } else {
    printf(" -> Pz,mu grid looks OK - \n\n");
  }

  /* Check wheter the DE,DPz range needs to be
     adjusted. */
  fctE=(depo_ptr->pdedp_maxDE - depo_ptr->pdedp_DEmax) / depo_ptr->pdedp_DEmax;
  fctPz=(depo_ptr->pdedp_maxDPz - depo_ptr->pdedp_DPzmax) / depo_ptr->pdedp_DPzmax;

  if(fabs(fctE) < dthres && fabs(fctPz) < dthres &&
     fctE <= 1  && fctPz <= 1){
    printf("  -> DE,DPz grid looks OK - \n\n");
  } else {
    printf("updating DE,DPz range \n");

    /* new values */
    depo_ptr->pdedp_DEmax = (1. + fctE) * depo_ptr->pdedp_DEmax;
    depo_ptr->pdedp_DPzmax = (1. + fctPz) * depo_ptr->pdedp_DPzmax;

    /* round off [doesn't need to preserve precision] */
    depo_ptr->pdedp_DEmax = 1E-2 * (ceil(1E2 * depo_ptr->pdedp_DEmax));
    depo_ptr->pdedp_DPzmax = 1E-4 * (ceil(1E4 * depo_ptr->pdedp_DPzmax));

    /* symmetric grid */
    depo_ptr->pdedp_DEmin = -depo_ptr->pdedp_DEmax;
    depo_ptr->pdedp_DPzmin = -depo_ptr->pdedp_DPzmax;


    /*  define new grid */
    /* Delta E */
    stp = 2. * depo_ptr->pdedp_DEmax / depo_ptr->pdedp_nbinDE;
    for(k=0; k < depo_ptr->pdedp_nbinDE; k++){
      depo_ptr->pdedp_varDE[k] = -depo_ptr->pdedp_DEmax + ((double)k) * stp + stp/2.;
    }
    /* Delta Pz */
    stp = 2. * depo_ptr->pdedp_DPzmax / depo_ptr->pdedp_nbinDPz;
    for(k=0; k<depo_ptr->pdedp_nbinDPz; k++){
      depo_ptr->pdedp_varDPz[k] = -depo_ptr->pdedp_DPzmax + ((double)k) * stp + stp/2.;
    }

    /* update flag */
    recompute = 1;
    printf("DE,DPz grid was recomputed\n");
  }
  /* If range of Pz, mu, DE, or DPz has changed, */
  /* update pDEDP computation */
  if(recompute != 0){
    /* reset pDEDP before resampling */
    memset(depo_ptr->pdedp_pdedp, 0, sizeof_pdedp(depo_ptr)*sizeof(double));  /* zero */
    /* compute pDEDP probability based on  updated range */
    pdedp_rcrd_resid(cfg_ptr, depo_ptr);
  }

  /* reset flag - no further optimizations */

  depo_ptr->pdedp_optimize = false;

  printf("- done.\n\n");

  return;
}

void deposition(Config_t* cfg_ptr, int irun_pdedp){
  if(get_do_modestep(cfg_ptr->ptrb_ptr) && irun_pdedp==0){
    printf("Reading deposition from file\n");
    sampledep(cfg_ptr);
  }
  else{  /* use Mario's deposition */
    if( irun_pdedp % 2  == 0){
      fulldepmp(cfg_ptr, cfg_ptr->depo_ptr);
    } else {
      fulldepmp_co(cfg_ptr, cfg_ptr->depo_ptr);
    }
  }
  return;
}

void fulldepmp(Config_t* cfg_ptr, Deposition_t* depo_ptr){
  /* loops over k, can live on device one day */

  int k, kd, np2;

  const double einj1 = depo_ptr->pdedp_Emin;  /* [keV] */
  const double einj2 = depo_ptr->pdedp_Emax;
  Particles_t* Ptcl = cfg_ptr->ptcl_ptr;
  double* pol = get_pol(Ptcl);
  double* thet = get_thet(Ptcl);
  double* pot = get_pot(Ptcl);
  double* dt = get_dt(Ptcl);
  double* ptch=get_ptch(Ptcl);
  const double pw = get_pw(cfg_ptr->eqlb_ptr);
  const double engn = get_engn(cfg_ptr);
  double* en = get_en(Ptcl);
  double* rho = get_rho(Ptcl);
  const double dt0 = cfg_ptr->dt0;
  double* zet = get_zet(Ptcl);
  double ekev = get_ekev(Ptcl);
  double* rmu = get_rmu(cfg_ptr->ptcl_ptr);
  double* b = get_b(cfg_ptr->ptcl_ptr);

  printf("fulldepmp::Emin %g Emax %g\n", einj1, einj2);

  /* Full deposition*/
  np2 = .5* cfg_ptr->nprt;

  /* outside-co moving */
  for(kd=0; kd < np2; kd++){
    ptch[kd] = rand_double();  /* rand */
    thet[kd] = 0.;
    dt[kd] = dt0;
    pol[kd] = .001 + .999 * rand_double();
    pol[kd] *= pw;
    zet[kd] = rand_double() * 2. * M_PI;
    /* kinetic energy */
    en[kd] = (einj1 + rand_double() * (einj2 - einj1)) * engn / ekev;
  }

  /* second half, should start at np2 */
  /* -inside-counter moving */
  for(kd=np2 ; kd < cfg_ptr->nprt; kd++){
    ptch[kd] = -rand_double();
    thet[kd] = M_PI;
    dt[kd] = dt0;
    pol[kd] = .001 + .999 * rand_double();
    pol[kd] *= pw;
    zet[kd] = rand_double() * 2.* M_PI;
    /* kinetic energy */
    en[kd] = (einj1 + rand_double() * (einj2 - einj1)) * engn / ekev;
  }

  printf(" fulldepmp deposit\n");
  for(k=0; k<cfg_ptr->nprt; k++){
    kfield(cfg_ptr, k);
    rho[k] = ptch[k]*sqrt(2. * en[k]) /b[k];
    rmu[k] = en[k] / b[k] -
        .5 * rho[k] * rho[k] * b[k];
    en[k] += pot[k];
  }

  /* re-sample lost particles */
  class_domain(cfg_ptr);
  fullredepmp(cfg_ptr, depo_ptr);

  return;
}

void fullredepmp(Config_t* cfg_ptr, Deposition_t* depo_ptr){
  int kd,np2,nlost,nmaxs,imaxs;

  const double einj1 = depo_ptr->pdedp_Emin;  /* [keV] */
  const double einj2 = depo_ptr->pdedp_Emax;
  Particles_t* Ptcl = cfg_ptr->ptcl_ptr;
  double* pol = get_pol(Ptcl);
  double* thet = get_thet(Ptcl);
  double* pot = get_pot(Ptcl);
  double* dt = get_dt(Ptcl);
  double* ptch=get_ptch(Ptcl);
  const double pw = get_pw(cfg_ptr->eqlb_ptr);
  const double engn = get_engn(cfg_ptr);
  double* en = get_en(Ptcl);
  double* rho = get_rho(Ptcl);
  const double dt0 = cfg_ptr->dt0;
  double* zet = get_zet(Ptcl);
  double ekev = get_ekev(Ptcl);
  double* rmu = get_rmu(Ptcl);
  double* b = get_b(Ptcl);
  int* otp = get_otp(Ptcl);

  nmaxs=5E3; /* max number of attemps to redeposit particle */
  /* -  Full deposition, */
  np2 = .5 * cfg_ptr->nprt;

  /* dont print this recursively */
  if(! depo_ptr->recursing){
    printf("Entering FULLREDEPMP...\n");
  }

  nlost=0;

  /* -outside-co moving */
  for(kd=0; kd < np2; kd++){
    imaxs=1;

    if(depo_ptr->pdedp_do_focusdep && imaxs < nmaxs){
      check_res_ptc(cfg_ptr, kd);
      imaxs += 1;
    }

    if(pol[kd] >= pw || otp[kd] == 2 ||
       otp[kd] == 4 || otp[kd] == 6){
      /* lost, replace it */
      /* printf("dbg lost1 %d\n", nlost); */
      /* if(otp[kd] == 2 || otp[kd] == 4 || otp[kd] == 6) printf("otp %d\n", otp[kd]); */

      nlost+=1;
      ptch[kd] = rand_double();
      thet[kd] = 0.;
      dt[kd] = dt0;
      pol[kd] =  .002*pw + .997*pw*rand_double();
      zet[kd] = 2. * M_PI * rand_double();
      /* kinetic energy*/
      en[kd] = (einj1 + rand_double() * (einj2 - einj1)) * engn  /  ekev;
      kfield(cfg_ptr, kd);
      rho[kd] = ptch[kd]*sqrt(2.*en[kd])/b[kd];
      rmu[kd] = en[kd]/b[kd] - .5*rho[kd]*rho[kd]*b[kd];
      en[kd] += pot[kd];
    }
  }

  /* -inside-counter moving */
  for(kd=np2; kd < cfg_ptr->nprt ; kd++)
  {
    imaxs=1;

    if(depo_ptr->pdedp_do_focusdep && imaxs < nmaxs){
      check_res_ptc(cfg_ptr, kd);
      imaxs += 1;
    }
    if(pol[kd] >= pw || otp[kd] == 2 ||
       otp[kd] == 4 || otp[kd] == 6){
      /* lost, replace it */
      /* printf("dbg lost2 %d\n", nlost);  */
      /* if(otp[kd] == 2 || otp[kd] == 4 || otp[kd] == 6){ */
      /*   printf("otp %d\n", otp[kd]); */
      /* }else{ */
      /*   printf("pol %g >= pw %g\n", pol[kd], pw); */
      /* } */

      nlost += 1;
      ptch[kd] = -rand_double();
      thet[kd] = M_PI;
      dt[kd] = dt0;
      pol[kd] = .002 * pw + .997 *pw * rand_double();
      zet[kd] =  2. * M_PI * rand_double();
      /* kinetic energy */
      en[kd] = (einj1 + rand_double() * (einj2-einj1)) * engn / ekev;
      kfield(cfg_ptr, kd);
      rho[kd] = ptch[kd]*sqrt(2.*en[kd])/b[kd];
      rmu[kd] = en[kd]/b[kd] - .5*rho[kd]*rho[kd]*b[kd];
      en[kd] += pot[kd];
    }
  }
  if (nlost > 0){
    printf("- Number of lost particles: %d -> iterate sampling...\n", nlost);
    class_domain(cfg_ptr);
    depo_ptr->recursing = true;
    fullredepmp(cfg_ptr, depo_ptr);   /* goto 11 */
  }
  /* the prodigal particle has returned */
  depo_ptr->recursing = false;

  return;
}

void check_res_ptc(Config_t* cfg_ptr, int kd){
  int j, k, ind;
  double ptot, pmax;
  double edum, pzdum, mudum;
  double tmp;

  Particles_t* Ptcl = cfg_ptr->ptcl_ptr;
  Deposition_t* depo_ptr = cfg_ptr->depo_ptr;
  double* pol = get_pol(Ptcl);
  const double pw = get_pw(cfg_ptr->eqlb_ptr);
  const double engn = get_engn(cfg_ptr);
  double* en = get_en(Ptcl);
  double* rho = get_rho(Ptcl);
  double ekev = get_ekev(Ptcl);
  double* rmu = get_rmu(Ptcl);
  double* g = get_g(Ptcl);

  ptot = 0.;
  pmax = 0.;

  edum = en[kd]*ekev/engn;
  pzdum=(g[kd]*rho[kd] - pol[kd])/pw;
  mudum=rmu[kd]/en[kd];

  const int iE = get_bin(edum, depo_ptr->pdedp_varE, depo_ptr->pdedp_nbinE);
  const int iPz = get_bin(pzdum, depo_ptr->pdedp_varPz, depo_ptr->pdedp_nbinPz);
  const int iMu = get_bin(mudum, depo_ptr->pdedp_varmu, depo_ptr->pdedp_nbinmu);

  if(iE < 0 || iE >= depo_ptr->pdedp_nbinE ||
     iPz < 0 || iPz >= depo_ptr->pdedp_nbinPz ||
     iMu < 0 || iMu >= depo_ptr->pdedp_nbinmu){
     /* Out of domain, exclude from computation */
    /* printf("DBG check_res_ptc bin bounds. %g iE=%d/%d %g iPz=%d/%d %g iMu=%d/%d\n",
           edum, iE, depo_ptr->pdedp_nbinE,
           pzdum, iPz, depo_ptr->pdedp_nbinPz,
           mudum, iMu, depo_ptr->pdedp_nbinmu); */
    pol[kd] = 2. * pw;
    return;
  }
  /* else, valid bin - proceed */
  for(j=0; j < depo_ptr->pdedp_nbinDE; j++){
    for(k=0; k < depo_ptr->pdedp_nbinDPz; k++){
      ind = get_pdedp_ind(depo_ptr, iE, iPz, iMu, j, k);
      tmp = depo_ptr->pdedp_pdedp[ind];
      ptot += tmp;

      /* update max val seen */
      if(pmax < tmp){
        pmax=tmp;
      }
    }
  }

  if(ptot <= pmax){
    /* throw away particle */
    pol[kd] = 2. * pw;
  }

  return;
}

void fulldepmp_co(Config_t* cfg_ptr, Deposition_t* depo_ptr){
  /*    all confined orbits, broad energy range, co- only */
  int k, kd;

  const double einj1 = depo_ptr->pdedp_Emin;  /* [keV] */
  const double einj2 = depo_ptr->pdedp_Emax;
  Particles_t* Ptcl = cfg_ptr->ptcl_ptr;
  double* ptch=get_ptch(Ptcl);
  double* thet = get_thet(Ptcl);
  double* dt = get_dt(Ptcl);
  double* pol = get_pol(Ptcl);
  double* pot = get_pot(Ptcl);
  double* zet = get_zet(Ptcl);
  double* en = get_en(Ptcl);
  const double engn = get_engn(cfg_ptr);
  double ekev = get_ekev(Ptcl);
  double* rho = get_rho(Ptcl);
  double* rmu = get_rmu(Ptcl);
  const double dt0 = cfg_ptr->dt0;
  const double pw = get_pw(cfg_ptr->eqlb_ptr);
  double* b = get_b(Ptcl);

  /*   Full deposition,
       outside-co moving */
  for(kd=0; kd < cfg_ptr->nprt; kd++){
    ptch[kd] = rand_double();
    thet[kd] = 0.;
    dt[kd] = dt0;
    pol[kd] =  .002*pw + .996 * pw * rand_double();
    zet[kd]=2. * M_PI * rand_double();    /*  random toroidal */
    /* kinetic energy */
    en[kd] = (einj1 + rand_double() * (einj2 - einj1)) * engn / ekev;
  }

  printf("%d fulldepmp_co deposit\n", kd);
  for(k=0; k < cfg_ptr->nprt; k++){
    kfield(cfg_ptr, k);
    rho[k] = ptch[k]*sqrt(2.*en[k])/b[k];
    rmu[k] = en[k]/b[k] - .5*rho[k]*rho[k]*b[k];
    en[k] += pot[k];
  }

  /* re-sample lost particles */
  class_domain(cfg_ptr);
  fullredepmp(cfg_ptr, depo_ptr);   /* goto 11 */

  return;
}

#ifdef __NVCC__
__host__ __device__
#endif
void rcrd_vararr(Config_t* cfg_ptr, int k, int step){
  /* for particle k at step , stash some values in res_id_arr */
  /* other variables saved in konestep*/
  Particles_t* ptcl_ptr = cfg_ptr->ptcl_ptr;
  Deposition_t* depo_ptr = cfg_ptr->depo_ptr;
  double* pol = get_pol(ptcl_ptr);
  const double pw = get_pw(cfg_ptr->eqlb_ptr);
  const double engn = get_engn(cfg_ptr);
  double* en = get_en(ptcl_ptr);
  double* rho = get_rho(ptcl_ptr);
  const double ekev = get_ekev(ptcl_ptr);
  double* rmu = get_rmu(ptcl_ptr);
  double* g = get_g(ptcl_ptr);
  double* res_id_arr = depo_ptr->res_id_arr;

  res_id_arr[get_res_id_ind(cfg_ptr, k, step, 0)] = en[k]*ekev/engn; /* E */
  res_id_arr[get_res_id_ind(cfg_ptr, k, step, 1)] = (g[k]*rho[k] - pol[k])/pw; /* Pz */
  res_id_arr[get_res_id_ind(cfg_ptr, k, step, 2)] = rmu[k]/en[k]; /* mu */
  res_id_arr[get_res_id_ind(cfg_ptr, k, step, 3)] = pol[k]/pw;

  return;
}
