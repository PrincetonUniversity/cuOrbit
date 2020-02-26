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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

/* for constants like pi */
#include <math.h>

#include "orbit_config_api.h"
#include "orbit_perturbation.h"
#include "orbit_particles.h"
#include "orbit_util.h"
#include "cuda_helpers.h"

const int NAMP_ = 155;

struct Perturb {
  int nflr;
  int lpt;
  int md1;
  int md2;
  int modes;
  double psi_solRZ;

  double falf;
  double ascale;
  double global_scaling_factor;
  double alimit;
  double freq_scaling_factor;
  bool do_modestep;
  double omeg0;
  //double sng;

  int *mmod;  /* pol mode numbers */
  int *nmod;  /* pol mode numbers */
  double *omegv;  /* mode frequency */
  double *amp;    /* amp */
  double *damp;
  double *harm;
  double *xx;
  double *phaz;
  /* displacements xi1-xi3 */
  double *xi1, *xi2, *xi3;
  /* alphas a1-a3*/
  double *a1, *a2, *a3;

  /* modestepping */
  double *damdt;
  double *dpzdt;
  double *dp;
  double *dadm;

};

Perturb_t* Perturb_ctor(){
  return (Perturb_t*)umacalloc(1, sizeof(Perturb_t));
}


void initialize_Perturb(Perturb_t* ptrb_ptr, Config_t* config_ptr,
                        Equilib_t* equilib_ptr, Particles_t* ptcl_ptr){


  /* first we set values expected from config */
  ptrb_ptr->nflr = 1;  /* xxx, is this configurable? */
  ptrb_ptr->falf = config_ptr->falf;
  ptrb_ptr->ascale = config_ptr->ascale;
  ptrb_ptr->alimit = config_ptr->alimit;
  ptrb_ptr->global_scaling_factor = config_ptr->global_scaling_factor;
  ptrb_ptr->freq_scaling_factor = config_ptr->freq_scaling_factor;
  ptrb_ptr->do_modestep = false;
  if(config_ptr->do_modestep) ptrb_ptr->do_modestep = true;
  //ptrb_ptr->sng = config_ptr->sng;

  /* in set1 */
  //ptrb_ptr->omeg0 = 9.58E6 * get_zprt(ptcl_ptr) * config_ptr->bkg / get_prot(ptcl_ptr);

  int j, k, m, n, ind, md;
  int idum;
  int nmd, mmin, mmax, ndum;
  double fkhz;
  double omrat;
  double xxmax;
  double px;
  int lptm1;
  size_t sz;

  FILE *ifp;
  const char *mode = "r";
  ifp = fopen(config_ptr->displ_file, mode);
  if (ifp == NULL) {
    fprintf(stderr, "Can't open input file %s!\n", config_ptr->displ_file);
    exit(1);
  }
  printf("Parsing Displacement file %s\n",  config_ptr->displ_file);

  /* file header lines */
  fscanf(ifp, "%*[^\n]\n");  /* skip 0 */
  fscanf(ifp, "%*[^\n]\n");  /* skip 1 */
  fscanf(ifp, "%*[^\n]\n");  /* skip 2 */

  fscanf(ifp, "%d %d %d %d %lf %d ", &(ptrb_ptr->lpt), &nmd, &mmin, &mmax, &omrat, &ndum);
  fscanf(ifp, "%*[^\n]\n");  /*  skip remaing part of line */
  fkhz = omrat * ptrb_ptr->falf;
  lptm1 = ptrb_ptr->lpt - 1;
  printf("lpt = %d\nnmd = %d\nmmin = %d\nmmax = %d\nomrat = %g\nfkhz = %g\n",
         ptrb_ptr->lpt, nmd, mmin, mmax, omrat, fkhz);

  fscanf(ifp, "%*[^\n]\n");  /* skip 4 */
  fscanf(ifp, "%d ", &idum);
  assert (idum == ptrb_ptr->lpt);
  for(j=0; j < ptrb_ptr->lpt; j++){
    /* zero md */ //xxx according to code, we can just skip these...overwrite below... confirm with Mario
    //fscanf(ifp, "%lf ", &(ptrb_ptr->xi1[j]));
    fscanf(ifp, "%*f ");
  }

  fscanf(ifp, "%*[^\n]\n");  /* skip*/
  fscanf(ifp, "%d %d ", &idum, &(ptrb_ptr->modes));
  assert(idum == ptrb_ptr->lpt);

  /* malloc xi */
  printf("Malloc'ing arrays for Perturbation\n");
  sz = (unsigned)(ptrb_ptr->lpt * ptrb_ptr->modes);
  ptrb_ptr->xi1 = (double*)umacalloc(sz, sizeof(double));
  ptrb_ptr->xi2 = (double*)umacalloc(sz, sizeof(double));
  ptrb_ptr->xi3 = (double*)umacalloc(sz, sizeof(double));
  ptrb_ptr->a1 = (double*)umacalloc(sz, sizeof(double));
  ptrb_ptr->a2 = (double*)umacalloc(sz, sizeof(double));
  ptrb_ptr->a3 = (double*)umacalloc(sz, sizeof(double));

  ptrb_ptr->xx = (double*)umacalloc((unsigned)ptrb_ptr->lpt, sizeof(double));

  ptrb_ptr->amp = (double*)umacalloc((unsigned)ptrb_ptr->modes, sizeof(double));
  ptrb_ptr->damp = (double*)umacalloc((unsigned)ptrb_ptr->modes, sizeof(double));
  ptrb_ptr->harm = (double*)umacalloc((unsigned)ptrb_ptr->modes, sizeof(double));
  ptrb_ptr->omegv = (double*)umacalloc((unsigned)ptrb_ptr->modes, sizeof(double));
  ptrb_ptr->phaz = (double*)umacalloc(
      (unsigned)ptrb_ptr->modes * (unsigned)config_ptr->nprt,
      sizeof(double));  /* xxx check size , and pmegv*/

  ptrb_ptr->mmod = (int*)umacalloc((unsigned)ptrb_ptr->modes, sizeof(int));
  ptrb_ptr->nmod = (int*)umacalloc((unsigned)ptrb_ptr->modes, sizeof(int));

  /* modestep */
  ptrb_ptr->damdt = (double*)umacalloc((unsigned)ptrb_ptr->modes, sizeof(double));
  ptrb_ptr->dpzdt = (double*)umacalloc((unsigned)ptrb_ptr->modes, sizeof(double));
  ptrb_ptr->dp = (double*)umacalloc((unsigned)ptrb_ptr->modes, sizeof(double));
  ptrb_ptr->dadm = (double*)umacalloc((unsigned)ptrb_ptr->modes, sizeof(double));


  for(md=0; md < ptrb_ptr->modes; md++){
    fscanf(ifp, "** m = %d ", &(ptrb_ptr->mmod[md]));
    assert(md == ptrb_ptr->mmod[md]);
    for(j=0; j < ptrb_ptr->lpt; j++){
      ind = ptrb_ptr->lpt * md + j;
      fscanf(ifp, "%lf ", &(ptrb_ptr->xi1[ind]));
    }  /* j */
  }  /* md */

  fclose(ifp);

  /* now from here we initialize remaining structure */

  for(md=0; md < ptrb_ptr->modes; md++){
    ptrb_ptr->harm[md] = 1;
    ptrb_ptr->omegv[md] = 2E3 * M_PI * fkhz / ptrb_ptr->omeg0;
    ptrb_ptr->nmod[md] = nmd;
  }

  printf("Change Amplitude = %g\nChange Frequency = %g\nLimit = %g\n",
         ptrb_ptr->global_scaling_factor, ptrb_ptr->freq_scaling_factor, ptrb_ptr->alimit);

  md = -1;
  for(k=0; k < ptrb_ptr->modes; k++){
    for(j=0; j < ptrb_ptr->lpt; j++){
      ind = ptrb_ptr->lpt * k + j;
      ptrb_ptr->xx[j] = fabs(ptrb_ptr->xi1[ind]);
    }
    xxmax = darray_max(ptrb_ptr->xx, (unsigned)ptrb_ptr->lpt);
    ptrb_ptr->damp[k] = xxmax;

    if (ptrb_ptr->damp[k] < ptrb_ptr->alimit){

      continue;
    }
    if (ptrb_ptr->mmod[k] == 0) continue;

    md++;
    ptrb_ptr->amp[md] = ptrb_ptr->ascale * xxmax * ptrb_ptr->global_scaling_factor; /* modify  scaling */
    ptrb_ptr->omegv[md] = ptrb_ptr->omegv[k] * ptrb_ptr->freq_scaling_factor; /* modify freq */
    ptrb_ptr->mmod[md] = ptrb_ptr->mmod[k];
    ptrb_ptr->nmod[md] = ptrb_ptr->nmod[k];
    m = ptrb_ptr->mmod[md];
    n = ptrb_ptr->nmod[md];
    for(j=0; j < ptrb_ptr->lpt; j++){
      ind = ptrb_ptr->lpt * md + j;
      /*  normalize, careful inds */
      ptrb_ptr->xi1[ind] = ptrb_ptr->xi1[ ptrb_ptr->lpt * k + j] /xxmax;

      /* compute alpha from disp */
      px = j * get_pw(equilib_ptr) / lptm1;
      ptrb_ptr->a1[ind] = ( m / qfun(equilib_ptr, px) - n) * ptrb_ptr->xi1[ind] / (
          m * gfun(equilib_ptr, px) + n * rifun(equilib_ptr, px));
    }  /* j */
    printf("%d %d %d %f %f %f\n",
           md, ptrb_ptr->nmod[md], ptrb_ptr->mmod[md], 1E5*ptrb_ptr->amp[md],
           ptrb_ptr->omegv[md] * (ptrb_ptr->omeg0) / 6280.,
           xxmax);
  }  /* k */

  /*  does this really need 1 based?, lets try without*/
  ptrb_ptr->md1 = 0;
  /* and do we need md2 at all? looks vestigile, lets try without*/
  ptrb_ptr->md2 = ptrb_ptr->modes;

  splna(ptrb_ptr, equilib_ptr, ptcl_ptr);
  splnx(ptrb_ptr, equilib_ptr, ptcl_ptr);
  return;
}


void splna(Perturb_t* ptrb_ptr, Equilib_t* equilib_ptr, Particles_t* ptcl_ptr){
  int ind, j, m, md;
  int jm, jp, jpp;
  const int lpt = ptrb_ptr->lpt;
  const int lptm = lpt - 1;
  const double dpx = get_pw(equilib_ptr) / (double)lptm;

  double* pol = get_pol(ptcl_ptr);

  const int lpx = 1;  /* mp change mar 2016 */

  for(md=0; md < ptrb_ptr->modes; md++){
    m = ptrb_ptr->mmod[md];
    ind = md * lpt;
    for(j=0; j < lpx; j++){
      ptrb_ptr->xi1[ind + j] = pow(pol[j], m) *
          ptrb_ptr->xi1[md*lpt + (lpx-1)] /
          pow(pol[lpx-1], m);
    }
    ptrb_ptr->a2[md*lpt] = (
        10. * ptrb_ptr->a1[md*lpt + 1] -
        7.  * ptrb_ptr->a1[md*lpt] -
        3.  * ptrb_ptr->a1[md*lpt + 2]
                            ) / (4. * dpx);
    if (m != 1){
      ptrb_ptr->a2[md*lpt] = 0;
    }
    for(j=1; j<lptm; j++){
      jm = j-1;
      jp = j+1;
      jpp = imin(j+2, lptm);

      ptrb_ptr->a2[ind + j] =
          -1. * ptrb_ptr->a2[ind + jm] +
          2. * (ptrb_ptr->a1[ind + j] - ptrb_ptr->a1[ind + jm]) / dpx;
      /* smooth a1 */

      ptrb_ptr->a1[ind + jp] =
          0.4 * dpx * ptrb_ptr->a2[ind + j] +
          0.3 * ptrb_ptr->a1[ind + jpp] +
          0.7 * ptrb_ptr->a1[ind + j];
    }  /* j */

    for(j=0; j<lptm; j++){
      ptrb_ptr->a3[ind + j] = (
          ptrb_ptr->a2[ind + j + 1] - ptrb_ptr->a2[ind + j] ) / (2. * dpx);
    }  /* j */
  }    /* md */
  return;
}

void splnx(Perturb_t* ptrb_ptr, Equilib_t* equilib_ptr, Particles_t* ptcl_ptr){
  int ind, j, m, md;
  int jm, jp, jpp;
  const int lpt = ptrb_ptr->lpt;
  const int lptm = lpt - 1;
  const double dpx = get_pw(equilib_ptr) / (double)lptm;

  double* pol = get_pol(ptcl_ptr);

  const int lpx = 1;  /* mp change mar 2016 */

  for(md=0; md<ptrb_ptr->modes; md++){
    m = ptrb_ptr->mmod[md];
    ind = md * lpt;
    for(j=0; j<lpx; j++){
      ptrb_ptr->xi1[ind + j] = pow(pol[j], m) * ptrb_ptr->xi1[ind + (lpx-1) ] /
          pow(pol[lpx-1], m);
    }
    ptrb_ptr->xi2[ind] =  (
        10. * ptrb_ptr->xi1[ind + 1] -
        7.  * ptrb_ptr->xi1[ind] -
        3.  * ptrb_ptr->xi1[ind + 2]
                           ) / ( 4. * dpx);
    if(m != 1){
      ptrb_ptr->xi2[ind] = 0;
    }

    for(j=1; j<lptm; j++){
      jm = j-1;
      jp = j+1;
      jpp = imin(j+2, lptm);

      ptrb_ptr->xi2[ind + j] =
          -1. * ptrb_ptr->xi2[ind + jm] +
          2. * (ptrb_ptr->xi1[ind + j] - ptrb_ptr->xi1[ind + jm]) / dpx;
      /* smooth a1 */

      ptrb_ptr->xi1[ind + jp] =
          0.4 * dpx * ptrb_ptr->xi2[ind + j] +
          0.3 * ptrb_ptr->xi1[ind + jpp] +
          0.7 * ptrb_ptr->xi1[ind + j];
    }  /* j */

    for(j=0; j<lptm; j++){
      ptrb_ptr->xi3[ind + j] = (
          ptrb_ptr->xi2[ind + j + 1] - ptrb_ptr->xi2[ind + j] ) / (2. * dpx);
    }  /* j */
  }    /* md */
  return;
}

void set_omeg0(Perturb_t* ptrb_ptr, double val){
  ptrb_ptr->omeg0 = val;
}

#ifdef __NVCC__
__host__ __device__
#endif
double get_omeg0(Perturb_t* ptrb_ptr){
    return ptrb_ptr->omeg0;
}

#ifdef __NVCC__
__host__ __device__
#endif
int get_nflr(Perturb_t* ptrb_ptr){
  return ptrb_ptr->nflr;
}

#ifdef __NVCC__
__host__ __device__
#endif
int get_lpt(Perturb_t* ptrb_ptr){
  return ptrb_ptr->lpt;
}

#ifdef __NVCC__
__host__ __device__
#endif
int get_md1(Perturb_t* ptrb_ptr){
  return ptrb_ptr->md1;
}

#ifdef __NVCC__
__host__ __device__
#endif
int get_md2(Perturb_t* ptrb_ptr){
  return ptrb_ptr->md2;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_omegv(Perturb_t* ptrb_ptr){
  return ptrb_ptr->omegv;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_phaz(Perturb_t* ptrb_ptr){
  return ptrb_ptr->phaz;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_a1(Perturb_t* ptrb_ptr){
  return ptrb_ptr->a1;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_a2(Perturb_t* ptrb_ptr){
  return ptrb_ptr->a2;
}
#ifdef __NVCC__
__host__ __device__
#endif
double* get_a3(Perturb_t* ptrb_ptr){
  return ptrb_ptr->a3;
}


double pol2pot(Config_t* cfg_ptr, double pdum){
  double potout;
  const double pw = get_pw(cfg_ptr->eqlb_ptr);
  const double pamp = get_pamp(cfg_ptr);
  const double rprof = get_rprof(cfg_ptr);
  potout = pamp * exp(-rprof * pdum / pw);
  return potout;
}


#ifdef __NVCC__
__host__ __device__
#endif
double* get_amp(Perturb_t* ptrb_ptr){
  return ptrb_ptr->amp;
}

#ifdef __NVCC__
__host__ __device__
#endif
int* get_mmod(Perturb_t* ptrb_ptr){
  return ptrb_ptr->mmod;
}

#ifdef __NVCC__
__host__ __device__
#endif
int* get_nmod(Perturb_t* ptrb_ptr){
  return ptrb_ptr->nmod;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_xi1(Perturb_t* ptrb_ptr){
  return ptrb_ptr->xi1;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_xi2(Perturb_t* ptrb_ptr){
  return ptrb_ptr->xi2;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_xi3(Perturb_t* ptrb_ptr){
  return ptrb_ptr->xi3;
}
#ifdef __NVCC__
__host__ __device__
#endif
void modestep(Config_t* cfg_ptr){
  /* interally this loops over particles already,
   so it is the invokers responsibility to ensure it is
  call only once. */
  /* modestep needs the whole kitchen */
  Equilib_t* eqlb_ptr = cfg_ptr->eqlb_ptr;
  Perturb_t* ptrb_ptr = cfg_ptr->ptrb_ptr;
  Particles_t* ptcl_ptr = cfg_ptr->ptcl_ptr;
  /* from struct */
  double *damdt = ptrb_ptr->damdt;
  double *dpzdt = ptrb_ptr->dpzdt;
  double *dp = ptrb_ptr->dp;
  double *dadm = ptrb_ptr->dadm;
  double *phaz = ptrb_ptr->phaz;
  double *amp = ptrb_ptr->amp;
  double *xi1 = ptrb_ptr->xi1;
  double *xi2 = ptrb_ptr->xi2;
  double *xi3 = ptrb_ptr->xi3;
  double *omegv = ptrb_ptr->omegv;
  int *nmod = ptrb_ptr->nmod;
  int *mmod = ptrb_ptr->mmod;
  int md1 = ptrb_ptr->md1;
  int md2 = ptrb_ptr->md2;
  double falf = ptrb_ptr->falf;

  /* from ptcl */
  double *zet = get_zet(ptcl_ptr);
  double *thet = get_thet(ptcl_ptr);
  double *g = get_g(ptcl_ptr);
  double *ri = get_ri(ptcl_ptr);
  double *q = get_q(ptcl_ptr);
  double *b = get_b(ptcl_ptr);
  double *rho = get_rho(ptcl_ptr);
  double* pol = get_pol(ptcl_ptr);
  double* wt = get_wt(ptcl_ptr);
  int idm = get_idm(ptcl_ptr);

  /* from eqlb */
  double pw = get_pw(eqlb_ptr);

  /* from  config */
  const double dt0 = get_dt0(cfg_ptr);
  int nprt = cfg_ptr->nprt;

  /* locals */
  int k;
  int ind;
  int md;
  int ncon;
  int jd;
  double agg;
  double ptav;
  double rbav;
  double pdum;
  double dpx;
  double dp2;
  double xinm;
  double alnm;
  double ptnm;
  double rbb;
  double dum;
  double dum1;
  double omg;
  double cnm;
  double snm;

  /* falf set from config */

  /* zero arrays */
  for(md=0; md < ptrb_ptr->modes; md++){
    damdt[md] = 0.;
    dpzdt[md] = 0.;
  }

  /* add stepping contributions, all particles and modes */
  for(md=md1; md<md2; md++){
    ncon = 0;
    ptav = 0.;
    rbav = 0.;

    for(k=0; k<nprt; k++){

      agg = nmod[md]*zet[k] - mmod[md]*thet[k];
      cnm = cos(agg - phaz[md*idm + k]);
      snm = sin(agg - phaz[md*idm + k]);
      pdum = pol[k];
      if(pdum < pw){
        ncon += 1;
        jd = compute_jd(eqlb_ptr, pdum);
        dpx = pdum - ((double)jd) * pw / (get_lsp(eqlb_ptr)-1);
        dp2 = dpx*dpx;
        ind = ptrb_ptr->lpt * md + k;
        xinm = amp[md]*(xi1[ind] + xi2[ind]*dpx + xi3[ind]*dp2);
        dum = 1./(mmod[md]*g[k] + nmod[md] * ri[k]);
        dum1 = mmod[md] - nmod[md] * q[k];
        alnm = xinm * dum1 * dum;
        ptnm = -(g[k]*q[k]+ri[k])*omegv[md]*xinm*dum;
        rbb = rho[k] * b[k] * b[k] *alnm - ptnm;
        rbav += agg;
        ptav += phaz[md*idm +k];
        dum = pdum / pw;
        /* XXX ind */
        damdt[md] += wt[ind] * rbb * cnm;  /* eq 6.80 */
        dpzdt[md] += wt[ind] * rbb * snm;  /* eq 6.81 */

        //dont think neededtd = nstep * dt0/tran;
      } /* pdum <pw  */
    }  /* k */
  }  /* md */

  for(md=md1; md<md2; md++){
    omg = omegv[md];
    dadm[md] = (falf*falf) * damdt[md] * dt0 / (
        ncon * omg * (amp[md]*amp[md]));
    dp[md] = -1. * (falf*falf) *dpzdt[md]*dt0/(
        ncon * omg * (amp[md]*amp[md]));
    amp[md] = fmax(amp[md], 1.E-9);
  }  /* md */

  for(md=md1; md<md2; md++){
    for(k=0; k<nprt; k++){
      phaz[md*idm + k] += dp[md];
    }  /* k */
  }    /* md */

  /* skipped "DATA OUT..." */
}
