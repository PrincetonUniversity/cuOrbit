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

#ifdef __NVCC__
#include <cuda.h>
#include <cuda_runtime.h>
#endif

#include "orbit_config_api.h"
#include "orbit_particles.h"
#include "orbit_util.h"
#include "orbit_constants.h"
#include "cuda_helpers.h"


struct Particles {
  int nprt;   /* number of particles (use in place of IDM*/
  double chrg;  /* 1 for ion, -1 for electron */
  double zprt;
  double prot;  /* proton mass in proton units */
  /* particle distribution */
  double ekev;
  double polo_scale;
  double p1_scale;
  double p2_scale;
  double pchi;

  size_t idm;
  int *otp;
  double *pol;
  double *zet;
  double *thet;
  double *rho;
  double *en;
  double *rmu;
  double *ptch;
  double *pot;
  double *time;  /* time step */
  double *g;
  double *gp;
  double *q;
  double *qp;
  double *b;  /* B, I associated with particle */
  double *ri;
  double *rip;
  double *w1, *w2, *w3;  /* particle stepping */
  double *dbdt;
  double *dbdp;
  double *dbdz;
  double *dbdpp;
  double *dbdpt;
  double *rpl;
  double *rplp;
  double *rplt;
  /* init during set1, not sure what used for yet*/
  double* dt;
  double* tim1;
  double* wt;
  /* used during kupdate */
  int* nout;
  double* e0;
  double* dptdp;
  double* dptdt;
  double* dptdz;
  double *alp;
  double *dadp;
  double *dadt;
  double *dadz;
  double *padt;

};

Particles_t* Particles_ctor(){
  return (Particles_t*)umacalloc(1, sizeof(Particles_t));
}


void initialize_Particles(Particles_t* ptcl_ptr, Config_t* cfg_ptr){

  /* from the config we can get number of particles
     which we will use to dimension our arrays
  */
  ptcl_ptr->nprt = cfg_ptr->nprt;
  ptcl_ptr->idm = (unsigned)ptcl_ptr->nprt;

  ptcl_ptr->chrg = cfg_ptr->chrg;
  ptcl_ptr->zprt = cfg_ptr->zprt;
  ptcl_ptr->prot = cfg_ptr->prot;
  ptcl_ptr->ekev = cfg_ptr->ekev;
  ptcl_ptr->polo_scale = cfg_ptr->polo_scale;
  ptcl_ptr->p1_scale = cfg_ptr->p1_scale;
  ptcl_ptr->p2_scale = cfg_ptr->p2_scale;
  ptcl_ptr->pchi = cfg_ptr->pchi;


  /* arrays */
  ptcl_ptr->otp = (int*)umacalloc(ptcl_ptr->idm, sizeof(int));

  ptcl_ptr->pol=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->zet=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->thet=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->rho=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->en=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->rmu=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->ptch=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->pot=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->time=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));  /* time step */
  ptcl_ptr->g=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->gp=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->q=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->qp=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->b=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));  /* B, I associated with particle */
  ptcl_ptr->ri=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->rip=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->w1=(double*)umacalloc(ptcl_ptr->idm, sizeof(double)); /* particle stepping */
  ptcl_ptr->w2=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->w3=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->dbdt=(double*)umacalloc(ptcl_ptr->idm, sizeof(double)); /* derivatives */
  ptcl_ptr->dbdp=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->dbdz=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->dbdpp=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->dbdpt=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));

  ptcl_ptr->rpl=(double*)umacalloc(ptcl_ptr->idm, sizeof(double)); /* derivatives */
  ptcl_ptr->rplp=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->rplt=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));

  /* stuff from set1 */
  ptcl_ptr->dt=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->tim1=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->wt=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));

  /* used during kupdate (stepping) */
  ptcl_ptr->nout=(int*)umacalloc(ptcl_ptr->idm, sizeof(int));
  ptcl_ptr->e0=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->dptdp=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->dptdt=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->dptdz=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->alp = (double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->dadp = (double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->dadt = (double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->dadz = (double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->padt = (double*)umacalloc(ptcl_ptr->idm, sizeof(double));


}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_b(Particles_t* ptcl_ptr){
  return ptcl_ptr->b;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_g(Particles_t* ptcl_ptr){
  return ptcl_ptr->g;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_q(Particles_t* ptcl_ptr){
  return ptcl_ptr->q;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_en(Particles_t* ptcl_ptr){
  return ptcl_ptr->en;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_ri(Particles_t* ptcl_ptr){
  return ptcl_ptr->ri;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_pol(Particles_t* ptcl_ptr){
  return ptcl_ptr->pol;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_rho(Particles_t* ptcl_ptr){
  return ptcl_ptr->rho;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_rmu(Particles_t* ptcl_ptr){
  return ptcl_ptr->rmu;
}

#ifdef __NVCC__
__host__ __device__
#endif
int* get_otp(Particles_t* ptcl_ptr){
  return ptcl_ptr->otp;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_ptch(Particles_t* ptcl_ptr){
  return ptcl_ptr->ptch;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_thet(Particles_t* ptcl_ptr){
  return ptcl_ptr->thet;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_pot(Particles_t* ptcl_ptr){
  return ptcl_ptr->pot;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_zet(Particles_t* ptcl_ptr){
  return ptcl_ptr->zet;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_time(Particles_t* ptcl_ptr){
  return ptcl_ptr->time;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_dt(Particles_t* ptcl_ptr){
  return ptcl_ptr->dt;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_tim1(Particles_t* ptcl_ptr){
  return ptcl_ptr->tim1;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_wt(Particles_t* ptcl_ptr){
  return ptcl_ptr->wt;
}


#ifdef __NVCC__
__host__ __device__
#endif
double get_zprt(Particles_t* ptcl_ptr){
  return ptcl_ptr->zprt;
}

#ifdef __NVCC__
__host__ __device__
#endif
double get_prot(Particles_t* ptcl_ptr){
  return ptcl_ptr->prot;
}

#ifdef __NVCC__
__host__ __device__
#endif
double get_ekev(Particles_t* ptcl_ptr){
  return ptcl_ptr->ekev;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_alp(Particles_t* ptcl_ptr){
  return ptcl_ptr->alp;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_dadp(Particles_t* ptcl_ptr){
  return ptcl_ptr->dadp;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_dadt(Particles_t* ptcl_ptr){
  return ptcl_ptr->dadt;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_dadz(Particles_t* ptcl_ptr){
  return ptcl_ptr->dadz;
}

#ifdef __NVCC__
__host__ __device__
#endif
double* get_padt(Particles_t* ptcl_ptr){
  return ptcl_ptr->padt;
}


void field(Config_t* cfg_ptr, const int nprts){
  int k;
  for(k=0; k < nprts; k++){
    kfield(cfg_ptr, k);
  }
  return;
}

#ifdef __NVCC__
__host__ __device__
#endif
void kfield(Config_t* cfg_ptr, int k){
  Particles_t* ptcl_ptr = cfg_ptr->ptcl_ptr;
  Equilib_t* const  Eq = cfg_ptr->eqlb_ptr;
  const int lsp = get_lsp(Eq);
  const double pw = get_pw(Eq);
  const double* thet = get_thet(ptcl_ptr);
  const double* zet = get_zet(ptcl_ptr);
  const int lst = get_lst(Eq);
  const int nrip = get_nrip(Eq);
  double* const rpl =   ptcl_ptr->rpl;
  double* const rplp =   ptcl_ptr->rplp;
  double* const rplt =   ptcl_ptr->rplt;
  double ** const R = get_R(Eq);
  double ** const B = get_B(Eq);
  double ** const QD = get_QD(Eq);
  double ** const GD = get_GD(Eq);
  double ** const RD = get_RD(Eq);
  double pdum, sdum, tdum;
  int idum, kd;
  double dpx, dp2;


  pdum = ptcl_ptr->pol[k];
  const int jd = compute_jd(Eq, pdum);
  dpx = pdum - ((double)jd) * pw / (lsp-1);
  dp2 = dpx*dpx;
  //  functions q and qp
  ptcl_ptr->q[k] = QD[0][jd] + QD[1][jd]*dpx + QD[2][jd]*dp2;
  ptcl_ptr->qp[k] = QD[1][jd] + 2*QD[2][jd]*dpx;
  ptcl_ptr->g[k] = GD[0][jd] + GD[1][jd]*dpx + GD[2][jd]*dp2;
  ptcl_ptr->gp[k] = GD[1][jd] + 2*GD[2][jd]*dpx;
  ptcl_ptr->ri[k] = RD[0][jd] + RD[1][jd]*dpx + RD[2][jd]*dp2;
  ptcl_ptr->rip[k] = RD[1][jd] + 2*RD[2][jd]*dpx;

  idum = thet[k] * pi2i;
  tdum = thet[k] - pi2*(idum-1);

  idum = tdum * pi2i;
  tdum = tdum - pi2*idum;

  kd = tdum*lst*pi2i;
  kd = imax(0,kd);
  kd = imin(lst-1,kd);
  const double dtx = tdum - kd*pi2/lst;
  const double dt2 = dtx*dtx;
  const int ind = jd*lst + kd;
  sdum = .5 - .5 * copysign(1. , jd - 0.5 );
  dpx = (1. - sdum) * dpx + sdum * sqrt( fmax(1.E-20, dpx));
  dp2 = dpx*dpx;

  if (get_krip(cfg_ptr->eqlb_ptr) != 0){
    rpl[k] = R[0][ind] + R[1][ind]*dpx + R[2][ind]*dp2
        + R[3][ind]*dtx + R[4][ind]*dpx*dtx + R[5][ind]*dtx*dp2
        + R[6][ind]*dt2 + R[7][ind]*dt2*dpx + R[8][ind]*dt2*dp2;

    rpl[k] *= .5 *(1. + copysign(1., pw - ptcl_ptr->pol[k]));

    rpl[k] = exp(rpl[k]);

    rplp[k] = R[1][ind] + 2*R[2][ind]*dpx + R[4][ind]*dtx
        + 2*R[5][ind]*dtx*dpx + R[7][ind]*dt2 + 2*R[8][ind]*dpx*dt2;

    rplp[k] *= rpl[k];

    rplt[k] = R[3][ind] + R[4][ind]*dpx + R[5][ind]*dp2
        + 2*R[6][ind]*dtx + 2*R[7][ind]*dtx*dpx + 2*R[8][ind]*dtx*dp2;

    rplt[k] *= rpl[k];

    ptcl_ptr->b[k] = B[0][ind] + B[1][ind]*dpx + B[2][ind]*dp2
        + B[3][ind]*dtx + B[4][ind]*dpx*dtx + B[5][ind]*dtx*dp2
        + B[6][ind]*dt2 + B[7][ind]*dt2*dpx + B[8][ind]*dt2*dp2
        + rpl[k] * sin(nrip * zet[k]);

    ptcl_ptr->dbdp[k] = B[1][ind] + 2*B[2][ind]*dpx + B[4][ind]*dtx
        + 2*B[5][ind]*dtx*dpx + B[7][ind]*dt2 + 2*B[8][ind]*dpx*dt2;

    ptcl_ptr->dbdp[k] /= (1. - sdum + 2 * sdum * dpx);

    ptcl_ptr->dbdt[k] = B[3][ind] + B[4][ind]*dpx + B[5][ind]*dp2
        + 2*B[6][ind]*dtx + 2*B[7][ind]*dtx*dpx + 2*B[8][ind] *dtx*dp2
        + rplt[k] * sin(nrip * zet[k]);

    ptcl_ptr->dbdz[k] = nrip * rpl[k] * cos(nrip * zet[k]);

    ptcl_ptr->dbdp[k] = ptcl_ptr->dbdp[k] + rplp[k] * sin(nrip * zet[k]);
  } else {
    /* krip is 0 */
    ptcl_ptr->b[k] = B[0][ind] + B[1][ind]*dpx + B[2][ind]*dp2
        + B[3][ind]*dtx + B[4][ind]*dpx*dtx + B[5][ind]*dtx*dp2
        + B[6][ind]*dt2 + B[7][ind]*dt2*dpx + B[8][ind]*dt2*dp2;

    ptcl_ptr->dbdp[k] = B[1][ind] + 2*B[2][ind]*dpx + B[4][ind]*dtx
        + 2*B[5][ind]*dtx*dpx + B[7][ind]*dt2 + 2*B[8][ind]*dpx*dt2;
    ptcl_ptr->dbdp[k] /= (1. - sdum + 2 * sdum * dpx);

    ptcl_ptr->dbdt[k] = B[3][ind] + B[4][ind]*dpx + B[5][ind]*dp2
        + 2*B[6][ind]*dtx + 2*B[7][ind]*dtx*dpx + 2*B[8][ind]*dtx*dp2;

    ptcl_ptr->dbdpp[k] = 2*B[2][ind] + 2*B[5][ind]*dtx + 2*B[8][ind]*dt2;
    ptcl_ptr->dbdpp[k] /= (1. - sdum + 2*sdum*dpx);

    ptcl_ptr->dbdpt[k] = B[4][ind] + 2*B[5][ind]*dpx + 2*B[7][ind]*dtx + 4*B[8][ind]*dpx*dtx;
    ptcl_ptr->dbdpt[k] /= (1. - sdum + 2*sdum*dpx);
  }

  if( get_npert(cfg_ptr) > 0) {
    ptrbak(cfg_ptr, k);
    ptrb2k(cfg_ptr, k);
  }

  return;

}

#ifdef __NVCC__
__host__ __device__
#endif
void kupdate(Config_t* cfg_ptr, int k){
   //locals
  int nesc;
  int md;
  Particles_t* Ptcl = cfg_ptr->ptcl_ptr;
  Perturb_t* Ptrb = cfg_ptr->ptrb_ptr;
  double* time = get_time(Ptcl);
  double* tim1 = get_tim1(Ptcl);
  const double pw = get_pw(cfg_ptr->eqlb_ptr);
  double* b = Ptcl->b;
  double* rmu = Ptcl->rmu;
  double* dt = Ptcl->dt;
  double* w1 = Ptcl->w1;
  double* w3 = Ptcl->w3;
  double* en = Ptcl->en;
  double* rho = Ptcl->rho;
  double* pot = Ptcl->pot;
  double* ptch = Ptcl->ptch;
  int* nout = Ptcl->nout;
  double* e0 = Ptcl->e0;
  double* phaz = get_phaz(Ptrb);
  double* omegv = get_omegv(Ptrb);
  const int md1 = get_md1(Ptrb);
  const int md2 = get_md2(Ptrb);
  const int idm = Ptcl->idm;

  kfield(cfg_ptr, k);

  nesc = 0;

  nout[k] = (int)(.6 *(1.  + copysign(1. , Ptcl->pol[k] - pw)));

  time[k] += dt[k]*(1 - nout[k]);
  tim1[k] += (1-nout[k])*dt[k];

  //  monitor the energy change
  //     and update the energy and pitch if not lost

  w3[k] = pot[k] + .5 *pow((b[k]*rho[k]),2) + rmu[k]*b[k];  /* def one of these */
  w1[k] = fabs((w3[k] - en[k])/e0[k]);
  en[k] = w3[k]*(1 - nout[k]) + en[k]*nout[k];
  ptch[k] *= nout[k] + (1  - nout[k])*rho[k]*b[k]/sqrt(2. * en[k] - 2. * pot[k]);
  nesc += nout[k];
  for(md=md1; md<md2; md++){// zero inds...
    phaz[md*idm + k] +=  omegv[md]*dt[k];
  }

  return;
}

#ifdef __NVCC__
__host__ __device__
#endif
void ptrb2k(Config_t* cfg_ptr, int k)
{
  double* dptdp = cfg_ptr->ptcl_ptr->dptdp;
  double* pot = cfg_ptr->ptcl_ptr->pot;
  double* pol = cfg_ptr->ptcl_ptr->pol;
  const double pamp = get_pamp(cfg_ptr);
  const double rprof = get_rprof(cfg_ptr);
  const double engn = get_engn(cfg_ptr);
  const double ekev = get_ekev(cfg_ptr->ptcl_ptr);
  const double pw = get_pw(cfg_ptr->eqlb_ptr);
  double dum;

  dum = pamp*engn/ekev;
  pot[k] += dum*exp(-rprof*pol[k]/pw);
  dptdp[k] -= rprof/pw*dum*exp(-rprof*pol[k]/pw);

  return;
}

#ifdef __NVCC__
__host__ __device__
#endif
void ptrbak(Config_t* cfg_ptr, int k)
{
  Equilib_t* Eqlb = cfg_ptr->eqlb_ptr;
  Particles_t* Ptcl = cfg_ptr->ptcl_ptr;
  Perturb_t* Ptrb = cfg_ptr->ptrb_ptr;
  double* alp = get_alp(Ptcl);
  double* dadp = get_dadp(Ptcl);
  double* dadt = get_dadt(Ptcl);
  double* dadz = get_dadz(Ptcl);
  double* padt = get_padt(Ptcl);
  double* pot = Ptcl->pot;
  double* dptdt = Ptcl->dptdt;
  double* dptdz = Ptcl->dptdz;
  double* dptdp = Ptcl->dptdp;
  int nflr = get_nflr(Ptrb);
  double* en = Ptcl->en;
  double* pol = Ptcl->pol;
  const double zprt = Ptcl->zprt;
  const double engn = get_engn(cfg_ptr);
  const double ekev = get_ekev(Ptcl);
  const double pw = get_pw(Eqlb);
  const int md1 = get_md1(Ptrb);
  const int md2 = get_md2(Ptrb);
  double* ptch = Ptcl->ptch;
  double prot = Ptcl->prot;
  double* amp = get_amp(Ptrb);
  double* a1 = get_a1(Ptrb);
  double* a2 = get_a2(Ptrb);
  double* a3 = get_a3(Ptrb);
  double* xi1 = get_xi1(Ptrb);
  double* xi2 = get_xi2(Ptrb);
  double* xi3 = get_xi3(Ptrb);
  int* nmod = get_nmod(Ptrb);
  int* mmod = get_mmod(Ptrb);
  const double* thet = get_thet(Ptcl);
  const double* zet = get_zet(Ptcl);
  const int idm = Ptcl->idm;
  double* b = Ptcl->b;
  double* g = Ptcl->g;
  double* gp = Ptcl->gp;
  double* q = Ptcl->q;
  double* qp = Ptcl->qp;
  double* ri = Ptcl->ri;
  double* rip = Ptcl->rip;
  const double bkg = get_bkg(cfg_ptr);
  double* phaz = get_phaz(Ptrb);
  double* omegv = get_omegv(Ptrb);
  int lpt = get_lpt(Ptrb);


  /* locals */
  int md,n,m,jd;
  int kf,nflr0;
  double pdum,dpx,dp2,agg;
  double alnm,alnmp,cnm,snm;
  double xinm,xinmp,gqi,gqip,gmni,gmnip;
  double edum,xisnm,xicnm,xipsnm;
  double x0p,x1p,rhol;
  double psi_flr,dpsi_flr,ph_flr,dph_flr;
  int ind;
  int lptm1;

  /* to appease compiler */
  dpsi_flr = 0.;
  dph_flr = 0.;

  alp[k] = 0.;
  dadp[k] = 0.;
  dadt[k] = 0.;
  dadz[k] = 0.;
  padt[k] = 0.;
  pot[k] = 0.;
  dptdt[k] = 0.;
  dptdz[k] = 0.;
  dptdp[k] = 0.;

  // make sure energy is finite
  nflr0=nflr;
  edum=en[k]*ekev/engn;
  if((edum-edum) != 0){
    nflr0=1 ;// skip FLR corrections if no finite energy
  }

  // compute GC particle position, Larmor radius
  if(nflr0 > 1){
    x0p = rpol(Eqlb, pol[k]); // gives rminor [m] from Psi_pol
    dph_flr = M_PI / nflr;
    rhol=1.02 * sqrt(prot) / zprt * (1.-ptch[k] * ptch[k]) *
      sqrt(1.E3 * en[k] * ekev / engn) / (1.E3 * bkg * b[k]); // [m] from NRL Plasma Formulary

    // look for max displacement perpendicular to flux surface
    x1p = x0p+rhol;
    psi_flr=polr_mp(Eqlb,  x1p, pol[k]+0.05);
    // get DPsi
    dpsi_flr = fabs(psi_flr-pol[k]);
  }

  //do  312 md = md1,md2 ! loop over harmonics
  for(md=md1; md<md2; md++){
    n = nmod[md];
    m = mmod[md];

    //! reset variables
    xisnm = 0.;
    xicnm = 0.;
    xipsnm = 0.;

    //do kf=1,nflr0 ! loop over gyro-orbit, average perturbation
    for(kf=0; kf<nflr0; kf++){
      if(nflr0 == 1){
        pdum = pol[k];
      }
      if(nflr0 > 1) {
        ph_flr= kf * dph_flr;
        pdum = pol[k] + dpsi_flr * cos(ph_flr);
      }

      lptm1 = lpt - 1;
      jd = pdum * lptm1 / pw;
      jd = imin(jd,lptm1);
      jd = imax(jd,0);
      dpx = pdum - jd *pw / lptm1;
      dp2 = dpx*dpx;

      ind = md*lpt + jd;

      alnm = amp[md]*(a1[ind] + a2[ind]*dpx + a3[ind]*dp2);
      alnmp = amp[md]*( a2[ind] + 2*a3[ind]*dpx);
      xinm = amp[md]*(xi1[ind] + xi2[ind]*dpx + xi3[ind]*dp2);
      xinmp = amp[md]*( xi2[ind] + 2*xi3[ind]*dpx);
      agg = n*zet[k] - m*thet[k];
      cnm = cos(agg - phaz[md*idm + k]);
      snm = sin(agg - phaz[md*idm + k]);
      alp[k] += alnm*snm/nflr0;
      dadp[k] += alnmp*snm/nflr0;
      dadt[k] -= m*alnm*cnm/nflr0;
      dadz[k] += alnm*n*cnm/nflr0;
      padt[k] -= omegv[md]*alnm*cnm/nflr0;

      xisnm += xinm*snm/nflr0;
      xicnm += xinm*cnm/nflr0;
      xipsnm += xinmp*snm/nflr0;

    }  //end loop - mock-up FLR effects

    // - compute the potential for each particle to make E parallel zero
    gqi = g[k]*q[k] + ri[k];
    gqip = gp[k]*q[k] + g[k]*qp[k] + rip[k];
    gmni = m*g[k] + n*ri[k];
    gmnip = m*gp[k] + n*rip[k];
    pot[k] -= xisnm*gqi*omegv[md]/gmni;
    dptdt[k] += m*xicnm*gqi*omegv[md]/gmni;
    dptdz[k] -= n*xicnm*gqi*omegv[md]/gmni;
    dptdp[k] = dptdp[k] - xipsnm*gqi*omegv[md]/gmni
      - xisnm*gqip*omegv[md]/gmni
      + xisnm*gqi*omegv[md]*gmnip/(gmni*gmni);
  } // loop over harmonics


  return;
}

/* For do_particles, we'll define different launch code for cpu and gpu,
 but use the same exact kernel code */

#ifdef __NVCC__
__host__ __device__
#endif
void do_particle_kernel(Config_t* cfg_ptr, int particle_id){
  int ktm;
  const int pdedp_tskip = get_pdedp_tskip(cfg_ptr->depo_ptr);

  for(ktm=1; ktm <= get_nstep_all(cfg_ptr); ktm++){
    konestep(cfg_ptr, particle_id);
    kupdate(cfg_ptr, particle_id);

    if(cfg_ptr->do_modestep && particle_id !=0){
      modestep(cfg_ptr);
    }

    if(compute_pdedp(cfg_ptr->depo_ptr) &&
       ktm >= pdedp_tskip &&
       ktm % pdedp_tskip == 0){
      rcrd_vararr(cfg_ptr, particle_id, ktm/pdedp_tskip - 1);
    }
  }
}


#ifdef __NVCC__
__global__
void do_particles_dev(Config_t* cfg_ptr){

  /* 1D grid of 1D blocks */
  const int particle_id = blockIdx.x * blockDim.x + threadIdx.x;
  if(particle_id >= cfg_ptr->nprt){
    return;
  }

  do_particle_kernel(cfg_ptr, particle_id);
}
#endif

void do_particles_host(Config_t* cfg_ptr){

  int particle_id;

  for(particle_id=0; particle_id < cfg_ptr->ptcl_ptr->nprt; particle_id++){
    do_particle_kernel(cfg_ptr, particle_id);
  }
}


void do_particles(Config_t* cfg_ptr){
#ifdef __NVCC__
  /* compute a reasonable launch size */
  int dimBlock=0;   /*  The launch configurator returned block size */
  int minGridSize=0; /*  The minimum grid size needed to achieve
                           the maximum occupancy for a full device launch */
  int dimGrid=0;    /*  The actual grid size needed, based on input size */

  int N = cfg_ptr->nprt;

  HANDLE_ERROR(cudaOccupancyMaxPotentialBlockSize(
      &minGridSize,
      &dimBlock,
      do_particles_dev,  /* kernel we're qrying */
      0,  /* smemsize, laffs, don't worry about it */
      0)); /* ubound block, 0 none */
  /* compute the grid size given recomended block */
  dimGrid = (N + dimBlock - 1) / dimBlock;
  printf("Launching do_particles as %d blocks of %d threads\n", dimGrid, dimBlock);

  /* launch it */
  do_particles_dev<<<(unsigned)dimGrid, (unsigned)dimBlock>>>(cfg_ptr);

  /* peek */
  HANDLE_ERROR(cudaPeekAtLastError());

  /* you need this for UVM... */
  HANDLE_ERROR(cudaDeviceSynchronize());

#else
  do_particles_host(cfg_ptr);
#endif

}


#ifdef __NVCC__
__host__ __device__
#endif
void konestep(Config_t* cfg_ptr, int k){

  Particles_t* Ptcl = cfg_ptr->ptcl_ptr;
  int* nout = Ptcl->nout;
  double* dt = Ptcl->dt;
  double* dptdp = cfg_ptr->ptcl_ptr->dptdp;
  double* pol = cfg_ptr->ptcl_ptr->pol;
  const double pw = get_pw(cfg_ptr->eqlb_ptr);
  int npert = get_npert(cfg_ptr);
  double* b = Ptcl->b;
  double* g = Ptcl->g;
  double* gp = Ptcl->gp;
  double* q = Ptcl->q;
  double* ri = Ptcl->ri;
  double* rip = Ptcl->rip;
  double* alp = Ptcl->alp;
  double* zet = Ptcl->zet;
  double* thet = Ptcl->thet;
  double* rho = Ptcl->rho;
  double* rmu = Ptcl->rmu;
  double* dadp = Ptcl->dadp;
  double* dbdp = Ptcl->dbdp;
  double* dadz = Ptcl->dadz;
  double* dbdz = Ptcl->dbdz;
  double* dadt = Ptcl->dadt;
  double* dbdt = Ptcl->dbdt;
  double* dptdt = Ptcl->dptdt;
  double* dptdz = Ptcl->dptdz;
  double* padt = Ptcl->padt;
  const double chrg = Ptcl->chrg;
  const double dt0 = get_dt0(cfg_ptr);

  int n1,j,i,ndum,ei;
  double xdum,ydum,rbb,dedb,deni,fac1,fac2,pdum,xdot,ydot;
  /* local temps, note compiler is suspicuous of those I have init 0..check */
  double y_[4];
  double d_[4];
  double e_[4];
  double a_[4]= {0.};
  double bx_[4] = {0.};
  double h_;
  int nt_;
  double c1_[4] = {0.};

  dt[k]=dt0;

  nout[k] = (int)(.6 * (1. + copysign(1., pol[k]- pw) ));
  nt_ = (int) (.6 * (1.  + copysign(1., pol[k]-.05 *pw)));
  dt[k] = nt_ * dt[k] + (1 - nt_) * dt0;

  xdum = sqrt(pol[k])*cos(thet[k]);
  ydum = sqrt(pol[k])*sin(thet[k]);

  y_[0] = nt_ * pol[k] + (1-nt_)*xdum;  /* already wrong here */
  y_[1] = nt_ * thet[k] + (1-nt_)*ydum;
  y_[2] = zet[k];
  y_[3] = rho[k];

  d_[0] = y_[0];
  d_[1] = y_[1];
  d_[2] = y_[2];
  d_[3] = y_[3];
  h_ = dt[k]/6.0 ;

  for(j=0; j < 4; j++){
    kfield(cfg_ptr, k);

    /* the arithmetic that follows here could use a lot of cleanup.
     i think simplification herewould help the optimizer... and the devs... */
    if(npert ==  0){
      //goto 61
      rbb = y_[3] * b[k]*b[k];  /*  rbb is off, y_[3] seems okay, when did b change */
      dedb = b[k] * y_[3]*y_[3] + rmu[k];
      deni = 1. / ( g[k]*q[k] + ri[k] + chrg * y_[3] * (g[k]*rip[k] - ri[k]*gp[k]));

      //   pol dot
      e_[0] = - chrg*g[k]*dedb*dbdt[k] * deni +
          chrg * ri[k] * dedb *dbdz[k] * deni;

      //   thet dot
      e_[1] = (chrg*dedb*dbdp[k]*g[k] +
                 rbb*(1. - chrg * gp[k] * y_[3])) * deni;

      //   zet dot
      e_[2] = (-chrg*dedb*dbdp[k]*ri[k] +
               rbb*(q[k] + chrg * y_[3] * rip[k]))*deni;

      //   rho dot   - term in  dbdz*dadp given by Mynick
      e_[3] = -dedb * (1. - chrg * y_[3] * gp[k] ) * dbdt[k] * deni -
          dedb * dbdz[k] * (q[k] + chrg * y_[3] * rip[k]) * deni;

      pdum = y_[0] * y_[0] + y_[1] * y_[1];
      xdot = .5 * y_[0] * e_[0] / pdum - y_[1] * e_[1];
      ydot = .5 *y_[1] * e_[0]/pdum + y_[0] *e_[1];
      e_[0] = nt_*e_[0] + (1-nt_)*xdot;
      e_[1] = nt_*e_[1] + (1-nt_)*ydot;

      for(ei=0; ei<4; ei++){
        e_[ei] *= (1-nout[k]);
      }

    } else{

      rbb = y_[3] * b[k]*b[k];
      dedb = b[k] * y_[3]*y_[3] + rmu[k];
      deni = 1. / ( g[k]*q[k] + ri[k] + (chrg * y_[3]+alp[k]) * (g[k]*rip[k] - ri[k]*gp[k]));
      fac1 = 1 - gp[k]*(chrg*y_[3] + alp[k]) - g[k]*dadp[k];
      fac2 = q[k] + rip[k]*(chrg*y_[3] + alp[k]) + ri[k]*dadp[k];

      //   pol dot
      e_[0] = (-ri[k]*rbb*dadz[k]
               - chrg*g[k]*dedb*dbdt[k]
               + g[k]*rbb*dadt[k]
               - chrg*g[k]*dptdt[k]
               + chrg*ri[k]*dptdz[k]
               + chrg*ri[k]*dedb*dbdz[k]) * deni;

      //   thet dot
      e_[1] = (chrg*dedb*dbdp[k]*g[k]+ rbb*fac1 + chrg*g[k]*dptdp[k])*deni;
      //   zet dot
      e_[2] = (-chrg*dedb*dbdp[k]*ri[k] + rbb*fac2 - chrg*ri[k]*dptdp[k])*deni;

      //   rho dot   - term in  dbdz*dadp given by Mynick
      e_[3] = ( -fac2 * (dedb*dbdz[k]+ dptdz[k])
                - fac1* (dedb*dbdt[k] + dptdt[k])
                + (dedb*dbdp[k]+dptdp[k])
                * (ri[k]*dadz[k] - g[k]*dadt[k])) *deni - chrg*padt[k] ;
      pdum = y_[0] * y_[0] + y_[1] * y_[1];
      xdot = .5 * y_[0] * e_[0] / pdum - y_[1] * e_[1];
      ydot = .5 *y_[1] * e_[0]/pdum + y_[0] *e_[1];
      e_[0] = nt_*e_[0] + (1-nt_)*xdot;
      e_[1] = nt_*e_[1] + (1-nt_)*ydot;
      for(ei=0; ei<4; ei++){
        e_[ei] *= (1-nout[k]);
      }
    }

    //goto 62, like 42, the answer to the universe and everything in it,  but twenty years older.
    n1=4;
    for(i=0; i<n1; i++){
      if(j==0){
        a_[i] = h_ * e_[i];
        y_[i] = d_[i] + 3 * a_[i];
      }
      if(j==1){
        bx_[i] = h_ * e_[i];
        y_[i] = d_[i] + 3 * bx_[i];
      }
      if(j==2){
        c1_[i] = h_ * e_[i];
        y_[i] = d_[i] + 6 * c1_[i];
      }
      if(j==3){
        y_[i] = d_[i] + a_[i] + 2 * bx_[i] + 2 * c1_[i] + h_ * e_[i];
      }
    }

    // 40
    ndum = (int)( .6 * (1 - copysign(1. , y_[0])));
    pol[k] = nt_*y_[0] + (1-nt_)*(y_[0]*y_[0] + y_[1]*y_[1]);
    thet[k] = nt_*y_[1] + (1-nt_)*(atan(y_[1]/y_[0]) + ndum * M_PI);
    zet[k] = y_[2];
    rho[k] = y_[3];
    nout[k] = .6 * (1. + copysign(1., pol[k]-pw));

  } //j
  return;
}

#ifdef __NVCC__
__host__ __device__
#endif
int get_idm(Particles_t* ptcl_ptr){
  return ptcl_ptr->idm;
}


void sampledep(Config_t* cfg_ptr){
  /* config */
  const double engn = get_engn(cfg_ptr);
  Particles_t* ptcl_ptr = cfg_ptr->ptcl_ptr;

  /* ptcl */
  double* b = ptcl_ptr->b;
  double* en = ptcl_ptr->en;
  double* rho = ptcl_ptr->rho;
  double* pot = ptcl_ptr->pot;
  double* pol = ptcl_ptr->pol;
  double* ptch = ptcl_ptr->ptch;
  double* zet = ptcl_ptr->zet;
  double* thet = ptcl_ptr->thet;
  double* rmu = ptcl_ptr->rmu;
  const double ekev = get_ekev(ptcl_ptr);


  /* eqlb */
  const double pw = get_pw(cfg_ptr->eqlb_ptr);

  /* locals */
  const char *mode = "r";
  FILE* ifp;
  int nheader;
  int nlines;
  int lineno;
  double rdum,zdum,ptdum,edum;
  double thetd, pold;
  int k;

  /* we need to open a file */
  ifp = fopen(cfg_ptr->fbmdata_file, mode);
  if (ifp == NULL) {
    fprintf(stderr, "Can't open input file %s!\n", cfg_ptr->fbmdata_file);
    exit(1);
  }
  printf("Parsing particle deposition file %s\n",  cfg_ptr->fbmdata_file);


  /* figure out how many lines of rz data */
  fscanf(ifp, "%d %d", &nheader, &nlines);
  printf("sampledep found %d header lines and %d data lines\n", nheader, nlines);

  /* sanity check */
  if(cfg_ptr->nprt > nlines){
    fprintf(stderr,
            "Error, you requested more particles(%d) than are in %s.\n" \
            "Decrease nprt <= %d and retry\n",
            cfg_ptr->nprt,
            cfg_ptr->fbmdata_file,
            nlines);
    exit(1);
  }


  /* init counters */
  k = 0;
  /* read in the data */
  /* note, already read first line */
  for(lineno=1; lineno<nlines; lineno++){
    /* skip the nheader additional header lines */
    if(lineno < nheader+1) fscanf(ifp, "%*[^\n]\n");

    /* stay in bounds */
    if(k >= ptcl_ptr->nprt){
      fprintf(stderr,
              "We have more particle data available than nprt." \
              "Warning, only reading first nprt (%d) particles from %s\n",
              cfg_ptr->nprt,
              cfg_ptr->fbmdata_file);
      break;
    }
    /* read data */
    fscanf(ifp, "%lf %lf %lf %lf ", &rdum, &zdum, &ptdum, &edum);

    if(edum * 1.E-3 < cfg_ptr->pdedp_Emin ) break;  /* limit energy */

    /* do some stuff, keepiong track of particles within limits */
    // stub, need to impliment
    //XXX gc_xz2pt(rdum, zdum, xc, rmaj, ped, pold, thetd, xhave, zhave, dist, ierr);
    fprintf(stderr, "ERROR: gc_xz2pt or equiv is not implimented yet!\n"); exit(1);

    if(thetd > M_PI) thetd -= pi2;

    if(pold >= pw) break;  /* outside */

    if(pold <= 1.E-6) break;  /* axis */

    /* looks like we have a keeper */
    thet[k] = thetd;
    zet[k] = 0.;
    pol[k] = pold;
    ptch[k] = ptdum;
    en[k] = edum * engn * 1.E-3 / ekev ;  /* kinetic */
    k += 1;

  }  /* line */

  /* presumably we then change nprt... */
  cfg_ptr->nprt = k;
  ptcl_ptr->nprt = k;

  /* step */
  field(cfg_ptr, ptcl_ptr->nprt);
  for(k=0; k<ptcl_ptr->nprt; k++){
    rho[k] = ptch[k] * sqrt(2. * en[k]) / b[k];
    rmu[k] = en[k] / b[k] - .5 * rho[k] * rho[k] * b[k];
    en[k] = en[k] + pot[k];
  }  /* k */

  fclose(ifp);
  printf("sampledep completed read of %d particles\n",  ptcl_ptr->nprt);

}
