#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "orbit_config_api.h"
#include "orbit_particles.h"
#include "orbit_util.h"

typedef struct Particles {
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
  double *rpl;
  double *rplp;
  double *rplt;


} Particles_t;

Particles_t* Particles_ctor(){
  return (Particles_t*)calloc(1, sizeof(Particles_t));
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
  ptcl_ptr->otp = (int*)calloc(ptcl_ptr->idm, sizeof(int));

  ptcl_ptr->pol=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->zet=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->thet=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->rho=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->en=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->rmu=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->ptch=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->pot=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->time=(double*)calloc(ptcl_ptr->idm, sizeof(double));  /* time step */
  ptcl_ptr->g=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->gp=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->q=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->qp=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->b=(double*)calloc(ptcl_ptr->idm, sizeof(double));  /* B, I associated with particle */
  ptcl_ptr->ri=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->rip=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->w1=(double*)calloc(ptcl_ptr->idm, sizeof(double)); /* particle stepping */
  ptcl_ptr->w2=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->w3=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->dbdt=(double*)calloc(ptcl_ptr->idm, sizeof(double)); /* derivatives */
  ptcl_ptr->dbdp=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->dbdz=(double*)calloc(ptcl_ptr->idm, sizeof(double));

  ptcl_ptr->rpl=(double*)calloc(ptcl_ptr->idm, sizeof(double)); /* derivatives */
  ptcl_ptr->rplp=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->rplt=(double*)calloc(ptcl_ptr->idm, sizeof(double));



}


double* get_q(Particles_t* ptcl_ptr){
  return ptcl_ptr->q;
}

double* get_pol(Particles_t* ptcl_ptr){
  return ptcl_ptr->pol;
}
double* get_thet(Particles_t* ptcl_ptr){
  return ptcl_ptr->thet;
}
double* get_zet(Particles_t* ptcl_ptr){
  return ptcl_ptr->zet;
}


double get_zprt(Particles_t* ptcl_ptr){
  return ptcl_ptr->zprt;
}
double get_prot(Particles_t* ptcl_ptr){
  return ptcl_ptr->prot;
}
double get_ekev(Particles_t* ptcl_ptr){
  return ptcl_ptr->ekev;
}


void field(Config_t* cfg_ptr, Particles_t* ptcl_ptr, const int nprts){
  int k;
  for(k=0; k < nprts; k++){
    kfield(cfg_ptr, ptcl_ptr, k);
  }
  return;
}

void kfield(Config_t* cfg_ptr, Particles_t* ptcl_ptr, int k){
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
  idum = tdum *pi2i;
  tdum = tdum - pi2*idum;

  kd = tdum*lst*pi2i;
  kd = imax(0,kd);
  kd = imin(lst-1,kd);
  const double dtx = tdum - (kd-1)*pi2/lst;
  const double dt2 = dtx*dtx;
  const int ind = jd*lst + kd;
  rpl[k] = R[0][ind] + R[1][ind]*dpx + R[2][ind]*dp2
      + R[3][ind]*dtx + R[4][ind]*dpx*dtx + R[5][ind]*dtx*dp2
      + R[6][ind]*dt2 + R[7][ind]*dt2*dpx + R[8][ind]*dt2*dp2;
  rpl[k] *= .5 *(1. + copysign(1., pw- ptcl_ptr->pol[k]));
  rpl[k] = exp(rpl[k]);
  rplp[k] = R[1][ind] + 2*R[2][ind]*dpx + R[4][ind]*dtx
      + 2*R[5][ind]*dtx*dpx + R[7][ind]*dt2 + 2*R[8][ind]*dpx*dt2;
  rplp[k] = rplp[k]*rpl[k];
  rplt[k] = R[3][ind] + R[4][ind]*dpx + R[5][ind]*dp2
      + 2*R[6][ind]*dtx + 2*R[7][ind]*dtx*dpx + 2*R[8][ind]*dtx*dp2;
  rplt[k] = rplt[k]*rpl[k];

  sdum = .5 - .5 * copysign(1. , jd - 1.5 ); //xxx this is wrong
  dpx = (1. - sdum) * dpx + sdum * sqrt( fmax(1.E-20, dpx));
  dp2 = dpx*dpx;
  ptcl_ptr->b[k] = B[0][ind] + B[1][ind]*dpx + B[2][ind]*dp2
      + B[3][ind]*dtx + B[4][ind]*dpx*dtx + B[5][ind]*dtx*dp2
      + B[4][ind]*dt2 + B[7][ind]*dt2*dpx + B[8][ind]*dt2*dp2
      + rpl[k] * sin(nrip * zet[k]);
  ptcl_ptr->dbdp[k] = B[1][ind] + 2*B[2][ind]*dpx + B[4][ind]*dtx
      + 2*B[5][ind]*dtx*dpx + B[7][ind]*dt2 + 2*B[8][ind]*dpx*dt2;
  ptcl_ptr->dbdp[k] = ptcl_ptr->dbdp[k]/(1. - sdum + 2 * sdum * dpx);
  ptcl_ptr->dbdt[k] = B[3][ind] + B[4][ind]*dpx + B[5][ind]*dp2
      + 2*B[6][ind]*dtx + 2*B[7][ind]*dtx*dpx + 2*B[8][ind] *dtx*dp2
      + rplt[k] * sin(nrip * zet[k]);
  ptcl_ptr->dbdz[k] = nrip * rpl[k] * cos(nrip * zet[k]);
  ptcl_ptr->dbdp[k] = ptcl_ptr->dbdp[k] + rplp[k] * sin(nrip * zet[k]);
  return;
}
