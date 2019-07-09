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
  double* nout;
  double* nfin;
  double* e0;

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
  ptcl_ptr->dbdpp=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->dbdpt=(double*)calloc(ptcl_ptr->idm, sizeof(double));

  ptcl_ptr->rpl=(double*)calloc(ptcl_ptr->idm, sizeof(double)); /* derivatives */
  ptcl_ptr->rplp=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->rplt=(double*)calloc(ptcl_ptr->idm, sizeof(double));

  /* stuff from set1 */
  ptcl_ptr->dt=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->tim1=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->wt=(double*)calloc(ptcl_ptr->idm, sizeof(double));

  /* used during kupdate (stepping) */
  ptcl_ptr->nout=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->nfin=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->e0=(double*)calloc(ptcl_ptr->idm, sizeof(double));

}

double* get_b(Particles_t* ptcl_ptr){
  return ptcl_ptr->b;
}

double* get_g(Particles_t* ptcl_ptr){
  return ptcl_ptr->g;
}

double* get_q(Particles_t* ptcl_ptr){
  return ptcl_ptr->q;
}

double* get_en(Particles_t* ptcl_ptr){
  return ptcl_ptr->en;
}

double* get_pol(Particles_t* ptcl_ptr){
  return ptcl_ptr->pol;
}
double* get_rho(Particles_t* ptcl_ptr){
  return ptcl_ptr->rho;
}
double* get_rmu(Particles_t* ptcl_ptr){
  return ptcl_ptr->rho;
}
int* get_otp(Particles_t* ptcl_ptr){
  return ptcl_ptr->otp;
}
double* get_ptch(Particles_t* ptcl_ptr){
  return ptcl_ptr->ptch;
}
double* get_thet(Particles_t* ptcl_ptr){
  return ptcl_ptr->thet;
}
double* get_pot(Particles_t* ptcl_ptr){
  return ptcl_ptr->pot;
}
double* get_zet(Particles_t* ptcl_ptr){
  return ptcl_ptr->zet;
}
double* get_time(Particles_t* ptcl_ptr){
  return ptcl_ptr->time;
}
double* get_dt(Particles_t* ptcl_ptr){
  return ptcl_ptr->dt;
}
double* get_tim1(Particles_t* ptcl_ptr){
  return ptcl_ptr->tim1;
}
double* get_wt(Particles_t* ptcl_ptr){
  return ptcl_ptr->wt;
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


void field(Config_t* cfg_ptr, const int nprts){
  int k;
  for(k=0; k < nprts; k++){
    kfield(cfg_ptr, k);
  }
  return;
}

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
  idum = tdum *pi2i;
  tdum = tdum - pi2*idum;

  kd = tdum*lst*pi2i;
  kd = imax(0,kd);
  kd = imin(lst-1,kd);
  const double dtx = tdum - (kd-1)*pi2/lst;
  const double dt2 = dtx*dtx;
  const int ind = jd*lst + kd;
  sdum = .5 - .5 * copysign(1. , jd - 1.5 ); //xxx this is wrong
  dpx = (1. - sdum) * dpx + sdum * sqrt( fmax(1.E-20, dpx));
  dp2 = dpx*dpx;

  if (get_krip(cfg_ptr->eqlb_ptr) != 0){
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
  } else {
    /* krip is 0 */
    ptcl_ptr->b[k] = B[0][ind] + B[1][ind]*dpx + B[2][ind]*dp2
        + B[3][ind]*dtx + B[4][ind]*dpx*dtx + B[5][ind]*dtx*dp2
        + B[4][ind]*dt2 + B[7][ind]*dt2*dpx + B[8][ind]*dt2*dp2;

    ptcl_ptr->dbdp[k] = B[1][ind] + 2*B[2][ind]*dpx + B[4][ind]*dtx
        + 2*B[5][ind]*dtx*dpx + B[7][ind]*dt2 + 2*B[8][ind]*dpx*dt2;

    ptcl_ptr->dbdp[k] = ptcl_ptr->dbdp[k]/(1. - sdum + 2 * sdum * dpx);

    ptcl_ptr->dbdt[k] = B[3][ind] + B[4][ind]*dpx + B[5][ind]*dp2
        + 2*B[6][ind]*dtx + 2*B[7][ind]*dtx*dpx + 2*B[8][ind] *dtx*dp2;

    ptcl_ptr->dbdt[k] = B[3][ind] + B[4][ind]*dpx + B[5][ind]*dp2
      + 2*B[6][ind]*dtx + 2*B[7][ind]*dtx*dpx + 2*B[8][ind]*dtx*dp2;

    ptcl_ptr->dbdpp[k] = 2*B[2][ind] + 2*B[5][ind]*dtx + 2*B[8][ind]*dt2;

    ptcl_ptr->dbdpp[k] = ptcl_ptr->dbdpp[k]/(1. - sdum + 2*sdum*dpx);

    ptcl_ptr->dbdpt[k] = B[4][ind] + 2*B[5][ind]*dpx + 2*B[7][ind]*dtx + 4*B[8][ind]*dpx*dtx;

    ptcl_ptr->dbdpt[k] = ptcl_ptr->dbdpt[k]/(1. - sdum + 2*sdum*dpx);
  }

  if( get_npert(cfg_ptr->ptrb_ptr) > 0) {
    /* ptrbak() */
    /* ptrb2k() */
  }

  return;

}

void kupdate(Config_t* cfg_ptr, int k){
   //locals
  int nesc;
  int ntim,md;
  Particles_t* Ptcl = cfg_ptr->ptcl_ptr;
  Perturb_t* Ptrb = cfg_ptr->ptrb_ptr;
  double* time = get_time(Ptcl);
  double* tim1 = get_tim1(Ptcl);
  const double pw = get_pw(cfg_ptr->eqlb_ptr);
  double* b = Ptcl->b;
  double* rmu = Ptcl->b;
  double* dt = Ptcl->dt;
  double* w1 = Ptcl->w1;
  double* w3 = Ptcl->w3;
  double* en = Ptcl->en;
  double* rho = Ptcl->rho;
  double* pot = Ptcl->pot;
  double* ptch = Ptcl->ptch;
  double* nout = Ptcl->nout;
  double* nfin = Ptcl->nfin;
  double* e0 = Ptcl->e0;
  double* phaz = get_phaz(Ptrb);
  double* omegv = get_omegv(Ptrb);
  const int md1 = get_md1(Ptrb);
  const int md2 = get_md2(Ptrb);
  const int idm = Ptcl->idm;
  double trun = get_trun(cfg_ptr);

  kfield(cfg_ptr, k);

  nesc = 0;
  ntim = 0;

  nout[k] = .6 *(1.  + copysign(1. , Ptcl->pol[k] - pw));
  nfin[k] =  .6 *(1.  + copysign(1. , time[k]-trun));

  time[k] = time[k] + dt[k]*(1 - nout[k])*(1-nfin[k]);
  tim1[k] = tim1[k] + (1-nout[k])*dt[k]*(1-nfin[k]);

  //  monitor the energy change
  //     and update the energy and pitch if not lost
  w3[k] = pot[k] + .5 *pow((b[k]*rho[k]),2) + rmu[k]*b[k];
  w1[k] = fabs((w3[k] - en[k])/e0[k]);
  en[k] = w3[k]*(1 - nout[k]) + en[k]*nout[k];
  ptch[k] = ptch[k] * nout[k] + (1  - nout[k])*rho[k]*b[k]/sqrt(2. * en[k] - 2. * pot[k]);
  nesc = nesc + nout[k];
  ntim = ntim + nfin[k];
  for(md=md1-1; md<md2; md++){// zero inds...
    phaz[md*idm + k] = phaz[md*idm + k] + omegv[md]*dt[k];
  }

  return;
}
