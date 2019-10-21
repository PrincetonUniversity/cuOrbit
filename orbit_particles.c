#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "orbit_config_api.h"
#include "orbit_particles.h"
#include "orbit_util.h"
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
  double* nout;
  double* nfin;
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
  ptcl_ptr->nout=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->nfin=(double*)umacalloc(ptcl_ptr->idm, sizeof(double));
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

double* get_alp(Particles_t* ptcl_ptr){
  return ptcl_ptr->alp;
}

double* get_dadp(Particles_t* ptcl_ptr){
  return ptcl_ptr->dadp;
}

double* get_dadt(Particles_t* ptcl_ptr){
  return ptcl_ptr->dadt;
}

double* get_dadz(Particles_t* ptcl_ptr){
  return ptcl_ptr->dadz;
}

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
    ptrbak(cfg_ptr, k);
    ptrb2k(cfg_ptr, k);
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
  double* rmu = Ptcl->rmu;
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
  for(md=md1; md<md2; md++){// zero inds...
    phaz[md*idm + k] = phaz[md*idm + k] + omegv[md]*dt[k];
  }

  return;
}

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
  pot[k] =  pot[k] + dum*exp(-rprof*pol[k]/pw);
  dptdp[k] = dptdp[k] - rprof/pw*dum*exp(-rprof*pol[k]/pw);

  return;
}

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
  const double pamp = get_pamp(cfg_ptr);
  const double rprof = get_rprof(cfg_ptr);
  const double engn = get_engn(cfg_ptr);
  const double ekev = get_ekev(Ptcl);
  const double pw = get_pw(Eqlb);
  const int md1 = get_md1(Ptrb);
  const int md2 = get_md2(Ptrb);
  double* ptch = Ptcl->ptch;
  double prot = Ptcl->prot;
  double* alfv = get_alfv(Ptrb);
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
  int nval;
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
    x0p=rpol(Eqlb, pol[k]); // gives rminor [m] from Psi_pol
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
    nval = alfv[md];  /* xxx does this need to be saved in a struct ?*/

    n = nmod[md];
    m = mmod[md];

    //! reset variables
    xisnm = 0.;
    xicnm = 0.;
    xipsnm = 0.;

    //do kf=1,nflr0 ! loop over gyro-orbit, average perturbation
    for(kf=1; kf<=nflr0; kf++){
      if(nflr0 == 1){
        pdum = pol[k];
      }
      if(nflr0 > 1) {
        ph_flr=(kf-1.)*dph_flr;
        pdum=pol[k]+dpsi_flr*cos(ph_flr);
      }

      lptm1 = (lpt) - 1;
      jd = pdum*(lptm1)/pw + 1;
      jd = imin(jd,lptm1);
      jd = imax(jd,1);
      dpx = pdum - (jd-1)*pw/(lptm1);
      dp2 = dpx*dpx;
      // zero inds
      jd--;
      ind = md*910 + jd;
      alnm = amp[nval]*(a1[ind] + a2[ind]*dpx + a3[ind]*dp2);
      alnmp = amp[nval]*( a2[ind] + 2*a3[ind]*dpx);
      xinm = amp[nval]*(xi1[ind] + xi2[ind]*dpx + xi3[ind]*dp2);
      xinmp = amp[nval]*( xi2[ind] + 2*xi3[ind]*dpx);
      agg = n*zet[k] - m*thet[k];
      printf("DBG md %d , idm %d , k %d , phaz[md*idm + k]\n",md, idm, k);
      cnm = cos(agg - phaz[md*idm + k]);
      snm = sin(agg - phaz[md*idm + k]);
      alp[k] = alp[k] + alnm*snm/nflr0;
      dadp[k] = dadp[k] + alnmp*snm/nflr0;
      dadt[k] = dadt[k] - m*alnm*cnm/nflr0;
      dadz[k] = dadz[k] + alnm*n*cnm/nflr0;
      padt[k] = padt[k] - omegv[md]*alnm*cnm/nflr0;

      xisnm = xisnm + xinm*snm/nflr0;
      xicnm = xicnm + xinm*cnm/nflr0;
      xipsnm = xipsnm + xinmp*snm/nflr0;

    }  //end loop - mock-up FLR effects

    // - compute the potential for each particle to make E parallel zero
    gqi = g[k]*q[k] + ri[k];
    gqip = gp[k]*q[k] + g[k]*qp[k] + rip[k];
    gmni = m*g[k] + n*ri[k];
    gmnip = m*gp[k] + n*rip[k];
    pot[k]= pot[k] - xisnm*gqi*omegv[md]/gmni;
    dptdt[k] = dptdt[k] + m*xicnm*gqi*omegv[md]/gmni;
    dptdz[k] = dptdz[k] - n*xicnm*gqi*omegv[md]/gmni;
    dptdp[k] = dptdp[k] - xipsnm*gqi*omegv[md]/gmni
      - xisnm*gqip*omegv[md]/gmni
      + xisnm*gqi*omegv[md]*gmnip/(gmni*gmni);
  } // loop over harmonics

  return;
}


void do_particles(Config_t* cfg_ptr){
    // before we begin, we need to compute "eps"
  const double eps = compute_eps(cfg_ptr->eqlb_ptr);

  int particle_id;
  int ktm;
  for(particle_id=0; particle_id < cfg_ptr->ptcl_ptr->nprt; particle_id++){
    for(ktm=1; ktm < get_nstep_all(cfg_ptr); ktm++){
      /* printf("DBG particle id %d ktm %d\n", particle_id, ktm); */

      konestep(cfg_ptr, particle_id);

      kupdate(cfg_ptr, particle_id);

    }
  }
}

void konestep(Config_t* cfg_ptr, int k){

  Particles_t* Ptcl = cfg_ptr->ptcl_ptr;
  Perturb_t* Ptrb = cfg_ptr->ptrb_ptr;
  double* nout = Ptcl->nout;
  double* nfin = Ptcl->nfin;
  double* dt = Ptcl->dt;
  double* dptdp = cfg_ptr->ptcl_ptr->dptdp;
  double* pol = cfg_ptr->ptcl_ptr->pol;
  const double pw = get_pw(cfg_ptr->eqlb_ptr);
  double* time = get_time(Ptcl);
  double* tim1 = get_tim1(Ptcl);
  double trun = get_trun(cfg_ptr);
  int npert = get_npert(Ptrb);
  double* b = Ptcl->b;
  double* g = Ptcl->g;
  double* gp = Ptcl->gp;
  double* q = Ptcl->q;
  double* qp = Ptcl->qp;
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

  int n1,j,i,ndum;
  double xdum,ydum,rbb,dedb,deni,fac1,fac2,pdum,xdot,ydot;
  /* local temps, note compiler is suspicuous of those I have init 0..check */
  double y_[4];
  double d_[4];
  double e_[4];
  double a_[4]= {0.};
  double bx_[4] = {0.};
  double h_;
  double nt_;
  double c1_[4] = {0.};
  double zdots_;

  n1=4;

  dt[k]=dt0;

  nout[k] = .6 * (1. + copysign(1., pol[k]- pw) );
  nfin[k] =  .6 * (1. + copysign(1., time[k] - trun));
  nt_ = .6 * (1.  + copysign(1., pol[k]-.05 *pw));
  dt[k] = nt_ * dt[k] + (1 - nt_) * dt0;

  xdum = sqrt(pol[k])*cos(thet[k]);
  ydum = sqrt(pol[k])*sin(thet[k]);


  y_[0] = nt_ * pol[k] + (1-nt_)*xdum;
  y_[1] = nt_ * thet[k] + (1-nt_)*ydum;
  y_[2] = zet[k];
  y_[3] = rho[k];

  d_[0] = y_[0];
  d_[1] = y_[1];
  d_[2] = y_[2];
  d_[3] = y_[3];
  h_ = dt[k]/6.0 ;

  for(j=0; j < 4; j++){
    kfield(cfg_ptr, 1);

    if(npert ==  0){
      //goto 61
      if(k == 0){
        printf("FAILURE: integ without perturbations is not yet implimented\n");
      }
      return;
    }

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
    e_[0] = e_[0] * (1-nout[k])*(1-nfin[k]);
    e_[1] = e_[1] * (1-nout[k])*(1-nfin[k]);
    e_[2] = e_[2] * (1-nout[k])*(1-nfin[k]);
    e_[3] = e_[3] * (1-nout[k])*(1-nfin[k]);

    //goto 62, like 42, the answere to the universe and everything in it,  but twenty years older.
    for(i=0; i<n1; i++){
      if(j==0){
        //for(i=0; i<n1; i++){
        a_[i] = h_ * e_[i];
        y_[i] = d_[i] + 3 * a_[i];
        //}
      }
      if(j==1){
        // for(i=0; i<n1; i++){
        bx_[i] = h_ * e_[i];
        y_[i] = d_[i] + 3 * bx_[i];
        //}
      }
      if(j==2){
        //for(i=0; i<n1; i++){
        c1_[i] = h_ * e_[i];
        y_[i] = d_[i] + 6 * c1_[i];
        //}
      }
      if(j==3){
        //for(i=0; i<n1; i++){
        y_[i] = d_[i] + a_[i] + 2 * bx_[i] + 2 * c1_[i] + h_ * e_[i];
        //}
      }
    }
    // 40
    ndum = .6 * (1 - copysign(1. , y_[0]));
    pol[k] = nt_*y_[0] + (1-nt_)*(y_[0]*y_[0] + y_[1]*y_[1]);
    thet[k] = nt_*y_[1] + (1-nt_)*(atan(y_[1]/y_[0]) + ndum * M_PI);
    zet[k] = y_[2];
    rho[k] = y_[3];
    nout[k] = .6 * (1. + copysign(1., pol[k]-pw));

  } //j
  return;
}
