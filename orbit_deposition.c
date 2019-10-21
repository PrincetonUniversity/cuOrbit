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

const size_t MAXFNAME=255;
const int ibin=40;  /* used in temp array in pdedp_finalize */

struct Deposition {
  /* pdedp */
  char* pdedp_file;
  char* bfield_file;
  int nruns;
  bool compute_pdedp;
  bool initial_update_pdedp;
  double deposit_on_bins_after_fraction;
  double pdedp_dtsamp;
  double pdedp_dtav;
  int pdedp_tskip;
  double pdedp_otpup;
  bool pdedp_focusdep;
  bool pdedp_optimize;

  bool pdedp_initialized;
  /* the 5d dist */
  double* pde_pdedp;
  double* res_id_arr;
  /* internal vars */
  int pde_nbinDE;
  int pde_nbinDPz;
  int pde_nbinE;
  int pde_nbinPz;
  int pde_nbinmu;
  double* pde_varDE;
  double* pde_varDPz;
  double* pde_varE;
  double* pde_varPz;
  double* pde_varmu;
  double pde_Emin;
  double pde_Emax;
  double pde_Pzmin;
  double pde_Pzmax;
  double pde_mumin;
  double pde_mumax;
  double pde_DEmin;
  double pde_DEmax;
  double pde_DPzmin;
  double pde_DPzmax;
  /* xxx do we really need both? */
  double pde_maxDE;
  double pde_maxDPz;
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
  return (Deposition_t*)calloc(1, sizeof(Deposition_t));
}

void initialize_Deposition(Deposition_t* Depo_ptr, Config_t* cfg_ptr){

  Depo_ptr->pdedp_file = strndup(cfg_ptr->pdedp_file, MAXFNAME);
  Depo_ptr->bfield_file = strndup(cfg_ptr->bfield_file, MAXFNAME);
  Depo_ptr->nruns = cfg_ptr->nruns;
  Depo_ptr->compute_pdedp = cfg_ptr->compute_pdedp;
  Depo_ptr->initial_update_pdedp = cfg_ptr->initial_update_pdedp;
  Depo_ptr->deposit_on_bins_after_fraction = cfg_ptr->deposit_on_bins_after_fraction;
  Depo_ptr->pdedp_dtsamp = cfg_ptr->pdedp_dtsamp;
  Depo_ptr->pdedp_dtav = cfg_ptr->pdedp_dtav;
  Depo_ptr->pdedp_tskip = cfg_ptr->pdedp_tskip;
  Depo_ptr->pdedp_otpup = cfg_ptr->pdedp_otpup;
  Depo_ptr->pdedp_focusdep = cfg_ptr->pdedp_focusdep;
  Depo_ptr->pdedp_optimize = cfg_ptr->pdedp_optimize;
  if(Depo_ptr->initial_update_pdedp){
    printf("Warning, found initial pdedp_update TRUE, setting optimize FALSE\n");
    Depo_ptr->pdedp_optimize = 0;
  }
  Depo_ptr->res_id_arr_j = 10000;
  Depo_ptr->res_id_arr_i = 4;
  Depo_ptr->res_id_arr = (double*)calloc((unsigned)(cfg_ptr->nprt *
                                         Depo_ptr->res_id_arr_j *
                                          Depo_ptr->res_id_arr_i),
                                         sizeof(double));  /* hardcoded? ask mp*/

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
  Depo_ptr->pde_varDE = (double*)calloc((unsigned)Depo_ptr->pde_nbinDE, sz);
  Depo_ptr->pde_varDPz = (double*)calloc((unsigned)Depo_ptr->pde_nbinDPz, sz);
  Depo_ptr->pde_varE = (double*)calloc((unsigned)Depo_ptr->pde_nbinE, sz);
  Depo_ptr->pde_varPz = (double*)calloc((unsigned)Depo_ptr->pde_nbinPz, sz);
  Depo_ptr->pde_varmu = (double*)calloc((unsigned)Depo_ptr->pde_nbinmu, sz);
  Depo_ptr->pde_pdedp = (double*)calloc(sizeof_pdedp(Depo_ptr), sz);
  Depo_ptr->pdedp_initialized = true;
}

bool compute_pdedp(Deposition_t* Depo_ptr){
  return Depo_ptr->compute_pdedp;
}

bool initial_update_pdedp(Deposition_t* Depo_ptr){
  return Depo_ptr->initial_update_pdedp;
}

bool pdedp_optimize(Deposition_t* Depo_ptr){
  return Depo_ptr->pdedp_optimize;
}


double get_pdedp_dtsamp(Deposition_t* Depo_ptr){
  return Depo_ptr->pdedp_dtsamp;
}

int get_pdedp_tskip(Deposition_t* Depo_ptr){
  return Depo_ptr->pdedp_tskip;
}
void set_pdedp_tskip(Deposition_t* Depo_ptr, double pdedp_tskip){
  Depo_ptr->pdedp_tskip = pdedp_tskip;
  return;
}

bool get_initial_update_pdedp(Deposition_t* Depo_ptr){
  return Depo_ptr->initial_update_pdedp;
}
void set_initial_update_pdedp(Deposition_t* Depo_ptr, bool val){
  Depo_ptr->initial_update_pdedp = val;
  return;
}

bool get_pdedp_focusdep(Deposition_t* Depo_ptr){
  return Depo_ptr->pdedp_focusdep;
}
void set_pdedp_focusdep(Deposition_t* Depo_ptr, bool val){
  Depo_ptr->pdedp_focusdep = val;
  return;
}

size_t sizeof_pdedp(Deposition_t* Depo_ptr){
  return (size_t) (Depo_ptr->pde_nbinE * Depo_ptr->pde_nbinPz * Depo_ptr->pde_nbinmu *
          Depo_ptr->pde_nbinDE * Depo_ptr->pde_nbinDPz);
}

static inline int get_pdedp_ind(Deposition_t* Depo_ptr, int iE, int iPz, int imu, int iDE, int iDPz){
  const int ind = Depo_ptr->pde_nbinPz * Depo_ptr->pde_nbinmu * Depo_ptr->pde_nbinDE * Depo_ptr->pde_nbinDPz * iE +
      Depo_ptr->pde_nbinmu * Depo_ptr->pde_nbinDE * Depo_ptr->pde_nbinDPz * iPz +
      Depo_ptr->pde_nbinDE * Depo_ptr->pde_nbinDPz * imu +
      Depo_ptr->pde_nbinDPz * iDE + iDPz;
  return ind;
}

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


void pdedp_read(Deposition_t* Depo_ptr, Config_t* cfg_ptr){
  /* this reads the probability distribution data
    for P(DE,DP| E,P,mu) from a text file (UFILE)*/
  int i;
  int je, jp, jmu, jde, jdp;
  size_t sz;
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
  Depo_ptr->pdedp_dtav = Depo_ptr->pdedp_dtsamp/50.;

  /* continue skipping header */
  fscanf(ifp, "%*[^\n]\n");  /* skip line */
  fscanf(ifp, "%*[^\n]\n");  /* skip line */

  fscanf(ifp, "%d %*[^\n]\n", &Depo_ptr->pde_nbinDE);
  fscanf(ifp, "%d %*[^\n]\n", &Depo_ptr->pde_nbinDPz);
  fscanf(ifp, "%d %*[^\n]\n", &Depo_ptr->pde_nbinE);
  fscanf(ifp, "%d %*[^\n]\n", &Depo_ptr->pde_nbinPz);
  fscanf(ifp, "%d %*[^\n]\n", &Depo_ptr->pde_nbinmu);

  printf("\nNumber of points for p(DE,DP|E,P,mu) function:\n");
  printf("\t E\t P\t  mu\t DE\t DP\n");
  printf("\t%d\t%d\t%d\t%d\t%d\n",
        Depo_ptr->pde_nbinE, Depo_ptr->pde_nbinPz, Depo_ptr->pde_nbinmu,
        Depo_ptr->pde_nbinDE, Depo_ptr->pde_nbinDPz);

  /* malloc */
  if(! Depo_ptr->pdedp_initialized){
    initialize_pdedp(Depo_ptr);
  }

  for(i=0; i<Depo_ptr->pde_nbinDE; i++){
    fscanf(ifp, "%lf ", &Depo_ptr->pde_varDE[i]);
  }
  for(i=0; i<Depo_ptr->pde_nbinDPz; i++){
    fscanf(ifp, "%lf ", &Depo_ptr->pde_varDPz[i]);
  }
  for(i=0; i<Depo_ptr->pde_nbinE; i++){
    fscanf(ifp, "%lf ", &Depo_ptr->pde_varE[i]);
  }
  for(i=0; i<Depo_ptr->pde_nbinPz; i++){
    fscanf(ifp, "%lf ", &Depo_ptr->pde_varPz[i]);
  }
  for(i=0; i<Depo_ptr->pde_nbinmu; i++){
    fscanf(ifp, "%lf ", &Depo_ptr->pde_varmu[i]);
  }

  //     Define boundaries of computing grid.
  Depo_ptr->pde_Emin = Depo_ptr->pde_varE[0] - .5 * (
      Depo_ptr->pde_varE[Depo_ptr->pde_nbinE-1] - Depo_ptr->pde_varE[0]) /(
          Depo_ptr->pde_nbinE - 1.);

  Depo_ptr->pde_Emax = Depo_ptr->pde_varE[Depo_ptr->pde_nbinE - 1] + .5*(
      Depo_ptr->pde_varE[Depo_ptr->pde_nbinE] - Depo_ptr->pde_varE[0])/(
          Depo_ptr->pde_nbinE - 1.);

  Depo_ptr->pde_Pzmin = Depo_ptr->pde_varPz[0] - .5*(
      Depo_ptr->pde_varPz[Depo_ptr->pde_nbinPz-1] - Depo_ptr->pde_varPz[0])/(
          Depo_ptr->pde_nbinPz - 1.);

  Depo_ptr->pde_Pzmax = Depo_ptr->pde_varPz[Depo_ptr->pde_nbinPz-1] + .5*(
      Depo_ptr->pde_varPz[Depo_ptr->pde_nbinPz] - Depo_ptr->pde_varPz[0])/(
          Depo_ptr->pde_nbinPz - 1.);

  Depo_ptr->pde_mumin = Depo_ptr->pde_varmu[0] - .5*(
      Depo_ptr->pde_varmu[Depo_ptr->pde_nbinmu - 1] - Depo_ptr->pde_varmu[0])/(
          Depo_ptr->pde_nbinmu - 1.);

  Depo_ptr->pde_mumax = Depo_ptr->pde_varmu[Depo_ptr->pde_nbinmu-1] + .5*(
      Depo_ptr->pde_varmu[Depo_ptr->pde_nbinmu - 1] - Depo_ptr->pde_varmu[0])/(
          Depo_ptr->pde_nbinmu - 1.);

  Depo_ptr->pde_DEmin = Depo_ptr->pde_varDE[0] - .5*(
      Depo_ptr->pde_varDE[Depo_ptr->pde_nbinDE-1] - Depo_ptr->pde_varDE[0])/(
          Depo_ptr->pde_nbinDE - 1.);

  Depo_ptr->pde_DEmax = Depo_ptr->pde_varDE[Depo_ptr->pde_nbinDE - 1] + .5*(
      Depo_ptr->pde_varDE[Depo_ptr->pde_nbinDE - 1] - Depo_ptr->pde_varDE[0])/(
          Depo_ptr->pde_nbinDE - 1.);

  Depo_ptr->pde_DPzmin = Depo_ptr->pde_varDPz[0] - .5*(
      Depo_ptr->pde_varDPz[Depo_ptr->pde_nbinDPz - 1] - Depo_ptr->pde_varDPz[0]) / (
          Depo_ptr->pde_nbinDPz - 1. );

  Depo_ptr->pde_DPzmax = Depo_ptr->pde_varDPz[Depo_ptr->pde_nbinDPz - 1] + .5*(
      Depo_ptr->pde_varDPz[Depo_ptr->pde_nbinDPz - 1] -
      Depo_ptr->pde_varDPz[0]) / (Depo_ptr->pde_nbinDPz - 1.);

  /* read the 5d distribution in */
  /*  careful, the loop order is trecherous */
  for(je=0; je < Depo_ptr->pde_nbinE; je++){
    for(jp=0; jp < Depo_ptr->pde_nbinPz; jp++){
      for(jmu=0; jmu < Depo_ptr->pde_nbinmu; jmu++){
        for(jde=0; jde < Depo_ptr->pde_nbinDE; jde++){
          for(jdp=0; jdp < Depo_ptr->pde_nbinDPz; jdp++){
            i = Depo_ptr->pde_nbinPz * Depo_ptr->pde_nbinE * Depo_ptr->pde_nbinDPz * Depo_ptr->pde_nbinDE * jmu +
                Depo_ptr->pde_nbinE * Depo_ptr->pde_nbinDPz * Depo_ptr->pde_nbinDE * jp +
                Depo_ptr->pde_nbinDPz * Depo_ptr->pde_nbinDE * je +
                Depo_ptr->pde_nbinDE * jdp + jde;
            fscanf(ifp, "%lf ", &Depo_ptr->pde_pdedp[i]);
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

  Depo_ptr->pde_Emax = 110.;
  Depo_ptr->pde_Pzmin = -1.2;
  Depo_ptr->pde_Pzmax =  0.7;
  Depo_ptr->pde_mumin = 0.;  /* muB0/E */
  Depo_ptr->pde_mumax = 1.4;
  Depo_ptr->pde_DEmin = -0.100000;
  Depo_ptr->pde_DEmax = 0.100000;
  Depo_ptr->pde_DPzmin = -0.00100000;
  Depo_ptr->pde_DPzmax = 0.00100000;
  Depo_ptr->pde_nbinE = 12;  /* number of bins */
  Depo_ptr->pde_nbinPz=40;
  Depo_ptr->pde_nbinmu = 16;
  Depo_ptr->pde_nbinDE=29;   /* this must be an ODD number */
  assert(Depo_ptr->pde_nbinDE % 2 == 1);
  Depo_ptr->pde_nbinDPz = 29;  /* this must be an ODD number; */
  assert(Depo_ptr->pde_nbinDE % 2 == 1);

  /* malloc */
  if(! Depo_ptr->pdedp_initialized){
    initialize_pdedp(Depo_ptr);
  }

  /* fill in grid vectors */
  /* Energy  */
  stp = (Depo_ptr->pde_Emax - Depo_ptr->pde_Emin) / Depo_ptr->pde_nbinE;
  for(k=0; k < Depo_ptr->pde_nbinE; k++){
    Depo_ptr->pde_varE[k] = Depo_ptr->pde_Emin + k * stp + stp/2.;
  }

  /* Pz */
  stp = (Depo_ptr->pde_Pzmax - Depo_ptr->pde_Pzmin) / Depo_ptr->pde_nbinPz;
  for(k=0; k < Depo_ptr->pde_nbinPz; k++){
    Depo_ptr->pde_varPz[k] = Depo_ptr->pde_Pzmin + k * stp + stp/2.;
  }

  /* mu Bo/E */
  stp = (Depo_ptr->pde_mumax - Depo_ptr->pde_mumin) / Depo_ptr->pde_nbinmu;
  for(int k=0; k < Depo_ptr->pde_nbinmu; k++){
    Depo_ptr->pde_varmu[k] = Depo_ptr->pde_mumin + k * stp + stp/2.;
  }

  /* Delta-Energy */
  stp = (Depo_ptr->pde_DEmax - Depo_ptr->pde_DEmin) / Depo_ptr->pde_nbinDE;
  for(k=0; k < Depo_ptr->pde_nbinDE; k++){
    Depo_ptr->pde_varDE[k] = Depo_ptr->pde_DEmin + k * stp + stp/2.;
  }

  /*Delta-Pz */
  stp = (Depo_ptr->pde_DPzmax - Depo_ptr->pde_DPzmin) / Depo_ptr->pde_nbinDPz;
  for(k=0; k < Depo_ptr->pde_nbinDPz; k++){
    Depo_ptr->pde_varDPz[k] = Depo_ptr->pde_DPzmin + k * stp + stp/2.;
  }

  Depo_ptr->pde_pdedp = calloc(sizeof_pdedp(Depo_ptr), sizeof(double));


  /*      -------------------------------------------------------
          Initialize variables to monitor the maximum
          kicks in energy and Pz. These are used to optimize
          the (DE,DPz) range on-the-fly if the flag
          pde_optimize is set true*/
  Depo_ptr->pde_maxDE = 0.;
  Depo_ptr->pde_maxDPz = 0.;

  return;
}

void class_kdomain(Config_t* cfg_ptr, int k){


  /*    -  Orbit classification
      otp=1, co-passing confined, otp = 2, co-passing lost
      otp=3, ctr-passing confined, otp = 4, ctr-passing lost
      otp=5, trapped confined, otp = 6, trapped lost
      otp=7, stagnation, otp = 8, conf potato, otp = 9, trap potato
  */
  int ndum,ntot,ntlos,km,ku,kv,kt,j;
  double xdum,zdum,edum,pdum,podum,dum,dum1,dum2,vol;
  double pzdum,elt,ert,eax,etp1,etp2,mu,mube;
  double evmin,evmax,ev0,mub,rhod,elt2;
  double Epot,Ekin, E_ax,E_min_lcfs,E_max_lcfs;
  double E_th0,E_thpi,mu2;
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
  mub = mu*bax;
  pdum = (g[k]*rho[k] - pol[k])/pw;     /* pz0[k] */
  edum = en[k]*ekev/engn;
  pzdum = pdum*pw;      /* pz0[k]*pw   Pz */
  rhod = (pzdum + pol[k])/g[k] ; /*  particle rho */
  mube = rmu[k]*bax/en[k];   /*  mu*Bax/E */
  dum = pow(pzdum+pw,2) * ekev / (
      engn * 2 * pow(gfun(cfg_ptr->eqlb_ptr, pw),2));

  /* get potential energies */
  Epot = pol2pot(cfg_ptr, pol[k]);
  E_ax = pol2pot(cfg_ptr, 0.);
  E_min_lcfs = pol2pot(cfg_ptr, pw);
  E_max_lcfs = pol2pot(cfg_ptr, -pzdum);
  E_th0 = pol2pot(cfg_ptr, -pzdum);
  E_thpi = pol2pot(cfg_ptr, -pzdum);

  ert =  pow(dum*bmin, 2) + mu*bmin + E_min_lcfs;
  elt =  pow(dum*bmax, 2) + mu*bmax + E_max_lcfs;
  dum = pow(rhod, 2) * ekev/(2*engn);
  elt2 =  pow(dum*bmax,2) + mu*bmax + E_max_lcfs;
  eax = pow(ekev*pzdum,2)* pow(bax,2)/(
      engn * 2 * pow(gfun(cfg_ptr->eqlb_ptr, 0.),2)) + mu*bax + E_ax;
  etp1 = mu*bfield(cfg_ptr->eqlb_ptr, -pzdum,0.) + E_th0;
  etp2 = mu*bfield(cfg_ptr->eqlb_ptr, -pzdum, M_PI) + E_thpi;
  evmin = mu*bmin + E_min_lcfs;
  evmax = mu*bmax + E_max_lcfs;
  ev0 = mu*bax + E_ax;

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

  if pde_pdedp lives on card, this is trivial there*/
  double* const pde_pdedp = Depo_ptr->pde_pdedp;
  double sum_p[40][40][40];

  printf("-> Finalize pDEDP computation ...\n");
  /* Get indexes of (DE,DPz)=(0,0) bin */
  int iDE0 = get_bin(0., Depo_ptr->pde_varDE, Depo_ptr->pde_nbinDE);
  int iDPz0 = get_bin(0., Depo_ptr->pde_varDPz, Depo_ptr->pde_nbinDPz);

  /* make sure center bin has exactly DE=0,DPz=0 */
  /* set to zero */
  Depo_ptr->pde_varDE[iDE0]=0.;
  Depo_ptr->pde_varDPz[iDPz0]=0.;

  /*       Get average number of counts/bin from non-empty bins
           and fill in empty bins */
  double  cnt_aver=0.;
  int cnt_;
  int ind;
  int nbins=0;
  for(int iE=0; iE < Depo_ptr->pde_nbinE; iE++){
    for(int iPz=0; iPz < Depo_ptr->pde_nbinPz; iPz++){
      for(int imu=0; imu < Depo_ptr->pde_nbinmu; imu++){
        cnt_=0.;
        sum_p[iE][iPz][imu]=0.;
        for(int iDE=0; iDE < Depo_ptr->pde_nbinDE; iDE++){
          for(int iDPz=0; iDPz < Depo_ptr->pde_nbinDPz; iDPz++){
            cnt_ += pde_pdedp[
                get_pdedp_ind(Depo_ptr, iE, iPz, imu, iDE, iDPz)];
          }
        }
        if(cnt_ > 0) {
          /* update*/
          cnt_aver += cnt_;
          sum_p[iE][iPz][imu]=cnt_;
          nbins += 1;
        } else {
          /* fill with 1 count at (DE,DPz)=(0,0)*/
          pde_pdedp[get_pdedp_ind(Depo_ptr, iE, iPz, imu, iDE0, iDPz0)] = 1.;
          sum_p[iE][iPz][imu]=1.;
        }
      }
    }
  }

  if(nbins > 0) cnt_aver /= nbins;

  /* Normalize */
  for(int iE=0; iE < Depo_ptr->pde_nbinE; iE++){
    for(int iPz=0; iPz < Depo_ptr->pde_nbinPz; iPz++){
      for(int imu=0; imu < Depo_ptr->pde_nbinmu; imu++){
        for(int iDE=0; iDE < Depo_ptr->pde_nbinDE; iDE++){
          for(int iDPz=0; iDPz < Depo_ptr->pde_nbinDPz; iDPz++){
            pde_pdedp[get_pdedp_ind(Depo_ptr, iE, iPz, imu, iDE, iDPz)] *=
                cnt_aver/sum_p[iE][iPz][imu];
          }
        }
      }
    }
  }
  printf("\n-> p(DE,DPz|E,Pz,mu) matrices normalized\n");
  printf("   - Average number of counts: %f\n\n", cnt_aver);

  return;
}



void pdedp_out(Deposition_t* Depo_ptr){
  /* writes out the dist file */
  int k;
  int ind;
  FILE *ofp;
  ofp = fopen(Depo_ptr->pdedp_file, "w");  /* add f status */
  if (ofp == NULL) {
    fprintf(stderr, "\nCan't open output file %s!\n", Depo_ptr->pdedp_file);
    exit(1);
  }
  printf("\nOutputting P(DE,DP) file %s\n",  Depo_ptr->pdedp_file);

  /* place holders, apparently this will come from transp sometime
     right now, just mimic old behavior. */
  const int lshot = 123456;
  const char *date = "Apr. 2018";
  const int nd = 5;     /* 5-D data */
  const int nq = 0;     /* unknown parameter */
  const int nr = 6;     /* #of decimal places f13.6 */
  const int np = 0;     /* process code */
  const int ns = 1;    /*  # of scalars */
  const char *dev = "D3D";
  const char *labelx = "DEstep                "; /* *20 */
  const char *unitsx = "   kev    ";             /* *10 */
  const char *labely = "DPsteps               "; /* *20 */
  const char *unitsy = "          ";             /* *10 */
  const char *labelu = "Evar                  "; /* *20 */
  const char *unitsu = "   keV    ";             /* *10 */
  const char *labelv = "Pvar                  "; /* *20 */
  const char *unitsv = "          ";             /* *10 */
  const char *labelw = "Muvar                ";  /* *20 */
  const char *unitsw = "          ";             /* *10 */
  const char *com = ";----END-OF-DATA-----------------COMMENTS:-----------";
  const char *com2 = "UFILE WRITTEN BY ORBIT, see PDEDP_OUT";
  const char *com3 = "SMOOTHING FACTORS, DELAY FACTORS:";
  const char *com4 = "       NONE";
  const char *com5 = "USER COMMENTS:";
  const char *com6 = "       TEST FILE";

  /* xxx, note I removed a lot of whitespace formatting,
     can add back in later*/

  fprintf(ofp, "%d %s %d %d %d ;-SHOT #- F(X) DATA -PDEDP_OUT- \n",
          lshot, dev, nd, nq, nr);

  fprintf(ofp, "%s ;-SHOT DATE-  UFILES ASCII FILE SYSTEM\n", date);

  fprintf(ofp, "%d ;-NUMBER OF ASSOCIATED SCALAR QUANTITIES-\n", ns);

  fprintf(ofp,"1.0000E+00                   ;-SCALAR, LABEL FOLLOWS:\n");

  fprintf(ofp, "%s %s ;-INDEPENDENT VARIABLE LABEL: X-\n", labelx, unitsx);

  fprintf(ofp, "%s %s ;-INDEPENDENT VARIABLE LABEL: Y-\n", labely, unitsy);

  fprintf(ofp, "%s %s ;-INDEPENDENT VARIABLE LABEL: U-\n",
          labelu, unitsu);

  fprintf(ofp, "%s %s ;-INDEPENDENT VARIABLE LABEL: V-\n", labelv, unitsv);

  fprintf(ofp, "%s %s ;-INDEPENDENT VARIABLE LABEL: W-\n", labelw, unitsw);

  fprintf(ofp, "%f ; TSTEPSIM  - TIME STEP USED IN SIMULATION [ms]\n",
          Depo_ptr->pdedp_dtsamp);

  fprintf(ofp, " PROBABILITY DATA              ;-DEPENDENT VARIABLE LABEL-\n");

  fprintf(ofp, "%d ;-PROC CODE- 0:RAW 1:AVG 2:SM 3:AVG+SM\n", np);

  fprintf(ofp, "%d ;-# OF X PTS-\n", Depo_ptr->pde_nbinDE);

  fprintf(ofp, "%d ;-# OF Y PTS-\n", Depo_ptr->pde_nbinDPz);

  fprintf(ofp, "%d ;-# OF U PTS-\n", Depo_ptr->pde_nbinE);

  fprintf(ofp, "%d ;-# OF V PTS-\n", Depo_ptr->pde_nbinPz);

  fprintf(ofp, "%d ;-# OF W PTS- X,Y,U,V,W,F(X,Y,U,V,W) DATA FOLLOW:\n",
          Depo_ptr->pde_nbinmu);

  /*       -------------------------------------------------------
         Make sure center bin has exactly DE=0,DPz=0

         Get indexes of (DE,DPz)=(0,0) bin */
  double valdum=0.;
  int iDE0 =  get_bin(valdum, Depo_ptr->pde_varDE, Depo_ptr->pde_nbinDE);
  int iDPz0 = get_bin(valdum, Depo_ptr->pde_varDPz, Depo_ptr->pde_nbinDPz);

   /* set to zero */
  Depo_ptr->pde_varDE[iDE0]=0.;
  Depo_ptr->pde_varDPz[iDPz0]=0.;

  /* Write grid vectors */

  /* variables DEbins,DPbins */
  for(k=0; k < Depo_ptr->pde_nbinDE; k++){
    fprintf(ofp, "%f ", Depo_ptr->pde_varDE[k]);
  }
  fprintf(ofp,"\n");

  for(k=0; k < Depo_ptr->pde_nbinDPz; k++){
    fprintf(ofp, "%f ", Depo_ptr->pde_varDPz[k]);
  }
  fprintf(ofp,"\n");

  /* variables Ebins,Pbins,Mubins */
  for(k=0; k < Depo_ptr->pde_nbinE; k++){
    fprintf(ofp, "%f ", Depo_ptr->pde_varE[k]);
  }
  fprintf(ofp,"\n");

  for(k=0; k < Depo_ptr->pde_nbinPz; k++){
    fprintf(ofp, "%f ", Depo_ptr->pde_varPz[k]);
  }
  fprintf(ofp,"\n");

  for(k=0; k < Depo_ptr->pde_nbinmu; k++){
    fprintf(ofp, "%f ", Depo_ptr->pde_varmu[k]);
  }
  fprintf(ofp,"\n");


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

  for(int iE=0; iE < Depo_ptr->pde_nbinE; iE++){
    for(int iPz=0; iPz < Depo_ptr->pde_nbinPz; iPz++){
      for(int imu=0; imu < Depo_ptr->pde_nbinmu; imu++){

        /* compute normalization factor for this bin
           !          pnorm=0.0
           !          do iDE=1,pde_nbinDE
           !            do iDPz=1,pde_nbinDPz
           !              pnorm=pnorm+pde_Pdedp(iDE,iDPz,iE,iPz,imu)
           !            enddo
           !          enddo
           !
           !          ! write normalized p(DE,DPz|E,Pz,mu) */

        for(int iDE=0; iDE < Depo_ptr->pde_nbinDE; iDE++){
          /*            ! normalize total probability to 1
              !            if(pnorm.gt.0.) then
              !              do iDPz=1,pde_nbinDPz
              !                pde_Pdedp(iDE,iDPz,iE,iPz,imu)=
              !     >             pde_Pdedp(iDE,iDPz,iE,iPz,imu)/pnorm
              !               enddo
              !            endif */
          for(int iDPz=0; iDPz < Depo_ptr->pde_nbinDPz; iDPz++){
            fprintf(ofp, "%f ", Depo_ptr->pde_pdedp[
                get_pdedp_ind(Depo_ptr, iE, iPz, imu, iDE, iDPz)]);
          }
          fprintf(ofp, "\n");
        }
      }
    }
  }

  fprintf(ofp, "\n%s\n", com);
  fprintf(ofp, "%s\n", com2);
  fprintf(ofp, "%s\n", com3);
  fprintf(ofp, "%s\n", com4);
  fprintf(ofp, "%s\n", com5);
  fprintf(ofp, "%s\n", com6);

  fclose(ofp);

  return;
}

static inline int get_res_id_ind(Config_t* cfg_ptr, int k, int j, int i){
  /*  fname indices, sorry */
  const int nj = cfg_ptr->depo_ptr->res_id_arr_j;
  const int ni = cfg_ptr->depo_ptr->res_id_arr_i;

  return nj*ni*k + ni*j + i;
}

void pdedp_rcrd_resid(Config_t* cfg_ptr, Deposition_t* Depo_ptr){
  /*      pde_tskip : reduce loops, skip time steps */
  int j, j2, j3, k, ind;
  double dtdum=Depo_ptr->pdedp_tskip *
      1.0E3 * cfg_ptr->dt0 / get_omeg0(cfg_ptr->ptrb_ptr); /* [ms] */
  int jstart = (int)(0.2 * Depo_ptr->pdedp_dtsamp / dtdum -1 );  /* zero ind */
  /*  in original code, this was nstep,
      but called at end of loop when it should have achieved nstep_all value,
      I think, which is not "private" */
  int nsamples= cfg_ptr->nstep_all / Depo_ptr->pdedp_tskip;
  double newE, newPz;
  double Eav, Pzav, Muav;
  double dedum, dpzdum;
  int iE, iDE, iPz, iDPz, imu;
  double * const res_id_arr = Depo_ptr->res_id_arr;


  printf("   ... computing p(DE,DPz) matrix ...\n");
  int nintv = (int)( Depo_ptr->pdedp_dtsamp / dtdum);
  int Nav = imax(1, (int)( Depo_ptr->pdedp_dtav / dtdum));

  for(k=0; k < cfg_ptr->nprt; k++){
    for(j=jstart; j < nsamples; j++){
      Eav = 0.;
      Pzav = 0.;
      Muav = 0.;
      for(j3=0; j3 <= Nav-1; j3++){
        /* j is already zero indices */
        Eav += res_id_arr[get_res_id_ind(cfg_ptr,k,j-j3,0)];  /* E */
        Pzav += res_id_arr[get_res_id_ind(cfg_ptr,k,j-j3,1)];  /* Pz */
        Muav += res_id_arr[get_res_id_ind(cfg_ptr,k,j-j3,2)];  /* mu */
      }

      Eav /= ((double)Nav);
      Pzav /= ((double)Nav);
      Muav /= ((double)Nav);

      iE =  get_bin(Eav, Depo_ptr->pde_varDE, Depo_ptr->pde_nbinE);
      iPz =  get_bin(Pzav, Depo_ptr->pde_varDPz, Depo_ptr->pde_nbinPz);
      imu =  get_bin(Muav, Depo_ptr->pde_varmu, Depo_ptr->pde_nbinmu);

      j2 = nintv;
      /* zero indices */
      ind = get_res_id_ind(cfg_ptr,k,j +j2 -j3,3);
      if (j + j2 +1 < nsamples &&
          res_id_arr[ind] > 0 &&
          res_id_arr[ind] < 1)
      {
        newE = 0.;
        newPz = 0.;
        for(j3=0; j3 < Nav-1; j3++){
          newE  += res_id_arr[get_res_id_ind(cfg_ptr, k, j+j2-j3, 0)] / ((double)Nav);
          newPz += res_id_arr[get_res_id_ind(cfg_ptr, k, j+j2-j3, 1)] / ((double)Nav);
        }
        dedum = newE - Eav;
        dpzdum = newPz - Pzav;

        iDE =  get_bin(dedum, Depo_ptr->pde_varDE, Depo_ptr->pde_nbinDE);
        iDPz =  get_bin(dedum, Depo_ptr->pde_varDPz, Depo_ptr->pde_nbinDPz);

        if (Depo_ptr->pdedp_optimize == 1 &&
            newE > 0 &&
            iE >= 1 &&
            iE <= Depo_ptr->pde_nbinE  &&
            iPz >= 1 &&
            iPz <= Depo_ptr->pde_nbinPz &&
            imu >= 1 &&
            imu <= Depo_ptr->pde_nbinmu){
          Depo_ptr->pde_maxDE = fmax(fabs(1.05 * dedum),
                                    Depo_ptr->pde_maxDE);
          Depo_ptr->pde_maxDPz = fmax(fabs(1.05 * dpzdum),
                                     Depo_ptr->pde_maxDPz);
        }
      }
    }  /* j */
  }    /* k */
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
  Particles_t* ptcl_ptr = cfg_ptr->ptcl_ptr;
  double** const GD = get_GD(eqlb_ptr);
  double** const B = get_B(eqlb_ptr);
  const double pw = get_pw(eqlb_ptr);
  const double engn = get_engn(cfg_ptr);

  /*don't touch grid if rel. difference is smaller than this */
  const double dthres = 1.;
  int recompute;
  int k;
  double mumin, mumax, Pzmin, Pzmax, pde_engn;
  double fctE, fctPz, fctMu;
  double demax_old, dPzmax_old, stp;

  printf("Checking pDEDP boundaries ...\n");

  /* First, find boundary for mu and Pz based on
     energy range used in pDEDP calculation */

  /* maximum energy, normalized units */

  pde_engn = depo_ptr->pde_Emax * 10.533 *
      get_prot(ptcl_ptr) * pow(GD[0][0],2) /
      (get_rmaj(eqlb_ptr) * get_zprt(ptcl_ptr) * get_bkg(cfg_ptr) * pow(B[0][0],2));
  /* cf. engn definition in initial.f */

  /* min/max B field */
  const double Bmn = bfield(eqlb_ptr, pw, 0.);
  const double Bmx = bfield(eqlb_ptr, pw, M_PI);

  /* redefine range of mu
     add buffer to the actual range */
  mumin = 0.;
  mumax = 1./Bmn * (depo_ptr->pde_nbinmu + 1.) / depo_ptr->pde_nbinmu;

  /* upper limit for Pz, add buffer */
  Pzmax = gfun(eqlb_ptr, 0)/pw*sqrt(2. * engn) *
      (depo_ptr->pde_nbinPz + 1.) / depo_ptr->pde_nbinPz;

  /* lower limit for Pz, add buffer */
  Pzmin = -1. - gfun(eqlb_ptr, pw)/pw*sqrt(2. * engn) / Bmx *
      (depo_ptr->pde_nbinPz + 1.) /depo_ptr->pde_nbinPz;


  /* Check wheter the Pz,mu range needs to be adjusted. */

  fctMu = (depo_ptr->pde_mumax - mumax ) / depo_ptr->pde_mumax;
  fctPz = fmax(fabs((depo_ptr->pde_Pzmax - Pzmax) / depo_ptr->pde_Pzmax),
              fabs((depo_ptr->pde_Pzmax - Pzmin) / depo_ptr->pde_Pzmin));

    if(fctMu > dthres || fctPz > dthres ||
       Pzmin < depo_ptr->pde_Pzmin || Pzmax > depo_ptr->pde_Pzmax ||
       mumax > depo_ptr->pde_mumax) {
      /* update */

      /* display info with updated grid */
      printf("  -> New Pz,mu grid computed:\n");
      printf("  original: Pz1= %f,    Pz2= %f,    mu=%f\n",
             depo_ptr->pde_Pzmin,
             depo_ptr->pde_Pzmax,
             depo_ptr->pde_mumax);
      printf("  updated: Pz1= %f,    Pz2= %f,    mu=%f\n",
             Pzmin, Pzmax, mumax);

      depo_ptr->pde_Pzmax = Pzmax;
      depo_ptr->pde_Pzmin = Pzmin;
      depo_ptr->pde_mumin = 0.;
      depo_ptr->pde_mumax = mumax;


      /* Pz */
      stp=(depo_ptr->pde_Pzmax - depo_ptr->pde_Pzmin)/depo_ptr->pde_nbinPz;
      for(k=0; k < depo_ptr->pde_nbinPz; k++){
        depo_ptr->pde_varPz[k] = depo_ptr->pde_Pzmin + ((double)k)* stp + stp/2.;
      }

      /* mu Bo/E */
      stp=(depo_ptr->pde_mumax - depo_ptr->pde_mumin) / depo_ptr->pde_nbinmu;
      for(k=0; k < depo_ptr->pde_nbinmu; k++){
        depo_ptr->pde_varmu[k] = depo_ptr->pde_mumin +((double)k) * stp + stp/2.;
      }

      recompute = 1;
    } else {
      printf(" -> Pz,mu grid looks OK - \n\n");
    }

  /* Check wheter the DE,DPz range needs to be
     adjusted. */
  fctE=(depo_ptr->pde_maxDE - depo_ptr->pde_DEmax) / depo_ptr->pde_DEmax;
  fctPz=(depo_ptr->pde_maxDPz - depo_ptr->pde_DPzmax) / depo_ptr->pde_DPzmax;

  if(fabs(fctE) < dthres && fabs(fctPz) < dthres &&
     fctE <= 1  && fctPz <= 1){
    printf("  -> DE,DPz grid looks OK - \n\n");
  } else {
    /* update range */
    /* keep copy  */
    double demax_old = depo_ptr->pde_DEmax;
    double dPzmax_old = depo_ptr->pde_DPzmax;


    /* new values */
    depo_ptr->pde_DEmax = (1. + fctE) * depo_ptr->pde_DEmax;
    depo_ptr->pde_DPzmax = (1. + fctPz) * depo_ptr->pde_DPzmax;


    /* round off [doesn't need to preserve precision] */
    depo_ptr->pde_DEmax = 1E-2 * (ceil(1E2 * depo_ptr->pde_DEmax));
    depo_ptr->pde_DPzmax = 1E-4 * (ceil(1E4 * depo_ptr->pde_DPzmax));

    /* symmetric grid */
    depo_ptr->pde_DEmin = -depo_ptr->pde_DEmax;
    depo_ptr->pde_DPzmin = -depo_ptr->pde_DPzmax;


    /*  define new grid */
    /* Delta E */
    stp = 2. * depo_ptr->pde_DEmax / depo_ptr->pde_nbinDE;
    for(k=0; k < depo_ptr->pde_nbinDE; k++){
      depo_ptr->pde_varDE[k] = -depo_ptr->pde_DEmax + ((double)k) * stp + stp/2.;
    }
    /* Delta Pz */
    stp = 2. * depo_ptr->pde_DPzmax / depo_ptr->pde_nbinDPz;
    for(k=0; k<depo_ptr->pde_nbinDPz; k++){
      depo_ptr->pde_varDPz[k] = -depo_ptr->pde_DPzmax + ((double)k) * stp + stp/2.;
    }

    /* update flag */
    recompute = 1;
  }
  /* If range of Pz, mu, DE, or DPz has changed, */
  /* update pDEDP computation */
  if(recompute != 0){
    /* reset pDEDP before resampling */
    memset(depo_ptr->pde_pdedp, 0, sizeof_pdedp(depo_ptr)*sizeof(double));  /* zero */
    /* compute pDEDP probability based on  updated range */
    pdedp_rcrd_resid(cfg_ptr, depo_ptr);
  }

  /* reset flag - no further optimizations */

  depo_ptr->pdedp_optimize = false;

  printf("- done.\n\n");

  return;
}

void fulldepmp(Config_t* cfg_ptr, Deposition_t* depo_ptr){
  /* loops over k and k/2, can live on device one day */

  int k, kd, np2, ier;
  double dum,xproj,tdum;
  double edum,pzdum,mudum,bdum;
  double nprt;

  const double einj1 = depo_ptr->pde_Emin;  /* [keV] */
  const double einj2 = depo_ptr->pde_Emax;
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
  double* g = get_g(cfg_ptr->ptcl_ptr);
  int* otp = get_otp(Ptcl);

  /* Full deposition*/
  np2 = .5* cfg_ptr->nprt;

  /* outside-co moving */
  for(kd=0; kd < np2; kd++){
    ptch[kd] = rand_double();  /* rand */
    thet[kd] = 0.;
    dt[kd] = dt0;
    pol[kd] =  .001 + .999 * rand_double();
    pol[kd] = pol[kd]*pw;
    zet[kd] = rand_double() * 2. * M_PI;
    en[kd] = (einj1 + rand_double() * (einj2 - einj1))
        * engn / ekev;   /* kinetic energy */
  }

  /* -inside-counter moving */
  for(kd=0; kd < np2; kd++){
    ptch[kd] = -rand_double();
    thet[kd] = M_PI;
    dt[kd] = dt0;
    pol[kd] = .001 + .999*rand_double();
    pol[kd] = pol[kd] * pw;
    zet[kd] = rand_double() * 2.* M_PI;
    en[kd] = (einj1+rand_double()*(einj2-einj1))*engn/ekev; /* kinetic energy */
  }

  nprt = np2;  /* xxx not sure about this? */
  printf(" fulldepmp deposit\n");
  /*!      call field(nprt) */
  for(k=0; k<nprt; k++){
    kfield(cfg_ptr, k);
    rho[k] = ptch[k]*sqrt(2. * en[k]) /b[k];
    rmu[k] = en[k] / b[k] -
        .5 * rho[k] * rho[k] * b[k];
    en[k] = en[k] + pot[k];

    /* DEBUG: */
    edum = en[k];    /* !*ekev/engn*/
    pzdum = (g[k]*rho[k] - pol[k]) / pw;
    bdum = bfield(cfg_ptr->eqlb_ptr, pol[k], thet[k]);
    mudum = rmu[k];  /* /en[k] */
    /*!       mudum=edum/bdum*(1.0-ptch[k]*ptch[k]) */
  }

  /* re-sample lost particles */
  class_domain(cfg_ptr);
  fullredepmp(cfg_ptr, depo_ptr);

  return;
}

void fullredepmp(Config_t* cfg_ptr, Deposition_t* depo_ptr){
  int k,kd,np2,nlost,nmaxs,imaxs;
  double dum,xproj,tdum;

  const double einj1 = depo_ptr->pde_Emin;  /* [keV] */
  const double einj2 = depo_ptr->pde_Emax;
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
  double* g = get_g(Ptcl);
  int* otp = get_otp(Ptcl);

  nmaxs=5E3; /* max number of attemps to redeposit particle */
  /* -  Full deposition, */
  np2 = .5 * cfg_ptr->nprt;
  /* np2 = cfg_ptr->nprt; */

  /* dont print this recursively */
  if(! depo_ptr->recursing){
    printf("Entering FULLREDEPMP...\n");
  }

  nlost=0;

  /* -outside-co moving */
  for(kd=0; kd < np2; kd++){
    imaxs=1;

    if(depo_ptr->pdedp_focusdep && imaxs < nmaxs){
      check_res_ptc(cfg_ptr, kd);
      imaxs += 1;
    }

    if(pol[kd] >= pw || otp[kd] == 2 ||
       otp[kd] == 4 || otp[kd] == 6){
      /* lost, replace it */

      nlost+=1;
      ptch[kd] = rand_double();
      thet[kd] = 0.;
      dt[kd] = dt0;
      pol[kd] =  .002*pw + .997*pw*rand_double();
      zet[kd] = 2. * M_PI * rand_double();
      /* kinetic energy*/
      en[kd] = (einj1 + rand_double() * (einj2-einj1)) * engn / ekev;
      kfield(cfg_ptr, kd);
      rho[kd] = ptch[kd]*sqrt(2.*en[kd])/b[kd];
      rmu[kd] = en[kd]/b[kd] - .5*rho[kd]*rho[kd]*b[kd];
      en[kd] = en[kd] + pot[kd];
    }
  }

  /* -inside-counter moving */
  for(kd=0; kd < np2; kd++)
  {
    imaxs=1;

    if(cfg_ptr->pdedp_focusdep && imaxs < nmaxs){
      check_res_ptc(cfg_ptr, kd);
      imaxs += 1;
    }
    if(pol[kd] >= pw || otp[kd] == 2 ||
       otp[kd] == 4 || otp[kd] == 6){
      /* lost, replace it */
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
      en[kd] = en[kd] + pot[kd];
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
  double edum, pzdum, mudum, zdum;
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
  double* b = get_b(Ptcl);
  double* g = get_g(Ptcl);

  ptot = 0.;
  pmax = 0.;

  edum = en[kd]*ekev/engn;
  pzdum=(g[kd]*rho[kd] - pol[kd])/pw;
  mudum=rmu[kd]/en[kd];

  const int iE = get_bin(edum, depo_ptr->pde_varE, depo_ptr->pde_nbinE);
  const int iPz = get_bin(pzdum, depo_ptr->pde_varPz, depo_ptr->pde_nbinPz);
  const int iMu = get_bin(mudum, depo_ptr->pde_varmu, depo_ptr->pde_nbinmu);

  if(iE <= 0 || iE >= depo_ptr->pde_nbinE ||
       iPz <= 0 || iPz >= depo_ptr->pde_nbinPz ||
     iMu <= 0 || iMu >= depo_ptr->pde_nbinmu){
    pol[kd] =2. * pw;
    return;
  }
  /* else, valid bin - proceed */
  for(j=0; j < depo_ptr->pde_nbinDE; j++){
    for(k=0; k < depo_ptr->pde_nbinPz; k++){
      /* (j,k,iE,iPz,iMu) */
      /* this loop looks pretty wrong, lets make them all the same yeah? */
      ind = get_pdedp_ind(depo_ptr, iE, iPz, iMu, k, j);
      tmp = depo_ptr->pde_pdedp[ind];
      ptot = ptot + tmp;
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
  int k, kd, np2, nprt0;
  double dum, xproj, tdum;

  const double einj1 = depo_ptr->pde_Emin;  /* [keV] */
  const double einj2 = depo_ptr->pde_Emax;
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
    en[kd] = (einj1 + rand_double()*(einj2-einj1)) * engn/ekev; /* kinetic energy */
  }

  cfg_ptr->nprt = kd;
  nprt0 = kd;
  printf("%d fulldepmp_co deposit\n", nprt0);
  for(k=0; k < cfg_ptr->nprt; k++){
    kfield(cfg_ptr, kd);
    rho[k] = ptch[k]*sqrt(2.*en[k])/b[k];
    rmu[k] = en[k]/b[k] - .5*rho[k]*rho[k]*b[k];
    en[k] = en[k] + pot[k];
  }

  /* re-sample lost particles */
  class_domain(cfg_ptr);
  fullredepmp(cfg_ptr, depo_ptr);   /* goto 11 */

  return;
}
