#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "inih/ini.h"

#include "orbit_config_api.h"
#include "orbit_util.h"

const int IDP=210;
const int IDT=150;
const double pi2 = 2. * M_PI;
const double pi2i = 1. / pi2;


Config_t* Config_ctor(){
  Config_t* Cfg = (Config_t*)calloc(1, sizeof(Config_t));
    /* Create the other model componenets */
  Cfg->eqlb_ptr = Equilib_ctor();
  Cfg->ptcl_ptr = Particles_ctor();
  Cfg->ptrb_ptr = Perturb_ctor();
  return Cfg;
}


static int config_handler(void* res_ptr, const char* section, const char* name,
                   const char* value){
    Config_t* pconfig = (Config_t*)res_ptr;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    /* MATCH takes a section and configuration variable name,
       returns associated string.  We control the parsing. */
    if (MATCH("model", "name")) {
      pconfig->name = strdup(value);
    } else if (MATCH("model", "rng_seed")) {
        pconfig->seed = atoi(value);
    } else if (MATCH("model", "ntor")) {
        pconfig->ntor = atof(value);
    } else if (MATCH("model", "bkg")) {
        pconfig->bkg = atof(value);
    } else if (MATCH("model", "pamp")) {
        pconfig->pamp = atof(value);
    } else if (MATCH("model", "rprof")) {
        pconfig->rprof = atof(value);
    }
    /* Inputs */
    else if (MATCH("input", "spdata_file")) {
      pconfig->spdata_file = strdup(value);
    } else if (MATCH("input", "alphas_file")) {
      pconfig->alphas_file = strdup(value);
    } else if (MATCH("input", "displ_file")) {
      pconfig->displ_file = strdup(value);
    }
    /* Output */
    else if (MATCH("output", "pdedp_file")) {
      pconfig->pdedp_file = strdup(value);
    }

    /* Perturb Config Vars*/
    else if (MATCH("perturbation", "falf")) {
      pconfig->falf = atof(value);
    } else if (MATCH("perturbation", "ascale")) {
      pconfig->ascale = atof(value);
    } else if (MATCH("perturbation", "alimit")) {
      pconfig->alimit = atof(value);
    } else if (MATCH("perturbation", "global_scaling_factor")) {
      pconfig->global_scaling_factor = atof(value);
    } else if (MATCH("perturbation", "freq_scaling_factor")) {
      pconfig->freq_scaling_factor = atof(value);
    }
    /* Particles Config Vars*/
    else if (MATCH("particle", "nprt")) {
        pconfig->nprt = atoi(value);
    } else if (MATCH("particle", "zprt")) {
      pconfig->zprt = atof(value);
    } else if (MATCH("particle", "chrg")) {
      pconfig->chrg = atof(value);
    } else if (MATCH("particle", "prot")) {
      pconfig->prot = atof(value);
    } else if (MATCH("particle", "ekev")) {
      pconfig->ekev = atof(value);
    } else if (MATCH("particle_distribution", "polo_scale")) {
      pconfig->polo_scale = atof(value);
    } else if (MATCH("particle_distribution", "p1_scale")) {
      pconfig->p1_scale = atof(value);
    } else if (MATCH("particle_distribution", "p2_scale")) {
      pconfig->p2_scale = atof(value);
    } else if (MATCH("particle_distribution", "pchi")) {
      pconfig->pchi = atof(value);
    }

    /* unknown section/name, error */
    else {
      return 0;
    }
    return 1;
}

static int config_file_handler(char* config_fname, Config_t* config){
  if (ini_parse(config_fname, config_handler, config) < 0) {
    fprintf(stderr, "Can't load '%s'\n", config_fname);
    return 1;
  }
  printf("Config loaded from '%s'\n", config_fname);
  return 0;
}


void initialize_Config(Config_t* cfg_ptr){


  config_file_handler("INPUT/config.ini", cfg_ptr);
  printf("rx'd nprt %d seed %d\n", cfg_ptr->nprt, cfg_ptr->seed);
  /* initialize the other model components */
  initialize_Equilib(cfg_ptr->eqlb_ptr, cfg_ptr);

  initialize_Particles(cfg_ptr->ptcl_ptr, cfg_ptr);

  set1(cfg_ptr);

  //xxxx I need to check if we need anything in set1 for particles to be correct
  //xxxxx set1 expecs pol initialize_Particles(cfg_ptr->ptcl_ptr, cfg_ptr);

  initialize_Perturb(cfg_ptr->ptrb_ptr, cfg_ptr, cfg_ptr->eqlb_ptr, cfg_ptr->ptcl_ptr);


}

void set1(Config_t* cfg_ptr){
  int j, k;
  int kdum;
  double dum;
  double pdum, qdum, tdum;
  double q0, qw;
  Equilib_t* Eq = cfg_ptr->eqlb_ptr;
  Perturb_t* Ptrb = cfg_ptr->ptrb_ptr;
  Particles_t* Ptcl = cfg_ptr->ptcl_ptr;
  double * const q = get_q(Ptcl);
  const int nprt = cfg_ptr->nprt;

  double psiwal;
  double pq1;
  double rq1;
  double pdp, pdq, pqm, pqp, qdm, qdp, rqm, rqp;
  double denom, bdum, rdum, ddum;
  double pz, tz;
  double xdum, zdum;
  double xl, xr, zb, zt, rmd;

  for(k=1; k<=1000; k++){
    pdum = 0.001 * k * get_pw(Eq);
    dum = dum + qfun(Eq, pdum);
  }
  psiwal = dum;
  printf("\t Toroidal psi wall %f\n", psiwal);

  set_omeg0(Ptrb, 9.58E6 * get_zprt(Ptcl) * cfg_ptr->bkg / get_prot(Ptcl));

  set_xc(cfg_ptr, xproj(Eq, 0., 0.));
  set_eps(cfg_ptr, ( xproj(Eq, get_ped(Eq) , 0.) - xproj(Eq, 0., 0.)) / get_xc(cfg_ptr) );

  cfg_ptr->bmin = bfield(Eq, get_pw(Eq), 0.);
  cfg_ptr->bmax = bfield(Eq, get_pw(Eq), M_PI);
  cfg_ptr->bax = bfield(Eq, 0., 0.);

  dum = fabs(cfg_ptr->bax - 1.);
  /* Sanity check */
  if (dum > 5E-3){
    fprintf(stderr, "Equlib is improperly initialized %f\n", dum);
    exit(1);
  }

  double* pol = get_pol(Ptcl);
  double* thet = get_thet(Ptcl);
  double* pot = get_pot(Ptcl);
  double* time = get_time(Ptcl);
  double* dt = get_dt(Ptcl);
  double* tim1 = get_tim1(Ptcl);
  double* wt = get_wt(Ptcl);

  double** const B = get_B(Eq);
  double** const QD = get_QD(Eq);
  double** const GD = get_GD(Eq);
  double** const RD = get_RD(Eq);
  pol[0] = 1E-10;
  thet[0] = 0;

  /* xxx not sure why exactly we do this... */
  field(cfg_ptr, Ptcl, 1);
  /* or this */
  q0 = q[0];
  pol[0] = get_pw(Eq);
  field(cfg_ptr, Ptcl, 1);
  qw = q[0];
  for(k=1; k<=200; k++){
    pol[k] = 0.005 * k * get_pw(Eq);
    thet[k] = 0;
  }
  field(cfg_ptr, Ptcl, 200);
  pq1 = 9.99E9;
  rq1 = 9.99E9;
  if (darray_max(q, (size_t)nprt)>1. || darray_min(q, (size_t)nprt)<1.){
    for(k=1; k<=200; k++){
      kdum = k;
      pq1 = pol[k-1];
      if(q[k] > 1.){
        break;
      }
    }
    pqp = 0.005 * kdum * get_pw(Eq);
    pqm = 0.005 * (kdum -1) * get_pw(Eq);
    rqp = xproj(Eq, pqp, 0.) - get_xc(cfg_ptr);
    rqm = xproj(Eq, pqm, 0.) - get_xc(cfg_ptr);
    qdp = q[kdum];
    if(kdum == 1){
      rq1 = xproj(Eq, pq1, 0.) - get_xc(cfg_ptr);
    }else{
      qdm = q[kdum-1];
      if(kdum == 1){
        qdum = q0;
      }
      rq1 = rqp - (rqp-rqm)*(qdp - 1.)/(qdp - qdm);
      pq1 = pqp - (pqp-pqm)*(qdp - 1.)/(qdp - qdm);
    }
  }
  cfg_ptr->engn = 10.533 * get_prot(Ptcl) * get_ekev(Ptcl) * pow(GD[0][0], 2) /
      pow( get_rmaj(Eq) * get_zprt(Ptcl) * cfg_ptr->bkg * B[0][0], 2);

  printf("set1 %f %f %f %f %f %f %f %f\n",
         get_prot(Ptcl), get_ekev(Ptcl), QD[0][0], get_rmaj(Eq),
         get_zprt(Ptcl), cfg_ptr->bkg, B[0][0], cfg_ptr->engn);

  denom = B[0][0] * sqrt(2. * cfg_ptr->engn) * QD[0][0];
  cfg_ptr->tran = 6.28 * (GD[0][0] * QD[0][0] + RD[0][0])/denom;
  cfg_ptr->trun = cfg_ptr->ntor * cfg_ptr->tran;
  cfg_ptr->dt0 = cfg_ptr->tran/200;
  cfg_ptr->nskip = cfg_ptr->trun/
      (200 * cfg_ptr->dt0) + 1;   /* for diffusion and bootstrap */

  cfg_ptr->bsum = 0.;
  cfg_ptr->dsum = 0.;
  cfg_ptr->esum = 0.;
  for(k=0; k<nprt; k++){
    dt[k] = cfg_ptr->dt0;
    tim1[k] = 0.;
    time[k] = 0.;
    pot[k] = 0.;
    wt[k] = 0.;
  }

  /*   find plasma volume, pdum = pw for total
       units of numerical equilibrium  */
  pdum = get_pw(Eq);
  cfg_ptr->pvol = 0.;
  const int nint0 = 100;
  for(k=1; k<=nint0; k++){
    for(j=1; j<=nint0; j++){
      pz = get_pw(Eq) - (j - .5) * pdum/((double)nint0);
      tz = k * 2 * M_PI / ((double)nint0);
      cfg_ptr->pvol += pdum * 2 * M_PI * giac(Eq, pz, tz) / (
          (double)(nint0*nint0));
    }  /* j */
  }    /* k */
  cfg_ptr->pvol *=  2 * M_PI;  /* 2*pi*R */
  printf("Plasma Volume %f\tbax %f\n", cfg_ptr->pvol, cfg_ptr->bax);

  /* /\* setup for plots, commenting for now *\/ */
  /* const int nbinx = 50; */
  /* double pop[nbinx], xpop[nbinx], pv[nbinx], zv[nbinx]; */
  /* int nbin; */
  /* for(nbin=1; nbin<=nbinx; nbin++){ */
  /*   pop[nbin] = 0.; */
  /*   xpop[nbin] = 0.; */
  /*   pv[nbin] = 0.; */
  /*   zv[nbin] = 0.; */
  /* }  /\* nbin *\/ */

  /* some sort of checking going on here */
  xl = 1.E6;
  xr = -1.E6;
  zb = 1.E6;
  zt = -1.E6;
  pz = 0.99 * get_pw(Eq);
  for(k=1; k<=200; k++){
    tdum = .0314 * k;
    xdum = xproj(Eq, pz, tdum) * get_rmaj(Eq) / get_xc(cfg_ptr);
    xl = fmin(xl, xdum);
    xr = fmax(xr, xdum);
    zdum = zproj(Eq, pz, tdum) *  get_rmaj(Eq) / get_xc(cfg_ptr);
    zb = fmin(zb, zdum);
    zt = fmax(zt, zdum);
  }  /* k */

  /* check vanishing of D at axis */
  const double gdum = GD[0][0];
  const double gdump = GD[1][0];
  const double ridum = RD[0][0];
  const double ridump = RD[1][0];
  qdum = QD[0][0];
  const double qdump = QD[1][0];
  bdum = cfg_ptr->bax;
  bool invalid=false;
  for(j=1; j<=10; j++){
    rmd = j * cfg_ptr->engn / (10 * bdum);
    dum = 2 * cfg_ptr->engn -
        2 * rmd * bdum;
    if(dum < 0) continue;
    rdum = sqrt(dum) / bdum;
    ddum = gdum * qdum + ridum + rdum * (gdum*ridump - ridum*gdump);
    if(dum < 0){
      invalid=true;
    };
    ddum = gdum*qdum + ridum - rdum * (gdum*ridump - ridum*gdump);
    if(dum < 0){
      invalid=true;
    };
  }
  if(invalid){
    fprintf(stderr, "D vanishes for this energy, some uses of code invalid");
    exit(1);
  }

  /* volume spline */
  vspline(Eq);

  return;
}


void set_xc(Config_t* cfg_ptr, double val){
  cfg_ptr->xc = val;
}

double get_xc(Config_t* cfg_ptr){
  return cfg_ptr->xc;
}


void set_eps(Config_t* cfg_ptr, double val){
  cfg_ptr->xc = val;
}

double get_eps(Config_t* cfg_ptr){
  return cfg_ptr->xc;
}

