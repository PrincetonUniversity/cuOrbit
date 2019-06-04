#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
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
        pconfig->pamp = atof(value);
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


  /* double trun; */
  /* double tran;  /\* transit time for particle at the mag axi with pitch=1 *\/ */
  /* double dt0; */
  /* double bax; */
  /* double* dwal;  /\* wall definition *\/ */


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
  int k;
  int kdum;
  double dum;
  double pdum, qdum;
  double q0, qw;
  Equilib_t* Eq = cfg_ptr->eqlb_ptr;
  Perturb_t* Ptrb = cfg_ptr->ptrb_ptr;
  Particles_t* Ptcl = cfg_ptr->ptcl_ptr;
  double * const q = get_q(Ptcl);
  const int nprt = cfg_ptr->nprt;

  double psiwal;
  double pq1;
  double rq1;
  double pdp, pdq, pqm, pqp, qdm, qdp;
  double rqm, rqp;

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
  cfg_ptr->engn = 10.533 * get_prot(Ptcl) * get_ekev(Ptcl) * pow(get_GD(Eq)[0][0], 2) /
      pow( get_rmaj(Eq) * get_zprt(Ptcl) * cfg_ptr->bkg * get_B(Eq)[0][0], 2);

  printf("%f %f %f %f %f %f %f %f\n",
         get_prot(Ptcl), get_ekev(Ptcl), get_QD(Eq)[0][0], get_rmaj(Eq),
         get_zprt(Ptcl), cfg_ptr->bkg, get_B(Eq)[0][0], cfg_ptr->engn);

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

