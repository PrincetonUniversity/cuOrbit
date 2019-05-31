#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "inih/ini.h"

#include "orbit_config_api.h"


const int IDP=210;
const int IDT=150;



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
    else if (MATCH("perturbation", "nprt")) {
        pconfig->nprt = atoi(value);
    } else if (MATCH("perturbation", "zprt")) {
      pconfig->zprt = atof(value);
    } else if (MATCH("perturbation", "chrg")) {
      pconfig->chrg = atof(value);
    } else if (MATCH("perturbation", "prot")) {
      pconfig->prot = atof(value);
    } else if (MATCH("perturbation", "ekev")) {
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

  /* engn = ; */
  /* double bmin; */
  /* double bmax; */
  /* double rmaj; */
  /* double trun; */
  /* double tran;  /\* transit time for particle at the mag axi with pitch=1 *\/ */
  /* double dt0; */
  /* double xc; */
  /* double eps; */
  /* double bax; */
  /* double* dwal;  /\* wall definition *\/ */


  config_file_handler("INPUT/config.ini", cfg_ptr);
  printf("rx'd nprt %d seed %d\n", cfg_ptr->nprt, cfg_ptr->seed);
  /* initialize the other model components */
  initialize_Equilib(cfg_ptr->eqlb_ptr, cfg_ptr);

  set1(cfg_ptr);

  initialize_Particles(cfg_ptr->ptcl_ptr, cfg_ptr);

  initialize_Perturb(cfg_ptr->ptrb_ptr, cfg_ptr, cfg_ptr->eqlb_ptr, cfg_ptr->ptcl_ptr);


}

void set1(Config_t* cfg_ptr){

  
  return;
}
