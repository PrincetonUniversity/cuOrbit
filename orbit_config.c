#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "inih/ini.h"

#include "orbit_config_api.h"


const int IDP=210;
const int IDT=150;
const int NTOR=5000;



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
    } else if (MATCH("model", "ekev")) {
        pconfig->ekev = atof(value);
    } else if (MATCH("model", "bkg")) {
        pconfig->bkg = atof(value);
    } else if (MATCH("model", "pamp")) {
        pconfig->pamp = atof(value);
    } else if (MATCH("model", "rprof")) {
        pconfig->rprof = atof(value);
    /* } else if (MATCH("model", "XXX")) { */
    /*     pconfig->xxx = atof(value); */

        /* Perturb Config Vars*/
    } else if (MATCH("model", "omeg0")) {
      pconfig->omeg0 = atof(value);
      
      /* Particles Config Vars*/
    } else if (MATCH("model", "nprt")) {
        pconfig->nprt = atoi(value);
    } else {
      return 0;  /* unknown section/name, error */
    }
    return 1;
}

static int config_file_handler(char* config_fname, Config_t* config){
  if (ini_parse(config_fname, config_handler, &config) < 0) {
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
  /* double omeg0; */
  /* double xc; */
  /* double eps; */
  /* double bax; */
  /* double* dwal;  /\* wall definition *\/ */


  config_file_handler("INPUT/config.ini", cfg_ptr);
  printf("rx'd nprt %d seed %d\n", cfg_ptr->nprt, cfg_ptr->seed);

  /* initialize the other model components */
  initialize_Equilib(cfg_ptr->eqlb_ptr);
  initialize_Particles(cfg_ptr->ptcl_ptr);
  initialize_Perturb(cfg_ptr->ptrb_ptr, cfg_ptr, cfg_ptr->eqlb_ptr, cfg_ptr->ptcl_ptr);


}
