#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "inih/ini.h"

#include "orbit_config.h"
#include "orbit_equilibrium.h"
#include "orbit_perturbation.h"
#include "orbit_particles.h"

const int IDP=210;
const int IDT=150;
const int NTOR=5000;



typedef struct Config {
  /* meta */
  char* name;

  /* model */
  int nmds;  /* probably dont need this */
  int nrmds;  /* dont remember this */
  int seed;   /* used by the RNG */
  double ekev;
  double engn;
  double bkg;
  double bmin;
  double bmax;
  double rmaj;
  double trun;
  double tran;  /* transit time for particle at the mag axi with pitch=1 */
  double dt0;
  double xc;
  double eps;
  double bax;
  double *dwal;  /* wall definition */
  double pamp;
  double rprof;

  /* inputs */
  char* spdata_file;
  char* alphas_file;
  char* displ_file;
  /* outputs */
  char* pdedp_file;

  Equilib_t* eqlb_ptr;

  Perturb_t* ptrb_ptr;
  double omeg0;

  /* particle */
  Particles_t* ptcl_ptr;
  int nprt;


} Config_t;


static int config_handler(void* res_ptr, const char* section, const char* name,
                   const char* value){
    Config_t* pconfig = (Config_t*)res_ptr;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    /* MATCH takes a section and configuration variable name,
       returns associated string.  We control the parsing. */
    if (MATCH("model", "name")) {
        pconfig->name = strdup(value);
    } else if (MATCH("model", "nprt")) {
        pconfig->nprt = atoi(value);
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
    } else {
        return 0;  /* unknown section/name, error */
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
  cfg_ptr = (Config_t*)calloc(1, sizeof(Config_t));

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

}
