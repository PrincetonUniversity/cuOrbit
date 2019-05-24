#include "orbit_structures.h"


void initialize_config(Config_t* cfg_ptr){

  Config_t cfg = *cfg_ptr;

  cfg.seed = 12345;
  cfg.ekev = 75;
  /* engn = ; */
  cfg.bkg = 5.5;
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

  cfg.pamp = 8.63;
  cfg.rprof = 3.765;

};

void initialize_particle(Particle_t* ptc_ptr){
  
};

void initialize_Perturb(Perturb_t* ptrb_ptr){

};
