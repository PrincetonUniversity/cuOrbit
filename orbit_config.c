#include <stdlib.h>
#include <stdio.h>
#include "orbit_config.h"

const int IDP=210;
const int IDT=150;
const int NTOR=5000;


void initialize_Config(Config_t* cfg_ptr){

  //Config_t cfg = *cfg_ptr;

  cfg_ptr->nprt = 2000;
  cfg_ptr->seed = 12345;
  cfg_ptr->ekev = 75;
  /* engn = ; */
  cfg_ptr->bkg = 5.5;
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

  cfg_ptr->pamp = 8.63;
  cfg_ptr->rprof = 3.765;

};
