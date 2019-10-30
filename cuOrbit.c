#include <stdio.h>
#include "orbit_config_api.h"


// test
int main(){

  /* Contruct a config object,
   which conttructs Perturb, Equlib, and Particles */
  orbit_Config_t* Cfg=orbit_Config_ctor();

  printf("Initialize Config\n");
  orbit_initialize_Config(Cfg);

  orbit_main_loop(Cfg);

  return 0;
}
