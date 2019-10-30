#include <stdio.h>
#include "orbit_config_api.h"


// test
int main(){

  /* Contruct a config object,
   which conttructs Perturb, Equlib, and Particles */
  Config_t* Cfg=Config_ctor();

  printf("Initialize Config\n");
  initialize_Config(Cfg);

  main_loop(Cfg);

  return 0;
}
