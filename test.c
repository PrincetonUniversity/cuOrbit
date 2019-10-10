
#include "orbit_config_api.h"
#include "orbit_main.h"


// test
int main(){

  /* Contruct a config object,
   which conttructs Perturb, Equlib, and Particles */
  Config_t* Cfg=Config_ctor();
  initialize_Config(Cfg);

  main_loop(Cfg);

  /* these will be managed by the config soon */
}
