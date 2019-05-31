
#include "orbit_config_api.h"


// test
int main(){

  /* Contruct a config object,
   which conttructs Perturb, Equlib, and Particles */
  Config_t* Cfg=Config_ctor();
  initialize_Config(Cfg);

  /* these will be managed by the config soon */
}
