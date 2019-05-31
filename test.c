
#include "orbit_config_api.h"
#include "orbit_equilibrium.h"
#include "orbit_perturbation.h"
#include "orbit_particles.h"

// test
int main(){

  /* Contruct a config object,
   which conttructs Perturb, Equlib, and Particles */
  Config_t* Cfg=Config_ctor();
  printf("XXX Cfg->seed=%d\n", Cfg->seed);
  initialize_Config(Cfg);

  /* these will be managed by the config soon */
}
