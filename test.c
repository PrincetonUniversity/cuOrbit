
#include "orbit_config.h"
#include "orbit_equilibrium.h"
#include "orbit_perturbation.h"
#include "orbit_particles.h"

// test
int main(){

  Config_t* Cfg=Config_ctor();
  initialize_Config(Cfg);

  /* these will be managed by the config soon */
  Equilib_t* Eq=Equilib_ctor();
  initialize_Equilib(Eq);

  Particles_t* Ptcl=Particles_ctor();
  initialize_Particles(Ptcl);

  Perturb_t* Ptrb=Perturb_ctor();
  initialize_Perturb(Ptrb, Cfg, Eq, Ptcl);
}
