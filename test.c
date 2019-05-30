
#include "orbit_config.h"
#include "orbit_equilibrium.h"
#include "orbit_perturbation.h"
#include "orbit_particles.h"

// test
int main(){

  Config_t* Cfg=NULL;
  initialize_Config(Cfg);

  /* these will be managed by the config soon */
  Equilib_t* Eq=NULL;
  initialize_Equilib(&Eq);

  Particles_t* Ptcl=NULL;
  initialize_Particles(&Ptcl);

  Perturb_t* Ptrb=NULL;
  initialize_Perturb(&Ptrb, Cfg, Eq, Ptcl);
}
