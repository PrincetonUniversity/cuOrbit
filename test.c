
#include "orbit_structures.h"
#include "orbit_equilibrium.h"
#include "orbit_perturbation.h"
#include "orbit_particles.h"

// test
int main(){

  Config_t Cfg;
  initialize_Config(&Cfg);

  Equilib_t Eq;
  initialize_Equilib(&Eq);

  Particles_t Ptcl;
  initialize_Particles(&Ptcl, &Cfg);

  Perturb_t Ptrb;
  initialize_Perturb(&Ptrb, &Cfg, &Eq, &Ptcl);
}
