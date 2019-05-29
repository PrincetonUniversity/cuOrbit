#include "orbit_structures.h"
#include "orbit_equilibrium.h"
#include "orbit_perturbation.h"

// test
int main(){

  Config_t Cfg;
  initialize_Config(&Cfg);

  Equilib_t Eq;
  initialize_Equilib(&Eq);

  Particle_t Ptcl;
  initialize_Particle(&Ptcl);

  Perturb_t Ptrb;
  initialize_Perturb(&Ptrb, &Cfg, &Eq, &Ptcl);
}
