#ifndef SET_ORBIT_PERTURBATION_H_
#define SET_ORBIT_PERTURBATION_H_

#include "orbit_config.h"
#include "orbit_particles.h"
#include "orbit_equilibrium.h"

/* we can pack this in the struct if we really need it... */
extern const int NAMP_;

typedef struct Perturb Perturb_t;

void initialize_Perturb(Perturb_t*, Config_t*, Equilib_t*, Particles_t*);
Perturb_t* Perturb_ctor();

void splna(Perturb_t* , Equilib_t*, Particles_t*);
void splnx(Perturb_t* , Equilib_t*, Particles_t*);

#endif
