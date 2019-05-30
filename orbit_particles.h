#ifndef SET_ORBIT_PARTICLES_H_
#define SET_ORBIT_PARTICLES_H_

#include <stdlib.h>

//not sure if this is better way yet
/* typedef struct Particle { */
/*   double chrg;  /\* 1 for ion, -1 for electron *\/ */
/*   double zprt; */
/*   double prot;  /\* proton mass in proton units *\/ */

/* } Particle_t; */

typedef struct Particles Particles_t;

void initialize_Particles(Particles_t**);

double* get_pol(Particles_t* ptcl_ptr);

#endif
