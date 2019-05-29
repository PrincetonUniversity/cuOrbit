#ifndef SET_ORBIT_PARTICLES_H_
#define SET_ORBIT_PARTICLES_H_

#include "orbit_structures.h"

//not sure if this is better way yet
/* typedef struct Particle { */
/*   double chrg;  /\* 1 for ion, -1 for electron *\/ */
/*   double zprt; */
/*   double prot;  /\* proton mass in proton units *\/ */

/* } Particle_t; */

typedef struct Particles {
  /* Particle_t* particle; */
  double chrg;  /* 1 for ion, -1 for electron */
  double zprt;
  double prot;  /* proton mass in proton units */

  size_t idm;
  int *otp;
  double *pol;
  double *zet;
  double *thet;
  double *rho;
  double *en;
  double *rmu;
  double *ptch;
  double *pot;
  double *time;  /* time step */
  double *g;
  double *gp;
  double *q;
  double *qp;
  double *b;  /* B, I associated with particle */
  double *ri;
  double *rip;
  double *w1, *w2, *w3;  /* particle stepping */
} Particles_t;

void initialize_Particles(Particles_t*, Config_t*);

#endif
