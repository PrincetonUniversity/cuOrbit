#ifndef SET_ORBIT_PERTURBATION_H_
#define SET_ORBIT_PERTURBATION_H_

#include "orbit_structures.h"
#include "orbit_equilibrium.h"

/* we can pack this in the struct if we really need it... */
extern const int NAMP_;

typedef struct Perturb {
  int npert;
  int nflr;
  int lpt;
  int md1;
  int md2;
  int modes;
  double psi_solRZ;
  int *mmod;  /* pol mode numbers */
  int *nmod;  /* pol mode numbers */
  double *omegv;  /* mode frequency */
  double *amp;    /* amp */
  double *damp;
  double *harm;
  double *xx;
  //unused? double *yy;
  double *alfv;
  /* derivatives */
  double *dbdt;
  double *dbdp;
  double *dbdz;
  double *alp;
  double *dadp;
  double *dadt;
  double *dadz;
  double *padt;
  double *phaz;
  /* displacements xi1-xi3 */
  double *xi1, *xi2, *xi3;
  /* alphas a1-a3*/
  double *a1, *a2, *a3;
} Perturb_t;

void initialize_Perturb(Perturb_t*, Config_t*, Equilib_t*, Particle_t*);

void splna(Perturb_t* , Equilib_t*, Particle_t*);
void splnx(Perturb_t* , Equilib_t*, Particle_t*);

#endif
