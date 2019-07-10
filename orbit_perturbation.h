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
void set_omeg0(Perturb_t*, double);
double get_omeg0(Perturb_t*);
int get_nflr(Perturb_t*);
int get_npert(Perturb_t*);
int get_lpt(Perturb_t*);
int get_md1(Perturb_t*);
int get_md2(Perturb_t*);

double* get_phaz(Perturb_t*);
double* get_omegv(Perturb_t*);

double* get_alfv(Perturb_t*);
double* get_amp(Perturb_t*);
int* get_mmod(Perturb_t*);
int* get_nmod(Perturb_t*);
double* get_xi1(Perturb_t*);
double* get_xi2(Perturb_t*);
double* get_xi3(Perturb_t*);
double* get_a1(Perturb_t*);
double* get_a2(Perturb_t*);
double* get_a3(Perturb_t*);

double pol2pot(Config_t*, double);
#endif
