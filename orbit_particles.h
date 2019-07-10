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

void initialize_Particles(Particles_t*, Config_t*);
Particles_t* Particles_ctor();

double* get_b(Particles_t*);
double* get_g(Particles_t*);
double* get_q(Particles_t*);
double* get_en(Particles_t*);
double* get_pol(Particles_t*);
double* get_rho(Particles_t*);
double* get_rmu(Particles_t*);
int* get_otp(Particles_t*);
double* get_ptch(Particles_t*);
double* get_thet(Particles_t*);
double* get_pot(Particles_t*);
double* get_time(Particles_t*);
double* get_dt(Particles_t*);
double* get_tim1(Particles_t*);
double* get_wt(Particles_t*);
double* get_zet(Particles_t*);

double get_zprt(Particles_t*);
double get_prot(Particles_t*);
double get_prot(Particles_t*);
double get_ekev(Particles_t*);

double* get_dadp(Particles_t*);
double* get_dadt(Particles_t*);
double* get_dadz(Particles_t*);
double* get_padt(Particles_t*);


void field(Config_t*, int);
void kfield(Config_t*, int);

void ptrbak(Config_t*, int);
void ptrb2k(Config_t*, int);

void do_particles(Config_t* );
void konestep(Config_t*, int);


#endif
