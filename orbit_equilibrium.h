#ifndef SET_ORBIT_EQUILIBRIUM_H_
#define SET_ORBIT_EQUILIBRIUM_H_

#include "orbit_config.h"

/* this is generally speaking the values taken in from "spdata" */
typedef struct Equilib Equilib_t;

void initialize_Equilib(Equilib_t*, Config_t*);
Equilib_t* Equilib_ctor();

double** get_B(Equilib_t*);
double** get_G(Equilib_t*);
double** get_R(Equilib_t*);
double** get_X(Equilib_t*);
double** get_Z(Equilib_t*);


double** get_QD(Equilib_t*);
double** get_RD(Equilib_t*);
double** get_GD(Equilib_t*);
double** get_PD(Equilib_t*);
double** get_PS(Equilib_t*);
double** get_RP(Equilib_t*);
double** get_VD(Equilib_t*);

double get_pw(Equilib_t*);
double get_rmaj(Equilib_t*);
double get_ped(Equilib_t*);
double get_lsp(Equilib_t*);
double get_lst(Equilib_t*);
double get_nrip(Equilib_t*);


/* utils */
int compute_jd(Equilib_t*, double);
double gfun(Equilib_t*, double);
double qfun(Equilib_t*, double);
double rifun(Equilib_t*, double);
double bfield(Equilib_t*, double, double);
double giac(Equilib_t*, double, double);
void vspline(Equilib_t*);

double xproj(Equilib_t*, double, double);
double zproj(Equilib_t*, double, double);

#endif
