#ifndef SET_ORBIT_EQUILIBRIUM_H_
#define SET_ORBIT_EQUILIBRIUM_H_

/* this is generally speaking the values taken in from "spdata" */
typedef struct Equilib Equilib_t;

void initialize_Equilib(Equilib_t**);

double get_pw(Equilib_t*);


/* utils */
double gfun(Equilib_t*, double);
double qfun(Equilib_t*, double);
double rifun(Equilib_t*, double);


#endif
