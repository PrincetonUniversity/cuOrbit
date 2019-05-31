#ifndef SET_ORBIT_CONFIG_API_H_
#define SET_ORBIT_CONFIG_API_H_

#include "orbit_config.h"
#include "orbit_perturbation.h"
#include "orbit_equilibrium.h"
#include "orbit_particles.h"

/* these are set in the .c file for now */
extern const int IDP;
extern const int IDT;
extern const int NTOR;

typedef struct Config {
  /* meta */
  char* name;

  /* model */
  int nmds;  /* probably dont need this */
  int nrmds;  /* dont remember this */
  int seed;   /* used by the RNG */
  double engn;
  double bkg;
  double bmin;
  double bmax;
  double rmaj;
  double trun;
  double tran;  /* transit time for particle at the mag axi with pitch=1 */
  double dt0;
  double xc;
  double eps;
  double bax;
  double *dwal;  /* wall definition */
  double pamp;
  double rprof;

  /* inputs */
  char* spdata_file;
  char* alphas_file;
  char* displ_file;
  /* outputs */
  char* pdedp_file;

  Equilib_t* eqlb_ptr;

  /* Pertubation */
  Perturb_t* ptrb_ptr;
  double falf;
  double ascale;
  double alimit;
  double global_scaling_factor;
  double freq_scaling_factor;
  double sng;

  /* particle */
  Particles_t* ptcl_ptr;
  int nprt;
  double zprt;
  double chrg;
  double prot;
  double ekev;
  double polo_scale;
  double p1_scale;
  double p2_scale;
  double pchi;


} Config_t;


Config_t* Config_ctor();
void initialize_Config(Config_t*);

#endif
