#ifndef SET_ORBIT_CONFIG_API_H_
#define SET_ORBIT_CONFIG_API_H_

#include <stdbool.h>
#include "orbit_config.h"
#include "orbit_perturbation.h"
#include "orbit_equilibrium.h"
#include "orbit_particles.h"
#include "orbit_deposition.h"

/* these are set in the .c file for now */
extern const int IDP;
extern const int IDT;
extern const double pi2;
extern const double pi2i;


typedef struct Config {
  /* meta */
  char* name;

  /* model */
  int nmds;  /* probably dont need this */
  int nrmds;  /* dont remember this */
  int seed;   /* used by the RNG */
  int ntor;
  int nstep_all;
  int nskip;   //xxx i suspect this will go to another struct...
  double bsum; //xxx i suspect this will go to another struct...
  double dsum; //xxx i suspect this will go to another struct...
  double esum; //xxx i suspect this will go to another struct...
  double pvol;  /* stash plasma volume, another struct? */
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
  char* bfield_file;

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

  /* Deposition*/
  Deposition_t* depo_ptr;
  /* xxx stochastic, does this belong here? */
  double mubk_scale;
  int emink;
  int emaxk;
  int dmubk;
  int nstoche;
  int nstochp;
  /* pdedp */
  int nruns;
  bool compute_pdedp;
  bool initial_update_pdedp;
  double deposit_on_bins_after_fraction;
  double pdedp_dtsamp;
  double pdedp_dtav;
  int pdedp_tskip;
  double pdedp_otpup;
  bool pdedp_focusdep;
  bool pdedp_optimize;

} Config_t;


Config_t* Config_ctor();
void initialize_Config(Config_t*);

void set_xc(Config_t*, double);
double get_xc(Config_t*);
double get_bkg(Config_t*);
void set_eps(Config_t*, double);
double get_eps(Config_t*);
double get_engn(Config_t*);
double get_bax(Config_t*);
double get_bmax(Config_t*);
double get_bmin(Config_t*);
double get_pamp(Config_t*);
double get_rprof(Config_t*);
double get_trun(Config_t*);
int get_nstep_all(Config_t*);
double get_dt0(Config_t*);

#endif
