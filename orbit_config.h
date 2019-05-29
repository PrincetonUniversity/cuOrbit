#ifndef SET_ORBIT_CONFIG_H_
#define SET_ORBIT_CONFIG_H_

/* these are set in the .c file for now */
extern const int IDP;
extern const int IDT;
extern const int NTOR;

typedef struct Config {

  int nprt;   /* number of particles (use in place of IDM*/
  int nmds;  /* probably dont need this */
  int nrmds;  /* dont remember this */
  int seed;   /* used by the RNG */
  double ekev;
  double engn;
  double bkg;
  double bmin;
  double bmax;
  double rmaj;
  double trun;
  double tran;  /* transit time for particle at the mag axi with pitch=1 */
  double dt0;
  double omeg0;
  double xc;
  double eps;
  double bax;
  double *dwal;  /* wall definition */

  double pamp;
  double rprof;

} Config_t;

#include <stdlib.h>



void initialize_Config(Config_t*);


#endif
