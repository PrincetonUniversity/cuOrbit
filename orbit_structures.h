#ifndef SET_ORBIT_STRUCTURES_H_
#define SET_ORBIT_STRUCTURES_H_

extern const int IDM;
extern const int IDP;
extern const int IDT;
extern const int NTOR;

typedef struct Config {

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

typedef struct Particle {
  double chrg;  /* 1 for ion, -1 for electron */
  double zprt;
  double prot;  /* proton mass in proton units */
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
} Particle_t;



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
  double *xx;
  double *yy;
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


void initialize_config(Config_t*);
void initialize_particle(Particle_t*);
void initialize_Perturb(Perturb_t*);



#endif
