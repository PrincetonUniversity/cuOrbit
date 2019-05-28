#ifndef SET_ORBIT_STRUCTURES_H_
#define SET_ORBIT_STRUCTURES_H_

const int IDM=2000;
const int IDP=210;
const int IDT=150;
const int NTOR=5000;


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


/* this is generally speaking the values taken in from "spdata" */
typedef struct Equilib {
  double pw;
  double ped;
  double pwd;
  int lsp;
  int lst;
  int lemax;
  int lrmax;
  int krip;
  int  nrip;
  double rmaj;
  double d0;
  double brip;
  double wrip;
  double xrip;

  /* arrays */
  // we can make this into a single buffer with "smarter" indexes
  // once we are confident code works... for now might help debug etc.
  double *b1, *b2, *b3, *b4, *b5, *b6, *b7, *b8, *b9;
  double *g1, *g2, *g3, *g4, *g5, *g6, *g7, *g8, *g9;
  double *r1, *r2, *r3, *r4, *r5, *r6, *r7, *r8, *r9;
  double *x1, *x2, *x3, *x4, *x5, *x6, *x7, *x8, *x9;
  double *z1, *z2, *z3, *z4, *z5, *z6, *z7, *z8, *z9;
  double *gd1, *gd2, *gd3;
  double *pd1, *pd2, *pd3;
  double *ps1, *ps2, *ps3;
  double *qd1, *qd2, *qd3;
  double *rd1, *rd2, *rd3;
  double *rp1, *rp2, *rp3;

  //xxx i think this goes into the perturb struct?
  int NA_;  /* 155 */
  double *amp;
} Equilib_t;

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
void initialize_Equilib(Equilib_t*);



#endif
