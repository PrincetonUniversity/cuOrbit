#include <stdlib.h>
#include <stdio.h>

#include "orbit_particles.h"

typedef struct Particles {
  int nprt;   /* number of particles (use in place of IDM*/
  double chrg;  /* 1 for ion, -1 for electron */
  double zprt;
  double prot;  /* proton mass in proton units */

  size_t idm;
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
} Particles_t;


void initialize_Particles(Particles_t** ptcl_ptr_ptr){
  *ptcl_ptr_ptr = (Particles_t*)calloc(1, sizeof(Particles_t));
  Particles_t* ptcl_ptr = *ptcl_ptr_ptr;

  /* from the config we can get number of particles
     which we will use to dimenion our arrays
  */
  ptcl_ptr->idm = (unsigned)ptcl_ptr->nprt;

  /* we'll code these in for now, and admit a file+struct based config for parameters soon */

  ptcl_ptr->chrg=1.;  /* 1 for ion, -1 for electron */
  ptcl_ptr->zprt=1.;
  ptcl_ptr->prot=2.;  /* proton mass in proton units */

  /* arrays */
  ptcl_ptr->otp = (int*)calloc(ptcl_ptr->idm, sizeof(int));

  ptcl_ptr->pol=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->zet=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->thet=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->rho=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->en=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->rmu=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->ptch=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->pot=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->time=(double*)calloc(ptcl_ptr->idm, sizeof(double));  /* time step */
  ptcl_ptr->g=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->gp=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->q=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->qp=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->b=(double*)calloc(ptcl_ptr->idm, sizeof(double));  /* B, I associated with particle */
  ptcl_ptr->ri=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->rip=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->w1=(double*)calloc(ptcl_ptr->idm, sizeof(double)); /* particle stepping */
  ptcl_ptr->w2=(double*)calloc(ptcl_ptr->idm, sizeof(double));
  ptcl_ptr->w3=(double*)calloc(ptcl_ptr->idm, sizeof(double));



}

double* get_pol(Particles_t* ptcl_ptr){
  return ptcl_ptr->pol;
}
