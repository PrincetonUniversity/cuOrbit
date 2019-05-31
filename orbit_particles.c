#include <stdlib.h>
#include <stdio.h>

#include "orbit_config_api.h"
#include "orbit_particles.h"

typedef struct Particles {
  int nprt;   /* number of particles (use in place of IDM*/
  double chrg;  /* 1 for ion, -1 for electron */
  double zprt;
  double prot;  /* proton mass in proton units */
  /* particle distribution */
  double ekev;
  double polo_scale;
  double p1_scale;
  double p2_scale;
  double pchi;

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

Particles_t* Particles_ctor(){
  return (Particles_t*)calloc(1, sizeof(Particles_t));
}


void initialize_Particles(Particles_t* ptcl_ptr, Config_t* cfg_ptr){

  /* from the config we can get number of particles
     which we will use to dimension our arrays
  */
  ptcl_ptr->nprt = cfg_ptr->nprt;
  ptcl_ptr->idm = (unsigned)ptcl_ptr->nprt;

  ptcl_ptr->chrg = cfg_ptr->chrg;
  ptcl_ptr->zprt = cfg_ptr->zprt;
  ptcl_ptr->prot = cfg_ptr->prot;
  ptcl_ptr->ekev = cfg_ptr->ekev;
  ptcl_ptr->polo_scale = cfg_ptr->polo_scale;
  ptcl_ptr->p1_scale = cfg_ptr->p1_scale;
  ptcl_ptr->p2_scale = cfg_ptr->p2_scale;
  ptcl_ptr->pchi = cfg_ptr->pchi;


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

double get_zprt(Particles_t* ptcl_ptr){
  return ptcl_ptr->zprt;
}

double get_prot(Particles_t* ptcl_ptr){
  return ptcl_ptr->zprt;
}
