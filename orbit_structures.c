#include <stdlib.h>
#include <stdio.h>
#include "orbit_structures.h"


void initialize_config(Config_t* cfg_ptr){

  //Config_t cfg = *cfg_ptr;

  cfg_ptr->seed = 12345;
  cfg_ptr->ekev = 75;
  /* engn = ; */
  cfg_ptr->bkg = 5.5;
  /* double bmin; */
  /* double bmax; */
  /* double rmaj; */
  /* double trun; */
  /* double tran;  /\* transit time for particle at the mag axi with pitch=1 *\/ */
  /* double dt0; */
  /* double omeg0; */
  /* double xc; */
  /* double eps; */
  /* double bax; */
  /* double* dwal;  /\* wall definition *\/ */

  cfg_ptr->pamp = 8.63;
  cfg_ptr->rprof = 3.765;

};

void initialize_particle(Particle_t* ptc_ptr){

};

void initialize_Perturb(Perturb_t* ptrb_ptr){

};

void initialize_Equilib(Equilib_t* Equilib_ptr){


  /* for now we'll load from file, same as the fortran code */
  FILE *ifp;
  const char *mode = "r";
  const char inputFilename[] = "INPUT/spdata";
  char buf[255];
  int j, l, ind;
  size_t sz;

  ifp = fopen(inputFilename, mode);
  if (ifp == NULL) {
    fprintf(stderr, "Can't open input file in.list!\n");
    exit(1);
  }

  /* file header line */
  fscanf(ifp, "%s", buf);
  printf("CDF File: %s\n", buf);
  fscanf(ifp, "%s ", buf);
  printf("CDF Version: %s\n", buf);

  /* equilib array structure line */
   fscanf(ifp, "%d %d %d %d ",
         &(Equilib_ptr->lsp),
         &(Equilib_ptr->lst),
         &(Equilib_ptr->lemax),
         &(Equilib_ptr->lrmax));
  printf("lsp: %d\nlst %d\nlemax %d\nlrmax %d\n",
         Equilib_ptr->lsp, Equilib_ptr->lst, Equilib_ptr->lemax, Equilib_ptr->lrmax);

  printf("Malloc'ing Equilib arrays\n");
  sz = (unsigned)(Equilib_ptr->lsp * Equilib_ptr->lst);
  Equilib_ptr->b1 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->b2 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->b3 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->b4 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->b5 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->b6 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->b7 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->b8 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->b9 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->g1 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->g2 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->g3 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->g4 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->g5 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->g6 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->g7 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->g8 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->g9 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->r1 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->r2 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->r3 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->r4 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->r5 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->r6 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->r7 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->r8 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->r9 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->x1 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->x2 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->x3 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->x4 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->x5 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->x6 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->x7 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->x8 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->x9 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->z1 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->z2 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->z3 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->z4 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->z5 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->z6 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->z7 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->z8 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->z9 = (double*)calloc(sz, sizeof(double));
  sz = (unsigned)Equilib_ptr->lsp;
  Equilib_ptr->qd1 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->qd2 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->qd3 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->gd1 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->gd2 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->gd3 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->rd1 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->rd2 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->rd3 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->pd1 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->pd2 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->pd3 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->rp1 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->rp2 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->rp3 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->ps1 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->ps2 = (double*)calloc(sz, sizeof(double));
  Equilib_ptr->ps3 = (double*)calloc(sz, sizeof(double));


  /* first data line */
  fscanf(ifp, "%lf %lf ",  /* read doubles as "long float", ugh */
         &(Equilib_ptr->pw),
         &(Equilib_ptr->ped));
  printf("pw: %f\nped %f\n",
         Equilib_ptr->pw, Equilib_ptr->ped);

  /* /\*  we unfortunately need to surmise values from the middle of the file *\/ */
  /* seek_pos = Equilib_ptr->lsp * Equilib_ptr->lst * 9 * 4 + */
  /*     Equilib_ptr->lsp * 3 * 6; */

  /* consume the B G R X Z and derivatives */
  for(j=0; j<Equilib_ptr->lsp; j++){
    ind = Equilib_ptr->lst*j;
    /* B */
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->b1[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->b2[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->b3[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->b4[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->b5[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->b6[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->b7[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->b8[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->b9[ind + l]));
    }
    /* G */
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->g1[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->g2[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->g3[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->g4[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->g5[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->g6[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->g7[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->g8[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->g9[ind + l]));
    }
    /* X */
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->x1[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->x2[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->x3[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->x4[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->x5[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->x6[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->x7[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->x8[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->x9[ind + l]));
    }
    /* Z */
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->z1[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->z2[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->z3[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->z4[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->z5[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->z6[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->z7[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->z8[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->z9[ind + l]));
    }

    /*  */
    fscanf(ifp, "%lf %lf %lf ",
           &(Equilib_ptr->qd1[j]),
           &(Equilib_ptr->qd2[j]),
           &(Equilib_ptr->qd3[j]));
    fscanf(ifp, "%lf %lf %lf ",
           &(Equilib_ptr->gd1[j]),
           &(Equilib_ptr->gd2[j]),
           &(Equilib_ptr->gd3[j]));
    fscanf(ifp, "%lf %lf %lf ",
           &(Equilib_ptr->rd1[j]),
           &(Equilib_ptr->rd2[j]),
           &(Equilib_ptr->rd3[j]));
    fscanf(ifp, "%lf %lf %lf ",
           &(Equilib_ptr->pd1[j]),
           &(Equilib_ptr->pd2[j]),
           &(Equilib_ptr->pd3[j]));
    fscanf(ifp, "%lf %lf %lf ",
           &(Equilib_ptr->rp1[j]),
           &(Equilib_ptr->rp2[j]),
           &(Equilib_ptr->rp3[j]));
    fscanf(ifp, "%lf %lf %lf ",
           &(Equilib_ptr->ps1[j]),
           &(Equilib_ptr->ps2[j]),
           &(Equilib_ptr->ps3[j]));
  }  /* j */
  fscanf(ifp, "%d %d ", &(Equilib_ptr->krip),  &(Equilib_ptr->nrip));
  fscanf(ifp, "%lf %lf %lf ", &(Equilib_ptr->rmaj),  &(Equilib_ptr->d0), &(Equilib_ptr->brip));
  fscanf(ifp, "%lf %lf ", &(Equilib_ptr->wrip),  &(Equilib_ptr->xrip));
  for(j=0; j<Equilib_ptr->lsp; j++){
    ind = Equilib_ptr->lst*j;
    /* R */
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->r1[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->r2[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->r3[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->r4[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->r5[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->r6[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->r7[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->r8[ind + l]));
    }
    for(l=0; l<Equilib_ptr->lst; l++){
      fscanf(ifp, "%lf ", &(Equilib_ptr->r9[ind + l]));
    }
  }  /* j */

};


// test
int main(){
  Equilib_t Eq;
  initialize_Equilib(&Eq);
}
