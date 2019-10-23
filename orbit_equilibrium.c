#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "orbit_equilibrium.h"
#include "orbit_config_api.h"
#include "orbit_constants.h"
#include "orbit_util.h"
#include "cuda_helpers.h"

struct Equilib {
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
  double *B[9], *G[9], *R[9], *X[9], *Z[9];
  double *b1, *b2, *b3, *b4, *b5, *b6, *b7, *b8, *b9;
  double *g1, *g2, *g3, *g4, *g5, *g6, *g7, *g8, *g9;
  double *r1, *r2, *r3, *r4, *r5, *r6, *r7, *r8, *r9;
  double *x1, *x2, *x3, *x4, *x5, *x6, *x7, *x8, *x9;
  double *z1, *z2, *z3, *z4, *z5, *z6, *z7, *z8, *z9;
  double *GD[3], *PD[3], *PS[3], *QD[3], *RD[3], *RP[3];
  double *gd1, *gd2, *gd3;
  double *pd1, *pd2, *pd3;
  double *ps1, *ps2, *ps3;
  double *qd1, *qd2, *qd3;
  double *rd1, *rd2, *rd3;
  double *rp1, *rp2, *rp3;

  /* xxx maybe belongs somewhere else */
  double* VD[3];
  double* vd1;
  double* vd2;
  double* vd3;

};

Equilib_t* Equilib_ctor(){
  return (Equilib_t*)umacalloc(1, sizeof(Equilib_t));
}


void initialize_Equilib(Equilib_t* Equilib_ptr, Config_t* cfg_ptr){

  FILE *ifp;
  const char *mode = "r";
  char buf[255];
  int j, l, ind;
  size_t sz;

  ifp = fopen(cfg_ptr->spdata_file, mode);
  if (ifp == NULL) {
    fprintf(stderr, "Can't open input file %s!\n", cfg_ptr->spdata_file);
    exit(1);
  }
  printf("Parsing Equilibrium file %s\n",  cfg_ptr->spdata_file);

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

  Equilib_ptr->b1 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->B[0] = Equilib_ptr->b1;
  Equilib_ptr->b2 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->B[1] = Equilib_ptr->b2;
  Equilib_ptr->b3 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->B[2] = Equilib_ptr->b3;
  Equilib_ptr->b4 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->B[3] = Equilib_ptr->b4;
  Equilib_ptr->b5 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->B[4] = Equilib_ptr->b5;
  Equilib_ptr->b6 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->B[5] = Equilib_ptr->b6;
  Equilib_ptr->b7 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->B[6] = Equilib_ptr->b7;
  Equilib_ptr->b8 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->B[7] = Equilib_ptr->b8;
  Equilib_ptr->b9 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->B[8] = Equilib_ptr->b9;
  Equilib_ptr->g1 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->G[0] = Equilib_ptr->g1;
  Equilib_ptr->g2 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->G[1] = Equilib_ptr->g2;
  Equilib_ptr->g3 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->G[2] = Equilib_ptr->g3;
  Equilib_ptr->g4 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->G[3] = Equilib_ptr->g4;
  Equilib_ptr->g5 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->G[4] = Equilib_ptr->g5;
  Equilib_ptr->g6 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->G[5] = Equilib_ptr->g6;
  Equilib_ptr->g7 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->G[6] = Equilib_ptr->g7;
  Equilib_ptr->g8 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->G[7] = Equilib_ptr->g8;
  Equilib_ptr->g9 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->G[8] = Equilib_ptr->g9;
  Equilib_ptr->r1 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->R[0] = Equilib_ptr->r1;
  Equilib_ptr->r2 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->R[1] = Equilib_ptr->r2;
  Equilib_ptr->r3 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->R[2] = Equilib_ptr->r3;
  Equilib_ptr->r4 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->R[3] = Equilib_ptr->r4;
  Equilib_ptr->r5 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->R[4] = Equilib_ptr->r5;
  Equilib_ptr->r6 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->R[5] = Equilib_ptr->r6;
  Equilib_ptr->r7 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->R[6] = Equilib_ptr->r7;
  Equilib_ptr->r8 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->R[7] = Equilib_ptr->r8;
  Equilib_ptr->r9 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->R[8] = Equilib_ptr->r9;
  Equilib_ptr->x1 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->X[0] = Equilib_ptr->x1;
  Equilib_ptr->x2 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->X[1] = Equilib_ptr->x2;
  Equilib_ptr->x3 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->X[2] = Equilib_ptr->x3;
  Equilib_ptr->x4 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->X[3] = Equilib_ptr->x4;
  Equilib_ptr->x5 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->X[4] = Equilib_ptr->x5;
  Equilib_ptr->x6 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->X[5] = Equilib_ptr->x6;
  Equilib_ptr->x7 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->X[6] = Equilib_ptr->x7;
  Equilib_ptr->x8 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->X[7] = Equilib_ptr->x8;
  Equilib_ptr->x9 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->X[8] = Equilib_ptr->x9;
  Equilib_ptr->z1 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->Z[0] = Equilib_ptr->z1;
  Equilib_ptr->z2 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->Z[1] = Equilib_ptr->z2;
  Equilib_ptr->z3 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->Z[2] = Equilib_ptr->z3;
  Equilib_ptr->z4 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->Z[3] = Equilib_ptr->z4;
  Equilib_ptr->z5 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->Z[4] = Equilib_ptr->z5;
  Equilib_ptr->z6 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->Z[5] = Equilib_ptr->z6;
  Equilib_ptr->z7 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->Z[6] = Equilib_ptr->z7;
  Equilib_ptr->z8 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->Z[7] = Equilib_ptr->z8;
  Equilib_ptr->z9 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->Z[8] = Equilib_ptr->z9;

  sz = (unsigned)Equilib_ptr->lsp;
  Equilib_ptr->qd1 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->QD[0] = Equilib_ptr->qd1;
  Equilib_ptr->qd2 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->QD[1] = Equilib_ptr->qd2;
  Equilib_ptr->qd3 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->QD[2] = Equilib_ptr->qd3;
  Equilib_ptr->gd1 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->GD[0] = Equilib_ptr->gd1;
  Equilib_ptr->gd2 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->GD[1] = Equilib_ptr->gd2;
  Equilib_ptr->gd3 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->GD[2] = Equilib_ptr->gd3;
  Equilib_ptr->rd1 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->RD[0] = Equilib_ptr->rd1;
  Equilib_ptr->rd2 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->RD[1] = Equilib_ptr->rd2;
  Equilib_ptr->rd3 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->RD[2] = Equilib_ptr->rd3;
  Equilib_ptr->pd1 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->PD[0] = Equilib_ptr->pd1;
  Equilib_ptr->pd2 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->PD[1] = Equilib_ptr->pd2;
  Equilib_ptr->pd3 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->PD[2] = Equilib_ptr->pd3;
  Equilib_ptr->rp1 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->RP[0] = Equilib_ptr->rp1;
  Equilib_ptr->rp2 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->RP[1] = Equilib_ptr->rp2;
  Equilib_ptr->rp3 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->RP[2] = Equilib_ptr->rp3;
  Equilib_ptr->ps1 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->PS[0] = Equilib_ptr->ps1;
  Equilib_ptr->ps2 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->PS[1] = Equilib_ptr->ps2;
  Equilib_ptr->ps3 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->PS[2] = Equilib_ptr->ps3;

  sz = (size_t)IDP;
  Equilib_ptr->vd1 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->VD[0] = Equilib_ptr->vd1;
  Equilib_ptr->vd2 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->VD[1] = Equilib_ptr->vd2;
  Equilib_ptr->vd3 = (double*)umacalloc(sz, sizeof(double));
  Equilib_ptr->VD[2] = Equilib_ptr->vd3;
  Equilib_ptr->rp1 = (double*)umacalloc(sz, sizeof(double));

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

  fclose(ifp);

  /* /\* dbg G *\/ */
  /* /\* g is lsp *lst *\/ */
  /* double sumofag= 0.; */
  /* ind = 0; */
  /* for(j=0; j< Equilib_ptr->lsp; j++){ */
  /*   for(l=0; l< Equilib_ptr->lst; l++){ */
  /*     sumofag += Equilib_ptr->g1[ind]; */
  /*     ind++; */
  /*   } */
  /* } */
  /* printf("g1[0] %f sum(g1) = %f\n", Equilib_ptr->g1[0], sumofag); */

}

#ifdef __NVCC__
__host__ __device__
#endif
double** get_B(Equilib_t* Eq_ptr){
  return Eq_ptr->B;
}

#ifdef __NVCC__
__host__ __device__
#endif
double** get_G(Equilib_t* Eq_ptr){
  return Eq_ptr->G;
}

#ifdef __NVCC__
__host__ __device__
#endif
double** get_R(Equilib_t* Eq_ptr){
  return Eq_ptr->R;
}

#ifdef __NVCC__
__host__ __device__
#endif
double** get_X(Equilib_t* Eq_ptr){
  return Eq_ptr->X;
}

#ifdef __NVCC__
__host__ __device__
#endif
double** get_Z(Equilib_t* Eq_ptr){
  return Eq_ptr->Z;
}


#ifdef __NVCC__
__host__ __device__
#endif
double** get_QD(Equilib_t* Eq_ptr){
  return Eq_ptr->QD;
}

#ifdef __NVCC__
__host__ __device__
#endif
double** get_RD(Equilib_t* Eq_ptr){
  return Eq_ptr->RD;
}

#ifdef __NVCC__
__host__ __device__
#endif
double** get_GD(Equilib_t* Eq_ptr){
  return Eq_ptr->GD;
}

#ifdef __NVCC__
__host__ __device__
#endif
double** get_PD(Equilib_t* Eq_ptr){
  return Eq_ptr->PD;
}

#ifdef __NVCC__
__host__ __device__
#endif
double** get_PS(Equilib_t* Eq_ptr){
  return Eq_ptr->PS;
}

#ifdef __NVCC__
__host__ __device__
#endif
double** get_RP(Equilib_t* Eq_ptr){
  return Eq_ptr->RP;
}

#ifdef __NVCC__
__host__ __device__
#endif
double** get_VD(Equilib_t* Eq_ptr){
  return Eq_ptr->VD;
}



#ifdef __NVCC__
__host__ __device__
#endif
double get_pw(Equilib_t* Eq_ptr){
  return Eq_ptr->pw;
}

#ifdef __NVCC__
__host__ __device__
#endif
double get_rmaj(Equilib_t* Eq_ptr){
  return Eq_ptr->rmaj;
}

#ifdef __NVCC__
__host__ __device__
#endif
int compute_jd(Equilib_t* Eq_ptr, double x){
  /* computes jd in zero inds */
  int jd;
  jd = (int) (x * (Eq_ptr->lsp - 1) / Eq_ptr->pw);
  jd = imin(jd, Eq_ptr->lsp - 2);
  jd = imax(jd, 0);
  return jd;
}

#ifdef __NVCC__
__host__ __device__
#endif
static int compute_kd(Equilib_t* Eq_ptr, double x){
  /* computes kd in zero inds */
  const double lst = Eq_ptr->lst;
  int kd;
  kd = (int)(x * lst/pi2 );
  kd = imax(0, kd);
  kd = imin(lst-1, kd);
  return kd;
}

double gfun(Equilib_t* Eq_ptr, double px){
  /* gives g as function of poloidal flux */
  const int jd = compute_jd(Eq_ptr, px);
  const double dpx = px - jd * Eq_ptr->pw / (Eq_ptr->lsp-1);
  const double dp2 = dpx*dpx;
  return Eq_ptr->gd1[jd] + Eq_ptr->gd2[jd] * dpx + Eq_ptr->gd3[jd] * dp2;
}
double qfun(Equilib_t* Eq_ptr, double px){
  /* gives q as function of poloidal flux */
  const int jd = compute_jd(Eq_ptr, px);
  const double dpx = px - jd * Eq_ptr->pw / (Eq_ptr->lsp-1);
  const double dp2 = dpx*dpx;
  return Eq_ptr->qd1[jd] + Eq_ptr->qd2[jd] * dpx + Eq_ptr->qd3[jd] * dp2;
}

double rifun(Equilib_t* Eq_ptr, double px){
  /* gives ri as function of poloidal flux */
  const int jd = compute_jd(Eq_ptr, px);
  const double dpx = px - jd * Eq_ptr->pw / (Eq_ptr->lsp-1);
  const double dp2 = dpx*dpx;
  return Eq_ptr->rd1[jd] + Eq_ptr->rd2[jd] * dpx + Eq_ptr->rd3[jd] * dp2;
}

double giac(Equilib_t* Eq_ptr, double px, double tx){
  /* jacobian as function of poloidal flux */
  const double lst = Eq_ptr->lst;
  int idum;
  double tdum;
  const int jd = compute_jd(Eq_ptr, px);
  const double dpx = px - jd * Eq_ptr->pw / (Eq_ptr->lsp-1);
  const double dp2 = dpx*dpx;
  idum = tx/pi2;
  tdum = tx - pi2 * (double)(idum-1);
  idum = tdum *.1591549431;
  tdum = tdum - pi2 * (double)idum;
  const int kd = compute_kd(Eq_ptr, tdum);
  const double dtx = tdum - (double)(kd) * pi2 / lst;
  const double dt2 = dtx*dtx;
  const int ind = jd*lst + kd;
  //printf("DBG dpx %f dtx %f g1 %f\n",dpx,dtx, Eq_ptr->g1[ind]);
  return Eq_ptr->g1[ind] + Eq_ptr->g2[ind]*dpx + Eq_ptr->g3[ind]*dp2
      + Eq_ptr->g4[ind]*dtx + Eq_ptr->g5[ind]*dpx*dtx + Eq_ptr->g6[ind]*dtx*dp2
      + Eq_ptr->g7[ind]*dt2 + Eq_ptr->g8[ind]*dt2*dpx + Eq_ptr->g9[ind]*dt2*dp2;
}

#ifdef __NVCC__
__host__ __device__
#endif
double xproj(Equilib_t* Eq_ptr, double px, double tx){
  /* major radius position as function of poloidal flux, theta  */
  const double lst = Eq_ptr->lst;
  int idum;
  double tdum;
  const int jd = compute_jd(Eq_ptr, px);
  double dpx = px - jd * Eq_ptr->pw / (Eq_ptr->lsp-1);

  const double sdum = 0.5 - 0.5* copysign(1.,  jd - 0.5);
  dpx =  (1. - sdum) * dpx + sdum * sqrt(fmax(1.E-20, dpx));
  const double dp2 = dpx*dpx;

  idum = tx/pi2;
  tdum = tx - pi2 * (double)(idum-1);
  idum = tdum *.1591549431;
  tdum = tdum - pi2 * (double)idum;
  const int kd = compute_kd(Eq_ptr, tdum);
  const double dtx = tdum - (double)(kd) * pi2 / lst;
  const double dt2 = dtx*dtx;
  const int ind = jd*lst + kd;
  return Eq_ptr->x1[ind] + Eq_ptr->x2[ind]*dpx + Eq_ptr->x3[ind]*dp2
      + Eq_ptr->x4[ind]*dtx + Eq_ptr->x5[ind]*dpx*dtx + Eq_ptr->x6[ind]*dtx*dp2
      + Eq_ptr->x7[ind]*dt2 + Eq_ptr->x8[ind]*dt2*dpx + Eq_ptr->x9[ind]*dt2*dp2;
}

double zproj(Equilib_t* Eq_ptr, double px, double tx){
  /* vertical position as function of poloidal flux, theta  */
  const double lst = Eq_ptr->lst;
  int idum;
  double tdum;
  const int jd = compute_jd(Eq_ptr, px);
  double dpx = px - jd * Eq_ptr->pw / (Eq_ptr->lsp-1);

  const double sdum = 0.5 - 0.5* copysign(1.,  jd - 0.5);
  dpx =  (1. - sdum) * dpx + sdum * sqrt(fmax(1.E-20, dpx));
  const double dp2 = dpx*dpx;

  idum = tx/pi2;
  tdum = tx - pi2 * (double)(idum-1);
  idum = tdum *.1591549431;
  tdum = tdum - pi2 * (double)idum;
  const int kd = compute_kd(Eq_ptr, tdum);
  const double dtx = tdum - (double)(kd) * pi2 / lst;
  const double dt2 = dtx*dtx;
  const int ind = jd*lst + kd;
  return Eq_ptr->z1[ind] + Eq_ptr->z2[ind]*dpx + Eq_ptr->z3[ind]*dp2
      + Eq_ptr->z4[ind]*dtx + Eq_ptr->z5[ind]*dpx*dtx + Eq_ptr->z6[ind]*dtx*dp2
      + Eq_ptr->z7[ind]*dt2 + Eq_ptr->z8[ind]*dt2*dpx + Eq_ptr->z9[ind]*dt2*dp2;
}


double bfield(Equilib_t* Eq_ptr, double px, double tx){
  double dpx;
  double tdum;
  int idum;
  int kd;

  const int lst = Eq_ptr->lst;
  const int jd = compute_jd(Eq_ptr, px);

  dpx = px - ((double)jd) * Eq_ptr->pw / (Eq_ptr->lsp-1);
  const double sdum = .5 - .5 * copysign(1., jd - 0.5);
  dpx = (1. - sdum) * dpx + sdum * sqrt( fmax(1E-20, dpx));
  const double dp2 = dpx*dpx;
  idum = (int)(tx * pi2i);
  tdum = (int)(tx - pi2 * (idum - 1));
  idum = (int)(tdum * pi2i);
  tdum = tdum - pi2*idum;
  kd = (int)(tdum * lst * pi2i );
  kd = imax(0, kd);
  kd = imax(lst - 1, kd);
  const double dtx = tdum - (kd-1) * pi2 / lst;
  const double dt2 = dtx*dtx;
  const int ind = jd * lst + kd;
  return Eq_ptr->b1[ind] + Eq_ptr->b2[ind]*dpx + Eq_ptr->b3[ind]*dp2
      + Eq_ptr->b4[ind]*dtx + Eq_ptr->b5[ind]*dpx*dtx + Eq_ptr->b6[ind]*dtx*dp2
      + Eq_ptr->b7[ind]*dt2 + Eq_ptr->b8[ind]*dt2*dpx + Eq_ptr->b9[ind]*dt2*dp2;
}

#ifdef __NVCC__
__host__ __device__
#endif
double get_ped(Equilib_t* Eq_ptr){
  return Eq_ptr->ped;
}

#ifdef __NVCC__
__host__ __device__
#endif
double get_lsp(Equilib_t* Eq_ptr){
  return Eq_ptr->lsp;
}

#ifdef __NVCC__
__host__ __device__
#endif
double get_lst(Equilib_t* Eq_ptr){
  return Eq_ptr->lst;
}

#ifdef __NVCC__
__host__ __device__
#endif
double get_nrip(Equilib_t* Eq_ptr){
  return Eq_ptr->nrip;
}

#ifdef __NVCC__
__host__ __device__
#endif
double get_krip(Equilib_t* Eq_ptr){
  return Eq_ptr->krip;
}


void vspline(Equilib_t* Eq_ptr){
  const int lsp = get_lsp(Eq_ptr);
  const int lspm = lsp - 1;
  double* const vd1 = get_VD(Eq_ptr)[0];
  double* const vd2 = get_VD(Eq_ptr)[1];
  double* const vd3 = get_VD(Eq_ptr)[1];
  double dpx, integ, pdum, tdum;
  int j, k, l;
  int jm, jp, jpp;

  dpx = 0.01 * get_pw(Eq_ptr) / lspm;
  /* vd1 */
  vd1[0] = 0.;
  integ = 0;
  pdum = 0.;
  for(j=1; j<lsp; j++){
    for(k=0; k<100; k++){
      pdum += dpx;
      for(l=1; l<=100; l++){
        tdum = .01 * l * pi2;
        integ = integ + .01 * pi2 * dpx * giac(Eq_ptr, pdum, tdum);
      }  /* l */
    }    /* k */
    vd1[j] = integ;
  }  /* j */
  for(j=0; j<lsp; j++){
    vd1[j] /= vd1[lsp];
  }

  dpx = get_pw(Eq_ptr) / lspm;
  /*  vd2 */
  vd2[0] = (10*vd1[1] - 7*vd1[0] - 3*vd1[2])/(4*dpx);
  for(j=1; j<lspm; j++){
    jm = j - 1;
    jp = j + 1;
    jpp = imin(j + 2, lspm);
    vd2[j] = -vd2[jm] + 2*(vd1[j]-vd1[jm])/dpx;
      /* smooth a1 */
    vd1[jp] = .4 * dpx * vd2[j] + .3 * vd1[jpp] + .7*vd1[j];
  }
  for(j=0; j<lspm; j++){
    jp = j+1;
    vd3[j] = (vd2[jp] - vd2[j])/(2.*dpx);
  }
  return;
}

#ifdef __NVCC__
__host__ __device__
#endif
double rpol(Equilib_t* Eq_ptr, double pdum){
  int jd;
  const int lspm1 = Eq_ptr->lsp - 1;
  const double pw = Eq_ptr->pw;
  jd = ((int) pdum * lspm1 / pw) + 1;
  jd = imin(jd, lspm1);
  jd = imax(jd, 1);
  const double dpx = pdum - (jd -1) * pw / (lspm1);
  const double dp2 = dpx * dpx;
  /* zero inds */
  jd--;
  return Eq_ptr->rp1[jd] + Eq_ptr->rp2[jd]*dpx + Eq_ptr->rp3[jd]*dp2;
}

#ifdef __NVCC__
__host__ __device__
#endif
double polr_mp(Equilib_t* Eq_ptr, double rdum, double pd){
  const double eps = compute_eps(Eq_ptr);
  if(rdum >= eps){
    return Eq_ptr->pw;
  }
  if(rdum <= 0.){
    return 1.E-8;
  }

  int nstp = 0;
  const double rtol = 1.E-3;

  double rz = rpol(Eq_ptr, pd);

  while( fabs(rz-rdum) > rtol && nstp<1){
    pd *= (rdum/rz);
    rz = rpol(Eq_ptr, pd);
  }

  return pd;
}

#ifdef __NVCC__
__host__ __device__
#endif
double compute_eps(Equilib_t* eqlb_ptr){
  const double xproj0 = xproj(eqlb_ptr, 0., 0.);
  const double xproj1 = xproj(eqlb_ptr, eqlb_ptr->ped, 0.);

  return  (xproj1 - xproj0)/xproj0;
}
