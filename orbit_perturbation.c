#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

/* for constants like pi */
#include <math.h>

#include "orbit_equilibrium.h"
#include "orbit_perturbation.h"
#include "orbit_particles.h"
#include "orbit_util.h"

const int NAMP_ = 155;

void initialize_Perturb(Perturb_t* ptrb_ptr, Config_t* config_ptr,
                        Equilib_t* equilib_ptr, Particles_t* ptcl_ptr){

  /* these values are to become part of config  */
  double falf = 13.3550;
  double ascale = 5.;
  double global_scaling_factor = 3E-4;
  double alimit = 0.;
  double freq_scaling_factor = 1.;

  /* for now we'll load from file, same as the fortran code */
  FILE *ifp;
  const char *mode = "r";
  const char inputFilename[] = "INPUT/displ_3_9.dat";
  int j, k, m, n, ind, md;
  int idum;
  int nmd, mmin, mmax, ndum;
  double fkhz;
  double omrat;
  double xxmax;
  double px;
  int lptm1;
  size_t sz;


  ifp = fopen(inputFilename, mode);
  if (ifp == NULL) {
    fprintf(stderr, "Can't open input file %s!\n", inputFilename);
    exit(1);
  }
  printf("Parsing Displacement file %s\n",  inputFilename);

  /* file header lines */
  fscanf(ifp, "%*[^\n]\n");  /* skip 0 */
  fscanf(ifp, "%*[^\n]\n");  /* skip 1 */
  fscanf(ifp, "%*[^\n]\n");  /* skip 2 */

  fscanf(ifp, "%d %d %d %d %lf %d ", &(ptrb_ptr->lpt), &nmd, &mmin, &mmax, &omrat, &ndum);
  fscanf(ifp, "%*[^\n]\n");  /*  skip remaing part of line */
  fkhz = omrat * falf;
  lptm1 = ptrb_ptr->lpt - 1;
  printf("lpt = %d\nnmd = %d\nmmin = %d\nmmax = %d\nomrat = %g\nfkhz = %g\n",
         ptrb_ptr->lpt, nmd, mmin, mmax, omrat, fkhz);

  fscanf(ifp, "%*[^\n]\n");  /* skip 4 */
  fscanf(ifp, "%d ", &idum);
  assert (idum == ptrb_ptr->lpt);
  for(j=0; j < ptrb_ptr->lpt; j++){
    /* zero md */ //xxx according to code, we can just skip these...overwrite below... confirm with Mario
    //fscanf(ifp, "%lf ", &(ptrb_ptr->xi1[j]));
    fscanf(ifp, "%*lf ");
  }

  fscanf(ifp, "%*[^\n]\n");  /* skip*/
  fscanf(ifp, "%d %d ", &idum, &(ptrb_ptr->modes));
  assert(idum == ptrb_ptr->lpt);

  /* malloc xi */
  printf("Malloc'ing arrays for Perturbation\n");
  sz = (unsigned)(ptrb_ptr->lpt * ptrb_ptr->modes);
  ptrb_ptr->xi1 = (double*)calloc(sz, sizeof(double));
  ptrb_ptr->xi2 = (double*)calloc(sz, sizeof(double));
  ptrb_ptr->xi3 = (double*)calloc(sz, sizeof(double));
  ptrb_ptr->a1 = (double*)calloc(sz, sizeof(double));
  ptrb_ptr->a2 = (double*)calloc(sz, sizeof(double));
  ptrb_ptr->a3 = (double*)calloc(sz, sizeof(double));

  ptrb_ptr->xx = (double*)calloc((unsigned)ptrb_ptr->lpt, sizeof(double));

  ptrb_ptr->alfv = (double*)calloc((unsigned)ptrb_ptr->modes, sizeof(double));
  ptrb_ptr->amp = (double*)calloc((unsigned)ptrb_ptr->modes, sizeof(double));
  ptrb_ptr->damp = (double*)calloc((unsigned)ptrb_ptr->modes, sizeof(double));
  ptrb_ptr->harm = (double*)calloc((unsigned)ptrb_ptr->modes, sizeof(double));
  ptrb_ptr->omegv = (double*)calloc((unsigned)ptrb_ptr->modes, sizeof(double));

  ptrb_ptr->mmod = (int*)calloc((unsigned)ptrb_ptr->modes, sizeof(int));
  ptrb_ptr->nmod = (int*)calloc((unsigned)ptrb_ptr->modes, sizeof(int));



  for(md=0; md < ptrb_ptr->modes; md++){
    fscanf(ifp, "** m = %d ", &(ptrb_ptr->mmod[md]));
    assert(md == ptrb_ptr->mmod[md]);
    for(j=0; j < ptrb_ptr->lpt; j++){
      ind = ptrb_ptr->lpt * md + j;
      fscanf(ifp, "%lf ", &(ptrb_ptr->xi1[ind]));
    }  /* j */
  }  /* md */

  fclose(ifp);

  /* now from here we initialize remaining structure */

  /* who do these work for */
  //nval = 1
  //nval1 = 1

  for(md=0; md < ptrb_ptr->modes; md++){
    ptrb_ptr->harm[md] = 1;
    ptrb_ptr->alfv[md] = md;
    /* this is tricky, including config, do we need to do it this way... */
    ptrb_ptr->omegv[md] = 2E3 * M_PI * fkhz / config_ptr->omeg0;
    ptrb_ptr->nmod[md] = nmd;
  }

  // who does this belong to
  //xxxnvalx = nval;

  printf("Change Amplitude = %g\nChange Frequency = %g\nLimit = %g\n",
         global_scaling_factor, freq_scaling_factor, alimit);

  md = -1;
  for(k=0; k < ptrb_ptr->modes; k++){
    for(j=0; j < ptrb_ptr->lpt; j++){
      ind = ptrb_ptr->lpt * k + j;
      ptrb_ptr->xx[j] = fabs(ptrb_ptr->xi1[ind]);
    }
    xxmax = darray_max(ptrb_ptr->xx, (unsigned)ptrb_ptr->lpt);
    ptrb_ptr->damp[k] = xxmax;

    if (ptrb_ptr->damp[k] < alimit) continue;
    if (ptrb_ptr->mmod[k] == 0) continue;

    md++;
    ptrb_ptr->amp[md] = ascale * xxmax * global_scaling_factor; /* modify  scaling */
    ptrb_ptr->omegv[md] = ptrb_ptr->omegv[k] * freq_scaling_factor; /* modify freq */
    ptrb_ptr->mmod[md] = ptrb_ptr->mmod[k];
    ptrb_ptr->nmod[md] = ptrb_ptr->nmod[k];
    m = ptrb_ptr->mmod[md];
    n = ptrb_ptr->nmod[md];

    for(j=0; j < ptrb_ptr->lpt; j++){
      ind = ptrb_ptr->lpt * md + j;
      /*  normalize, careful inds */
      ptrb_ptr->xi1[ind] = ptrb_ptr->xi1[ ptrb_ptr->lpt * k + j] /xxmax;

      /* compute alpha from disp */
      px = j * equilib_ptr->pw / lptm1;
      ptrb_ptr->a1[ind] = ( m / qfun(equilib_ptr, px) - n) * ptrb_ptr->xi1[ind] / (
          m * gfun(equilib_ptr, px) + n * rifun(equilib_ptr, px));
    }  /* j */
    printf("%d %d %d %f %f %f\n",
           md, ptrb_ptr->nmod[md], ptrb_ptr->mmod[md], 1E5*ptrb_ptr->amp[md],
           ptrb_ptr->omegv[md] * (config_ptr->omeg0) / 6280.,
           xxmax);
  }  /* k */

  /*  does this really need 1 based?, lets try without*/
  //ptrb_ptr->md1 = 0;
  /* and do we need md2 at all? looks vestigile, lets try without*/
  //ptrb_ptr->md2 = ptrb_ptr->modes;

  splna(ptrb_ptr, equilib_ptr, ptcl_ptr);
  splnx(ptrb_ptr, equilib_ptr, ptcl_ptr);
  return;
};


void splna(Perturb_t* ptrb_ptr, Equilib_t* equilib_ptr, Particles_t* ptcl_ptr){
  int ind, j, m, md;
  int jm, jp, jpp;
  const int lpt = ptrb_ptr->lpt;
  const int lptm = lpt - 1;
  const double dpx = equilib_ptr->pw / (double)lptm;

  const int lpx = 1;  /* mp change mar 2016 */

  for(md=0; md < ptrb_ptr->modes; md++){
    m = ptrb_ptr->mmod[md];
    ind = md * lpt;
    for(j=0; j < lpx; j++){
      ptrb_ptr->xi1[ind + j] = pow(ptcl_ptr->pol[j], m) *
          ptrb_ptr->xi1[md*lpt + (lpx-1)] /
          pow(ptcl_ptr->pol[lpx-1], m);
    }
    ptrb_ptr->a2[md*lpt] = (
        10. * ptrb_ptr->a1[md*lpt + 1] -
        7.  * ptrb_ptr->a1[md*lpt] -
        3.  * ptrb_ptr->a1[md*lpt + 2]
                            ) / (4. * dpx);
    if (m != 1){
      ptrb_ptr->a2[md*lpt] = 0;
    }
    for(j=1; j<lptm; j++){
      jm = j-1;
      jp = j+1;
      jpp = imin(j+2, lpt);

      ptrb_ptr->a2[ind + j] =
          -1. * ptrb_ptr->a2[ind + jm] +
          2. * (ptrb_ptr->a1[ind + j] - ptrb_ptr->a1[ind + jm]) / dpx;
      /* smooth a1 */

      ptrb_ptr->a1[ind + jp] =
          0.4 * dpx * ptrb_ptr->a2[ind + j] +
          0.3 * ptrb_ptr->a1[ind + jpp] +
          0.7 * ptrb_ptr->a1[ind + j];
    }  /* j */

    for(j=0; j<lptm; j++){
      ptrb_ptr->a3[ind + j] = (
          ptrb_ptr->a2[ind + j + 1] - ptrb_ptr->a2[ind + j] ) / (2. * dpx);
    }  /* j */
  }    /* md */
  return;
}

void splnx(Perturb_t* ptrb_ptr, Equilib_t* equilib_ptr, Particles_t* ptcl_ptr){
  int ind, j, m, md;
  int jm, jp, jpp;
  const int lpt = ptrb_ptr->lpt;
  const int lptm = lpt - 1;
  const double dpx = equilib_ptr->pw / (double)lptm;

  const int lpx = 1;  /* mp change mar 2016 */

  for(md=0; md<ptrb_ptr->modes; md++){
    m = ptrb_ptr->mmod[md];
    ind = md * lpt;
    for(j=0; j<lpx; j++){
      ptrb_ptr->xi1[ind + j] = pow(ptcl_ptr->pol[j], m) * ptrb_ptr->xi1[ind + (lpx-1) ] /
          pow(ptcl_ptr->pol[lpx-1], m);
    }
    ptrb_ptr->xi2[ind] =  (
        10. * ptrb_ptr->xi1[ind + 1] -
        7.  * ptrb_ptr->xi1[ind] -
        3.  * ptrb_ptr->xi1[ind + 2]
                           ) / ( 4. * dpx);
    if(m != 1){
      ptrb_ptr->xi2[ind] = 0;
    }

    for(j=1; j<lptm; j++){
      jm = j-1;
      jp = j+1;
      jpp = imin(j+2, lpt);

      ptrb_ptr->xi2[ind + j] =
          -1. * ptrb_ptr->xi2[ind + jm] +
          2. * (ptrb_ptr->xi1[ind + j] - ptrb_ptr->xi1[ind + jm]) / dpx;
      /* smooth a1 */

      ptrb_ptr->xi1[ind + jp] =
          0.4 * dpx * ptrb_ptr->xi2[ind + j] +
          0.3 * ptrb_ptr->xi1[ind + jpp] +
          0.7 * ptrb_ptr->xi1[ind + j];
    }  /* j */

    for(j=0; j<lptm; j++){
      ptrb_ptr->xi3[ind + j] = (
          ptrb_ptr->xi2[ind + j + 1] - ptrb_ptr->xi2[ind + j] ) / (2. * dpx);
    }  /* j */
  }    /* md */
  return;
}

