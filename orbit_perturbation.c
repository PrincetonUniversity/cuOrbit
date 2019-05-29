#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "orbit_perturbation.h"

const int NAMP_ = 155;

void initialize_Perturb(Perturb_t* ptrb_ptr){

  /* these values are to become part of config  */
  double falf = 13.3550;


  /* for now we'll load from file, same as the fortran code */
  FILE *ifp;
  const char *mode = "r";
  const char inputFilename[] = "INPUT/displ_3_9.dat";
  int j, ind, md;
  int idum;
  int nmd, mmin, mmax, ndum;
  double fkhz;
  double omrat;
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

  for(md=0; md < ptrb_ptr->modes; md++){
    fscanf(ifp, "%*s=%d", &(ptrb_ptr->mmod[md]));
    for(j=0; j < ptrb_ptr->lpt; j++){
      ind = ptrb_ptr->lpt * md + j;
      fscanf(ifp, "%lf ", &(ptrb_ptr->xi1[ind]));
    }  /* j */
  }  /* md */

  fclose(ifp);

  /* now from here we initialize remaining structure */

  

};

