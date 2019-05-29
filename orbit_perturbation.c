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
  char buf[255];
  int j, l, ind;
  size_t sz;
  int nmd, mmin, mmax, ndum;
  double fkhz;
  double omrat;
  int lptm1;


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
  fkhz = omrat * falf;
  lptm1 = ptrb_ptr->lpt - 1;
  printf("ltp = %d\nnmd = %d\nmmin = %d\nmmax = %d\nomrat = %g\nfkhz = %g\n",
         ptrb_ptr->lpt, nmd, mmin, mmax, omrat, fkhz);


};

