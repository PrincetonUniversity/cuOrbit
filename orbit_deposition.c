#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "orbit_config_api.h"
#include "orbit_deposition.h"
#include "orbit_particles.h"  /* ekev */

const size_t MAXFNAME=255;

typedef struct Deposition {
  /* pdedp */
  char* pdedp_file;
  int nruns;
  bool compute_pdedp;
  bool initial_update_pdedp;
  double deposit_on_bins_after_fraction;
  double pdedp_dtsamp;
  double pdedp_dtav;
  int pdedp_tskip;
  double pdedp_otpup;
  bool pdedp_focusdep;
  bool pdedp_optimize;

  /* the 5d dist */
  double* pde_pdedp;

  /* internal vars */
  int pde_nbinDE;
  int pde_nbinDPz;
  int pde_nbinE;
  int pde_nbinPz;
  int pde_nbinmu;
  double* pde_varDE;
  double* pde_varDPz;
  double* pde_varE;
  double* pde_varPz;
  double* pde_varmu;
  double pde_Emin;
  double pde_Emax;
  double pde_Pzmin;
  double pde_Pzmax;
  double pde_mumin;
  double pde_mumax;
  double pde_DEmin;
  double pde_DEmax;
  double pde_DPzmin;
  double pde_DPzmax;


  /* xxx stochastic, does this belong here? */
  double mubk;
  int emink;
  int emaxk;
  int dmubk;
  int nstoche;
  int nstochp;

} Deposition_t;

Deposition_t* Deposition_ctor(){
  return (Deposition_t*)calloc(1, sizeof(Deposition_t));
}

void initialize_Deposition(Deposition_t* Depo_ptr, Config_t* cfg_ptr){

  Depo_ptr->pdedp_file = strndup(cfg_ptr->pdedp_file, MAXFNAME);
  Depo_ptr->nruns = cfg_ptr->nruns;
  Depo_ptr->compute_pdedp = cfg_ptr->compute_pdedp;
  Depo_ptr->initial_update_pdedp = cfg_ptr->initial_update_pdedp;
  Depo_ptr->deposit_on_bins_after_fraction = cfg_ptr->deposit_on_bins_after_fraction;
  Depo_ptr->pdedp_dtsamp = cfg_ptr->pdedp_dtsamp;
  Depo_ptr->pdedp_dtav = cfg_ptr->pdedp_dtav;
  Depo_ptr->pdedp_tskip = cfg_ptr->pdedp_tskip;
  Depo_ptr->pdedp_otpup = cfg_ptr->pdedp_otpup;
  Depo_ptr->pdedp_focusdep = cfg_ptr->pdedp_focusdep;
  Depo_ptr->pdedp_optimize = cfg_ptr->pdedp_optimize;

  /* xxx stochastic,  does this actually belong here*/
  Depo_ptr->mubk = cfg_ptr->mubk_scale *  get_ekev(cfg_ptr->ptcl_ptr);
  Depo_ptr->emink = cfg_ptr->emink;
  Depo_ptr->emaxk = cfg_ptr->emaxk;
  Depo_ptr->dmubk = cfg_ptr->dmubk;
  Depo_ptr->nstoche = cfg_ptr->nstoche;
  Depo_ptr->nstochp = cfg_ptr->nstochp;

  return;
}

bool compute_pdedp(Deposition_t* Depo_ptr){
  return Depo_ptr->compute_pdedp;
}

bool initial_update_pdedp(Deposition_t* Depo_ptr){
  return Depo_ptr->initial_update_pdedp;
}

bool pdedp_optimize(Deposition_t* Depo_ptr){
  return Depo_ptr->pdedp_optimize;
}


int get_pdedp_dtsamp(Deposition_t* Depo_ptr){
  return Depo_ptr->pdedp_dtsamp;
}

int get_pdedp_tskip(Deposition_t* Depo_ptr){
  return Depo_ptr->pdedp_tskip;
}
void set_pdedp_tskip(Deposition_t* Depo_ptr, double pdedp_tskip){
  Depo_ptr->pdedp_tskip = pdedp_tskip;
  return;
}



void pdedp_read(Deposition_t* Depo_ptr){
  /* this reads the probability distribution data
    for P(DE,DP| E,P,mu) from a text file (UFILE)*/
  int i;
  int je, jp, jmu, jde, jdp;
  size_t sz;
  FILE *ifp;
  const char *mode = "r";

  /* xxx do we need this, skip for now */
  class_domain(Depo_ptr);

  ifp = fopen(Depo_ptr->pdedp_file, mode);
  if (ifp == NULL) {
    fprintf(stderr, "\nCan't open input file %s!\n", Depo_ptr->pdedp_file);
    exit(1);
  }
  printf("\nParsing P(DE,DP) file %s\n",  Depo_ptr->pdedp_file);

  /* parse the file */
  for (i=0; i<9; i++){
    fscanf(ifp, "%*[^\n]\n");  /* skip line */
  }

  /* read in the time step */
  fscanf(ifp, "%lf %*[^\n]\n", &Depo_ptr->pdedp_dtsamp);
  Depo_ptr->pdedp_dtav = Depo_ptr->pdedp_dtsamp/50.;

  /* continue skipping header */
  fscanf(ifp, "%*[^\n]\n");  /* skip line */
  fscanf(ifp, "%*[^\n]\n");  /* skip line */

  fscanf(ifp, "%d %*[^\n]\n", &Depo_ptr->pde_nbinDE);
  fscanf(ifp, "%d %*[^\n]\n", &Depo_ptr->pde_nbinDPz);
  fscanf(ifp, "%d %*[^\n]\n", &Depo_ptr->pde_nbinE);
  fscanf(ifp, "%d %*[^\n]\n", &Depo_ptr->pde_nbinPz);
  fscanf(ifp, "%d %*[^\n]\n", &Depo_ptr->pde_nbinmu);

  printf("\nNumber of points for p(DE,DP|E,P,mu) function:\n");
  printf("\t E\t P\t  mu\t DE\t DP\n");
  printf("\t%d\t%d\t%d\t%d\t%d\n",
        Depo_ptr->pde_nbinE, Depo_ptr->pde_nbinPz, Depo_ptr->pde_nbinmu,
        Depo_ptr->pde_nbinDE, Depo_ptr->pde_nbinDPz);

  Depo_ptr->pde_varDE = (double*)calloc((unsigned)Depo_ptr->pde_nbinDE, sizeof(double));
  for(i=0; i<Depo_ptr->pde_nbinDE; i++){
    fscanf(ifp, "%lf ", &Depo_ptr->pde_varDE[i]);
  }
  Depo_ptr->pde_varDPz = (double*)calloc((unsigned)Depo_ptr->pde_nbinDPz, sizeof(double));
  for(i=0; i<Depo_ptr->pde_nbinDPz; i++){
    fscanf(ifp, "%lf ", &Depo_ptr->pde_varDPz[i]);
  }
  Depo_ptr->pde_varE = (double*)calloc((unsigned)Depo_ptr->pde_nbinE, sizeof(double));
  for(i=0; i<Depo_ptr->pde_nbinE; i++){
    fscanf(ifp, "%lf ", &Depo_ptr->pde_varE[i]);
  }
  Depo_ptr->pde_varPz = (double*)calloc((unsigned)Depo_ptr->pde_nbinPz, sizeof(double));
  for(i=0; i<Depo_ptr->pde_nbinPz; i++){
    fscanf(ifp, "%lf ", &Depo_ptr->pde_varPz[i]);
  }
  Depo_ptr->pde_varmu = (double*)calloc((unsigned)Depo_ptr->pde_nbinmu, sizeof(double));
  for(i=0; i<Depo_ptr->pde_nbinmu; i++){
    fscanf(ifp, "%lf ", &Depo_ptr->pde_varmu[i]);
  }

  //     Define boundaries of computing grid.
  Depo_ptr->pde_Emin = Depo_ptr->pde_varE[0] - .5 * (
      Depo_ptr->pde_varE[Depo_ptr->pde_nbinE-1] - Depo_ptr->pde_varE[0]) /(
          Depo_ptr->pde_nbinE - 1.);

  Depo_ptr->pde_Emax = Depo_ptr->pde_varE[Depo_ptr->pde_nbinE - 1] + .5*(
      Depo_ptr->pde_varE[Depo_ptr->pde_nbinE] - Depo_ptr->pde_varE[0])/(
          Depo_ptr->pde_nbinE - 1.);

  Depo_ptr->pde_Pzmin = Depo_ptr->pde_varPz[0] - .5*(
      Depo_ptr->pde_varPz[Depo_ptr->pde_nbinPz-1] - Depo_ptr->pde_varPz[0])/(
          Depo_ptr->pde_nbinPz - 1.);

  Depo_ptr->pde_Pzmax = Depo_ptr->pde_varPz[Depo_ptr->pde_nbinPz-1] + .5*(
      Depo_ptr->pde_varPz[Depo_ptr->pde_nbinPz] - Depo_ptr->pde_varPz[0])/(
          Depo_ptr->pde_nbinPz - 1.);

  Depo_ptr->pde_mumin = Depo_ptr->pde_varmu[0] - .5*(
      Depo_ptr->pde_varmu[Depo_ptr->pde_nbinmu - 1] - Depo_ptr->pde_varmu[0])/(
          Depo_ptr->pde_nbinmu - 1.);

  Depo_ptr->pde_mumax = Depo_ptr->pde_varmu[Depo_ptr->pde_nbinmu-1] + .5*(
      Depo_ptr->pde_varmu[Depo_ptr->pde_nbinmu - 1] - Depo_ptr->pde_varmu[0])/(
          Depo_ptr->pde_nbinmu - 1.);

  Depo_ptr->pde_DEmin = Depo_ptr->pde_varDE[0] - .5*(
      Depo_ptr->pde_varDE[Depo_ptr->pde_nbinDE-1] - Depo_ptr->pde_varDE[0])/(
          Depo_ptr->pde_nbinDE - 1.);

  Depo_ptr->pde_DEmax = Depo_ptr->pde_varDE[Depo_ptr->pde_nbinDE - 1] + .5*(
      Depo_ptr->pde_varDE[Depo_ptr->pde_nbinDE - 1] - Depo_ptr->pde_varDE[0])/(
          Depo_ptr->pde_nbinDE - 1.);

  Depo_ptr->pde_DPzmin = Depo_ptr->pde_varDPz[0] - .5*(
      Depo_ptr->pde_varDPz[Depo_ptr->pde_nbinDPz - 1] - Depo_ptr->pde_varDPz[0]) / (
          Depo_ptr->pde_nbinDPz - 1. );

  Depo_ptr->pde_DPzmax = Depo_ptr->pde_varDPz[Depo_ptr->pde_nbinDPz - 1] + .5*(
      Depo_ptr->pde_varDPz[Depo_ptr->pde_nbinDPz - 1] -
      Depo_ptr->pde_varDPz[0]) / (Depo_ptr->pde_nbinDPz - 1.);

  /* malloc and read the 5d distribution in */
  sz = Depo_ptr->pde_nbinE * Depo_ptr->pde_nbinPz * Depo_ptr->pde_nbinmu *
      Depo_ptr->pde_nbinDE * Depo_ptr->pde_nbinDPz;
  Depo_ptr->pde_pdedp = (double*)calloc(sz, sizeof(double));
  /*  the loop order is trecherous */
  for(je=0; je < Depo_ptr->pde_nbinE; je++){
    for(jp=0; jp < Depo_ptr->pde_nbinPz; jp++){
      for(jmu=0; jmu < Depo_ptr->pde_nbinmu; jmu++){
        for(jde=0; jde < Depo_ptr->pde_nbinDE; jde++){
          for(jdp=0; jdp < Depo_ptr->pde_nbinDPz; jdp++){
            i = Depo_ptr->pde_nbinPz * Depo_ptr->pde_nbinE * Depo_ptr->pde_nbinDPz * Depo_ptr->pde_nbinDE * jmu +
                Depo_ptr->pde_nbinE * Depo_ptr->pde_nbinDPz * Depo_ptr->pde_nbinDE * jp +
                Depo_ptr->pde_nbinDPz * Depo_ptr->pde_nbinDE * je +
                Depo_ptr->pde_nbinDE * jdp + jde;
            fscanf(ifp, "%lf ", &Depo_ptr->pde_pdedp[i]);
          }
        }
      }
    }
  }


  fclose(ifp);

  return;
}

void class_domain(Deposition_t* Depo_ptr){
  printf("class domain not implimented yet\n");
  return;
}

void pdedp_checkbdry(Deposition_t* Depo_ptr){
  return;
}

void pdedp_finalize(Deposition_t* Depo_ptr){
  return;
}

void pdedp_out(Deposition_t* Depo_ptr){
  return;
}

void pdedp_rcrd_resid(Deposition_t* Depo_ptr){
  return;
}

void rcrd_bfield(Deposition_t* Depo_ptr){
  return;
}

