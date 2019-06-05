#include <stdlib.h>
#include <stdbool.h>

#include "orbit_deposition.h"

typedef struct Deposition {
  /* xxx stochastic, does this belong here? */
  double mubk;
  int emink;
  int emaxk;
  int dmubk;
  int nstoche;
  int nstochp;

  /* pdedp */
  int nruns;
  bool compute_pdedp;
  bool initial_update_pdedp;
  double deposit_on_bins_after_fraction;
  double pde_dtsamp;
  double pde_dtav;
  int pde_tskip;
  double pde_otpup;
  bool pde_focusdep;
  bool pde_optimize;

} Deposition_t;

Deposition_t* Deposition_ctor(){
  return (Deposition_t*)calloc(1, sizeof(Deposition_t));
}

void initialize_Deposition(Deposition_t* Depo_ptr, Config_t* config_ptr){
  return;
}

