#include <stdlib.h>
#include <stdbool.h>

#include "orbit_config_api.h"
#include "orbit_deposition.h"
#include "orbit_particles.h"  /* ekev */

typedef struct Deposition {
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

  Depo_ptr->nruns = cfg_ptr->nruns;
  Depo_ptr->compute_pdedp = cfg_ptr->compute_pdedp;
  Depo_ptr->initial_update_pdedp = cfg_ptr->initial_update_pdedp;
  Depo_ptr->deposit_on_bins_after_fraction = cfg_ptr->deposit_on_bins_after_fraction;
  Depo_ptr->pde_dtsamp = cfg_ptr->pde_dtsamp;
  Depo_ptr->pde_dtav = cfg_ptr->pde_dtav;
  Depo_ptr->pde_tskip = cfg_ptr->pde_tskip;
  Depo_ptr->pde_otpup = cfg_ptr->pde_otpup;
  Depo_ptr->pde_focusdep = cfg_ptr->pde_focusdep;
  Depo_ptr->pde_optimize = cfg_ptr->pde_optimize;

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
  return Depo_ptr->pde_optimize;
}


int get_pde_dtsamp(Deposition_t* Depo_ptr){
  return Depo_ptr->pde_dtsamp;
}

int get_pde_tskip(Deposition_t* Depo_ptr){
  return Depo_ptr->pde_tskip;
}
void set_pde_tskip(Deposition_t* Depo_ptr, double pde_tskip){
  Depo_ptr->pde_tskip = pde_tskip;
  return;
}



void pdedp_read(Deposition_t* Depo_ptr){

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

