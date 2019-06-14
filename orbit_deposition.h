#ifndef SET_ORBIT_DEPOSITION_H_
#define SET_ORBIT_DEPOSITION_H_

#include "orbit_config.h"

typedef struct Deposition Deposition_t;
void initialize_Deposition(Deposition_t*, Config_t*);
Deposition_t* Deposition_ctor();

bool compute_pdedp(Deposition_t*);
bool initial_update_pdedp(Deposition_t*);
bool pdedp_optimize(Deposition_t*);

int get_pde_dtsamp(Deposition_t*);
int get_pde_tskip(Deposition_t*);
void set_pde_tskip(Deposition_t*, double);

void pdedp_read(Deposition_t*);
void pdedp_checkbdry(Deposition_t*);
void pdedp_finalize(Deposition_t*);
void pdedp_out(Deposition_t*);
void pdedp_rcrd_resid(Deposition_t*);
void rcrd_bfield(Deposition_t*);

#endif
