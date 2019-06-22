#ifndef SET_ORBIT_DEPOSITION_H_
#define SET_ORBIT_DEPOSITION_H_

#include "orbit_config.h"

typedef struct Deposition Deposition_t;
void initialize_Deposition(Deposition_t*, Config_t*);
Deposition_t* Deposition_ctor();
void initialize_pdedp(Deposition_t*);

bool compute_pdedp(Deposition_t*);
bool initial_update_pdedp(Deposition_t*);
bool pdedp_optimize(Deposition_t*);

int get_pdedp_dtsamp(Deposition_t*);
int get_pdedp_tskip(Deposition_t*);
void set_pdedp_tskip(Deposition_t*, double);
bool get_initial_update_pdedp(Deposition_t*);
void set_initial_update_pdedp(Deposition_t*, bool);
bool get_pdedp_focusdep(Deposition_t*);
void set_pdedp_focusdep(Deposition_t*, bool);

size_t sizeof_pdedp(Deposition_t*);

void pdedp_read(Deposition_t*, Config_t*);
void pdedp_init(Deposition_t*);
void fulldepmp(Config_t*, Deposition_t*);
void fulldepmp_co(Config_t*, Deposition_t*);
void fullredepmp(Config_t* , Deposition_t*);
void class_kdomain(Config_t*, int);
void class_domain(Config_t*);
void pdedp_checkbdry(Config_t*, Deposition_t*);
void pdedp_finalize(Deposition_t*);
void pdedp_out(Deposition_t*);
void pdedp_rcrd_resid(Config_t*, Deposition_t*);
void rcrd_bfield(Config_t*, Deposition_t*);
void check_res_ptc(Config_t*, int);

#endif
