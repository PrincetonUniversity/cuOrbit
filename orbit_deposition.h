#ifndef SET_ORBIT_DEPOSITION_H_
#define SET_ORBIT_DEPOSITION_H_

#include "orbit_config.h"

typedef struct Deposition Deposition_t;
void initialize_Deposition(Deposition_t*, Config_t*);
Deposition_t* Deposition_ctor();


#endif
