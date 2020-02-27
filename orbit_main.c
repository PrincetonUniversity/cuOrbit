/*  Copyright 2019, 2020 Garrett Wright, Princeton Plasma Physic Lab

    This file is part of CuOrbit.

    CuOrbit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CuOrbit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CuOrbit.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

#ifdef __NVCC__
#include <cuda.h>
#include <cuda_runtime.h>
#endif

#include "orbit_config_api.h"
#include "orbit_util.h"
#include "orbit_perturbation.h"
#include "orbit_equilibrium.h"
#include "orbit_particles.h"
#include "orbit_deposition.h"
#include "orbit_main.h"
#include "cuda_helpers.h"

void orbit_main_loop(orbit_Config_t* cfg_ptr){
  int irun_pdedp;

  Deposition_t* depo_ptr = cfg_ptr->depo_ptr;

  for(irun_pdedp=0; irun_pdedp < cfg_ptr->nruns; irun_pdedp++){
    if(irun_pdedp>1 &&
       irun_pdedp > cfg_ptr->nruns/2 &&
       get_pdedp_focusdep(cfg_ptr->depo_ptr))
    {
      printf("Enabling pdedp_do_focusdep\n");
      set_pdedp_do_focusdep(cfg_ptr->depo_ptr, true);
    }

    /* if(irun_pdedp>0){ */
    /*   printf("DBG exit after initial batch\n"); */
    /*   exit(0); */
    /* } */

    /* if we are computing pDEDP and it is the first run */
    if(compute_pdedp(cfg_ptr->depo_ptr) && irun_pdedp == 0){
      if(get_initial_update_pdedp_from_file(cfg_ptr->depo_ptr)){
        printf(" ... reading pDEDP from ufile...\n");
        pdedp_read(cfg_ptr->depo_ptr, cfg_ptr);
        printf(" ... done.\n");
      } else{
        /* compute new p(DE,DP), zero inits pdedp et al */
        pdedp_init(cfg_ptr->depo_ptr);
      }
    }

    /* Call Mario's deposition, or Roscoe's file read..*/
    deposition(cfg_ptr, irun_pdedp);

    
    /* MP: I moved the following definitions to orbit_deposition.c */
    /* double dum = 1E3 * cfg_ptr->dt0 / get_omeg0(cfg_ptr->ptrb_ptr);
    const int nstep_all = round(get_pdedp_dtrun(depo_ptr) / dum);
    cfg_ptr->nstep_all = nstep_all;
    set_pdedp_tskip(depo_ptr,
                    imax(round(nstep_all / 2E4) + 1,
                         1));  */  /* stay in array bounds */  
    
    printf("\n\n --- Start main run --- \n" );
    printf("\t no. of particles \t: %d\n", cfg_ptr->nprt);
    printf("\t no. of time steps \t: %d\n", cfg_ptr->nstep_all);
    printf("\t sim. time [ms] \t:  %f\n", cfg_ptr->pdedp_dtrun);
    printf("\t time step [us] \t:  %f\n",
           1E6 * cfg_ptr->dt0 / get_omeg0(cfg_ptr->ptrb_ptr));
    printf("\t time steps to skip \t:  %d\n\n", get_pdedp_tskip(depo_ptr));

    /* Check Run parameters, exit if issues are detected */
    if( cfg_ptr->nstep_all <= 1){
      printf("Error: number of simulation steps is %i," \
             " Please set pdedp_dtrun in config.ini. Aborting\n", \
	     cfg_ptr->nstep_all);
      exit(1);
    }
    
    /* launch the stepping functions */
    do_particles(cfg_ptr);

    /* end of main run*/

    if (compute_pdedp(depo_ptr)){
      /* this code was in the loop, needto  investigate how to run all at once */
      /* apparently this needs triple testing */
      pdedp_rcrd_resid(cfg_ptr, depo_ptr);
      if (pdedp_optimize(depo_ptr)){
        pdedp_checkbdry(cfg_ptr, depo_ptr);
      }

      // debug out
      //pdedp_out(depo_ptr);
    }
    printf("- p(DE,DP) calculations: done %d / %d \n", irun_pdedp+1, cfg_ptr->nruns);
    printf("\n ------------------------------ \n\n");
  }  /* irun_pdedp */

  /* this was out of loop */
  if (compute_pdedp(depo_ptr)){
    /* normalize and fill empty bins */
    pdedp_finalize(depo_ptr);
    /* write output file *AEP */
    pdedp_out(depo_ptr);
  }

  rcrd_bfield(cfg_ptr, depo_ptr);

  /* i think this is just diagnostic output */
  /* wrt6();  */
  /* phys6(); */

}
