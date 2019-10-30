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

void main_loop(Config_t* cfg_ptr){
  int irun_pdedp;

  Deposition_t* depo_ptr = cfg_ptr->depo_ptr;

  for(irun_pdedp=0; irun_pdedp < cfg_ptr->nruns; irun_pdedp++){
    if(irun_pdedp>0){
      set_initial_update_pdedp(cfg_ptr->depo_ptr, true);
    }
    if(irun_pdedp>1 && irun_pdedp > cfg_ptr->nruns/2){
      set_pdedp_focusdep(cfg_ptr->depo_ptr, true);
    }


    if(compute_pdedp(cfg_ptr->depo_ptr)){
      if(initial_update_pdedp(cfg_ptr->depo_ptr) && (irun_pdedp == 0)){
        printf(" ... reading pDEDP from ufile...\n");
        pdedp_read(cfg_ptr->depo_ptr, cfg_ptr);
        printf(" ... done.\n");
      } else {
        /* compute new p(DE,DP) */
        pdedp_init(cfg_ptr->depo_ptr);
      }
    }

    /* XXXXX, yeah something off here */
    /* for now we just  use the fulldepmp per mario */
    /* if( irun_pdedp % 2  == 0){ */
    /*   fulldepmp(cfg_ptr, depo_ptr); */
    /* } else { */
    /*   fulldepmp_co(cfg_ptr, depo_ptr); */
    /* } */
    fulldepmp(cfg_ptr, depo_ptr);


    double dum = 1E3 * cfg_ptr->dt0 / get_omeg0(cfg_ptr->ptrb_ptr);
    const int nstep_all = round(10. * get_pdedp_dtsamp(depo_ptr) / dum) + 1;
    cfg_ptr->nstep_all = nstep_all;
    set_pdedp_tskip(depo_ptr,
                    imax(round(nstep_all / 1E4) + 1,
                         5));    /* stay in array bounds */

    printf("\n\n --- Start main run --- \n" );
    printf("\t no. of particles \t: %d\n", cfg_ptr->nprt);
    printf("\t no. of time steps \t: %d\n", nstep_all);
    printf("\t sim. time [ms] \t:  %f\n", nstep_all*dum);
    printf("\t time step [us] \t:  %f\n\n",
           1E6 * cfg_ptr->dt0 / get_omeg0(cfg_ptr->ptrb_ptr));
    ///XXXX time steps to skip

    /* launch the stepping functions */
    do_particles(cfg_ptr);

    /* end of main run*/

    if (compute_pdedp(depo_ptr)){

      /* this code was a in the loop, need to  investigate how to run all at once */
      /* apprently this needs triple testing */
      pdedp_rcrd_resid(cfg_ptr, depo_ptr);
      if (pdedp_optimize(depo_ptr)){
        pdedp_checkbdry(cfg_ptr, depo_ptr);
      }

      pdedp_out(depo_ptr);
      printf("- p(DE,DP) calculations: done\n");
    }
    printf("- p(DE,DP) calculations: done %d / %d \n", irun_pdedp+1, cfg_ptr->nruns);

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
