#include "orbit_config_api.h"


void do_particles(Config_t* cfg_ptr){
    // before we begin, we need to compute "eps"
  const double eps = compute_eps(ped, *pw,
                                 *lsp, *lst, *idp,
                                 x1, x2, x3,
                                 x4, x5, x6,
                                 x7, x8, x9);



  int particle_id;
  int ktm;

  for(particle_id=0; particle_id < nprt; particle_id++){
    for(ktm=1; ktm < nstep_all; ktm++){

      konestep(particle_id, cfg_ptr);

      kupdate(particle_id, cfg_ptr);

    }
  }
}

void konestep(int k, Config_t cfg_ptr){

  int n1,j,i,ndum;
  double xdum,ydum,rbb,dedb,deni,fac1,fac2,pdum,xdot,ydot;
  /* local temps */
  double y_[4];
  double d_[4];
  double e_[4];
  double a_[4];
  double bx_[4];
  double h_;
  double c1_[4];
  double zdots_;

  n1=4;

  dt[k]=dt0;

  nout[k] = .6 * (1. + sign(1., pol[k]- pw) );
  nfin[k] =  .6 * (1. + sign(1., time[k] - trun));
  nt_[k] = .6 * (1.  + sign(1., pol[k]-.05 *pw));
  dt[k] = nt_[k] * dt[k] + (1 - nt_[k]) * dt0;

  xdum = sqrt(pol[k])*cos(thet[k]);
  ydum = sqrt(pol[k])*sin(thet[k]);


  y_[0] = nt_[k] * pol[k] + (1-nt_[k])*xdum;
  y_[1] = nt_[k] * thet[k] + (1-nt_[k])*ydum;
  y_[2] = zet[k];
  y_[3] = rho[k];

  // XXXX gbw, check this, looks like i typo'd
  d_[0] = y_[0];
  d_[1] = y_[1];
  d_[2] = y_[2];
  d_[3] = y_[3];
  h_ = dt[k]/6.0 ;

  for(j=0; j < 4; j++){
    kfield(cfg_ptr, Ptcl, 1);


    if(npert ==  0){
      //goto 61
      if(k == 0){
        printf("FAILURE: integ without perturbations is not yet implimented\n");
      }
      return;
    }

    rbb = y_[3] * b[k]*b[k];
    dedb = b[k] * y_[3]*y_[3] + rmu[k];
    deni = 1. / ( g[k]*q[k] + ri[k] + (chrg * y_[3]+alp[k]) * (g[k]*rip[k] - ri[k]*gp[k]));
    fac1 = 1 - gp[k]*(chrg*y_[3] + alp[k]) - g[k]*dadp[k];
    fac2 = q[k] + rip[k]*(chrg*y_[3] + alp[k]) + ri[k]*dadp[k];
    //   pol dot


    e_[0] = (-ri[k]*rbb*dadz[k]
             - chrg*g[k]*dedb*dbdt[k]
             + g[k]*rbb*dadt[k]
             - chrg*g[k]*dptdt[k]
             + chrg*ri[k]*dptdz[k]
             + chrg*ri[k]*dedb*dbdz[k]) * deni;
    //   thet dot
    e_[1] = (chrg*dedb*dbdp[k]*g[k]+ rbb*fac1 + chrg*g[k]*dptdp[k])*deni;
    //   zet dot
    e_[2] = (-chrg*dedb*dbdp[k]*ri[k] + rbb*fac2 - chrg*ri[k]*dptdp[k])*deni;


    ydot = .5 *y_[1] * e_[0]/pdum + y_[0] *e_[1];
    e_[0] = nt_[k]*e_[0] + (1-nt_[k])*xdot;
    e_[1] = nt_[k]*e_[1] + (1-nt_[k])*ydot;
    e_[0] = e_[0] * (1-nout[k])*(1-nfin[k]);
    e_[1] = e_[1] * (1-nout[k])*(1-nfin[k]);
    e_[2] = e_[2] * (1-nout[k])*(1-nfin[k]);
    e_[3] = e_[3] * (1-nout[k])*(1-nfin[k]);

    //goto 62, like 42, the answere to the universe and everything in it,  but twenty years older.
    for(i=0; i<n1; i++){
      if(j==0){
        //for(i=0; i<n1; i++){
        a_[i] = h_[k] * e_[i];
        y_[i] = d_[i] + 3 * a_[i];
        //}
      }
      if(j==1){
        // for(i=0; i<n1; i++){
        bx_[i] = h_[k] * e_[i];
        y_[i] = d_[i] + 3 * bx_[i];
        //}
      }
      if(j==2){
        //for(i=0; i<n1; i++){
        c1_[i] = h_[k] * e_[i];
        y_[i] = d_[i] + 6 * c1_[i];
        //}
      }
      if(j==3){
        //for(i=0; i<n1; i++){
        y_[i] = d_[i] + a_[i] + 2 * bx_[i] + 2 * c1_[i] + h_[k] * e_[i];
        //}
      }
    }
    // 40
    ndum = .6 * (1 - sign(1. , y_[0]));
    pol[k] = nt_[k]*y_[0] + (1-nt_[k])*(y_[0]*y_[0] + y_[1]*y_[1]);
    thet[k] =nt_[k]*y_[1] + (1-nt_[k])*(atan(y_[1]/y_[0]) + ndum * pi);
    zet[k] = y_[2];
    rho[k] = y_[3];
    nout[k] = .6 * (1. + sign(1., pol[k]-pw));

  } //j
  return;
}

