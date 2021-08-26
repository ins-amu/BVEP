data {
  int nn;  //number of brain regions
  int nt;  //number of data points per sensor
  int ns;  //number of sensors

  real dt;    //time step in Euler integration
  real tau0;  //time scale in Eipleptor model
  real I1;    //input current in in Eipleptor model
  real Ks;    //global connectivity parameter

  matrix[nt, ns] Obs_seeg;      // observation at sensor-level
  matrix[nn, nn] SC;           // brain structural connectivity
  matrix[ns, nn] Gr;           // gain (projection matrix)
  matrix[nn, nn] eigen_vec;   // eigen vectors of dot(Gr.T, Gr)
 
  real zlim[2]; 
  real xlim[2]; 

}

transformed data {
  vector[nt*ns] Obs_seeg_vect;
  Obs_seeg_vect=to_vector(Obs_seeg);
}

parameters {
  vector[nn] x_init_star;
  vector[nn] z_init_star;
  vector[nn] eta_star;
  real <lower=-1.0> K_star;
  real <lower=-1.0> amplitude_star;
  vector [ns] offset_star;
  real eps_star;
}

transformed parameters {
  vector[nn] eta;
  vector[nn] x_init;
  vector[nn] z_init;
  real K;
  real amplitude;
  vector[ns] offset;
  real eps;
  
  eta =  -3.5+0.1*eigen_vec*eta_star;
  x_init= -2.0 +0.1*eigen_vec*x_init_star;
  z_init= 5.0+0.1*eigen_vec*z_init_star;
  K =  Ks + 0.1*K_star;
  amplitude = 1.0 + amplitude_star;
  offset=offset_star+10.0;
  eps=exp(0.33*eps_star);
}

model {
  matrix[nn, nt] x;
  matrix[nn, nt] z;
  matrix[nt, ns] Seeg;
  vector[nt*ns] Seeg_vect;
  
  real dx;
  real dz;
  real gx;
  
  /* priors*/
  
  eta_star ~ normal(0.,1.0);  
  x_init_star ~ normal(0.,1.0);
  z_init_star ~ normal(0, 1.0);
  K_star ~ normal(0., 1.);  
  amplitude_star ~ normal(0. ,1.);
  offset_star ~ normal(0.,1.);
  eps_star ~ normal(0.,1.);

  
  /* integrate & predict */
    
    for (i in 1:nn) {
      x[i, 1] =  x_init[i];
      z[i, 1] =  z_init[i];
    } 
  
  for (t in 1:(nt-1)) {
    for (i in 1:nn) {
      gx = 0;
      for (j in 1:nn)
        gx = gx + SC[i, j]*(x[j, t] - x[i, t]);
      dx = 1.0 - x[i, t]*x[i, t]*x[i, t] - 2.0*x[i, t]*x[i, t] - z[i, t] + I1;
      dz = (1/tau0)*(4*(x[i, t] - eta[i]) - z[i, t] - K*gx);
      x[i, t+1] = x[i, t] + dt*dx ; 
      z[i, t+1] = z[i, t] + dt*dz ; 
    }
  }   
  
  for(i in 1:nt){
    Seeg[i,]=to_row_vector(amplitude*(Gr*x[,i])+offset);
  }
  
  Seeg_vect=to_vector(Seeg);
  
  for(i in (10*ns+1):(nt*ns)){
    target+=  normal_lpdf(Obs_seeg_vect[i]| Seeg_vect[i], eps);  
  }
}  

generated quantities {
  
  matrix[nn, nt] x;
  matrix[nn, nt] z;
  
  matrix[nt,ns] Seeg_qqc;
  vector[nt*ns] Seeg_qqc_vect;
  vector[nt*ns] Seeg_ppc;
  vector[nt*ns] log_lik;
  vector[nt] log_lik_sum = rep_vector(0,nt);  
  
  real gx;
  real dx;
  real dz;
  
  
  for (i in 1:nn) {
    x[i, 1] = x_init[i];
    z[i, 1] = z_init[i];
  } 
  
  for (t in 1:(nt-1)) {
    for (i in 1:nn) {
      gx = 0;
      for (j in 1:nn)
        gx = gx + SC[i, j]*(x[j, t] - x[i, t]);
      dx = 1.0 - x[i, t]*x[i, t]*x[i, t] - 2.0*x[i, t]*x[i, t] - z[i, t] + I1;
      dz = (1/tau0)*(4*(x[i, t] - eta[i]) - z[i, t] - K*gx);
      x[i, t+1] = x[i, t] + dt*dx ; 
      z[i, t+1] = z[i, t] + dt*dz ; 
    }
  }  
  
  
  
  for(i in 1:nt){
    Seeg_qqc[i,]=to_row_vector(amplitude*(Gr*x[,i])+offset);
  }
  
  
  Seeg_qqc_vect=to_vector(Seeg_qqc);

  for (i in 1:(nt*ns)){
        Seeg_ppc[i] = normal_rng(Seeg_qqc_vect[i], eps);
      }


  for (i in 1:(nt*ns)){
        log_lik[i] = normal_lpdf(Obs_seeg_vect[i]| Seeg_qqc_vect[i], eps);
      }  
  

  for (i in 1:nt){
        for (j in 1:ns)
            log_lik_sum[i] += normal_lpdf(Obs_seeg[i,j]| Seeg_qqc[i,j], eps);
      }   

}
