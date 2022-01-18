functions {
  vector ode_rhs(real time, vector xz, matrix SC, real I1, real tau0, real K, vector eta) {
    int nn = rows(xz)/2;
    vector[nn] x = xz[1:nn];
    vector[nn] z = xz[nn+1:2*nn];
    vector[nn] gx = SC * x;
    vector[nn] dx = 1.0 - x.*x.*x - 2.0*x.*x - z + I1;
    vector[nn] dz = (1/tau0)*(4*(x - eta) - z - K*gx);
    return append_row(dx, dz);
  }
}

data {
  int nn;  //number of brain regions
  int nt;  //number of data points per sensor
  int ns;  //number of sensors

  real dt;    //time step in Euler integration
  real tau0;  //time scale in Eipleptor model
  real I1;    //input current in in Eipleptor model
  real Ks;    //global connectivity parameter

  matrix[nt, ns] Obs_seeg;      // observation at sensor-level
  matrix[nn, nn] SC_star;           // brain structural connectivity
  matrix[ns, nn] Gain;           // gain (projection matrix)
  matrix[nn, nn] eigen_vec;   // eigen vectors of dot(Gain.T, Gain)

  real amplitude_fixed;
  real offset_offset;
  real offset_tol;
  real eps_tol;
  real eta_spread;
  real eta_offset; 
 
  // ode solving
  int ode_solver;
  real ode_rtol;
  real ode_atol;
  int ode_maxstep;
}

transformed data {
  vector[nt*ns] Obs_seeg_vect;
  Obs_seeg_vect=to_vector(Obs_seeg);
  vector[2*nn] ode_atol_fwd = rep_vector(ode_atol, 2*nn);
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
  
  eta =  eta_offset+eta_spread*eigen_vec*eta_star; 
  x_init= -2.0 +0.1*eigen_vec*x_init_star;
  z_init= 5.0+0.1*eigen_vec*z_init_star;
  K =  Ks + 0.1*K_star;
  amplitude = amplitude_fixed; // 1.0 + amplitude_star;
  offset=offset_tol * offset_star + offset_offset;
  eps=exp(eps_tol*eps_star); 
}

model {
  matrix[nn, nt] x;
  matrix[nn, nt] z;
  matrix[nt, ns] Seeg;
  vector[nt*ns] Seeg_vect;
  
  vector[nn] dx;
  vector[nn] dz;
  vector[nn] gx;
  
  /* priors*/
  
  eta_star ~ normal(0.,1.0);  
  x_init_star ~ normal(0.,1.0);
  z_init_star ~ normal(0, 1.0);
  K_star ~ normal(0., 1.);  
  amplitude_star ~ normal(0. ,1.);
  offset_star ~ normal(0.,1.);
  eps_star ~ normal(0.,1.);

  
  x[,1] = x_init;
  z[,1] = z_init;
  if (ode_solver == 0) {
    // this is Euler. 75% of the grad eval time
    for (t in 1:(nt - 1)) {
      gx = SC_star * x[,t];
      dx = 1.0 - x[,t].*x[,t].*x[,t] - 2.0*x[,t].*x[,t] - z[,t] + I1;
      dz = (1/tau0)*(4*(x[,t] - eta) - z[,t] - K*gx);
      x[,t+1] = x[,t] + dt*dx ; 
      z[,t+1] = z[,t] + dt*dz ; 
    }
  } else if (ode_solver == 1) {
    // this is Heun
    for (t in 1:(nt - 1)) {
      gx = SC_star * x[,t];
      dx = 1.0 - x[,t].*x[,t].*x[,t] - 2.0*x[,t].*x[,t] - z[,t] + I1;
      dz = (1/tau0)*(4*(x[,t] - eta) - z[,t] - K*gx);
      vector[nn] xi = x[,t] + dt*dx;
      vector[nn] zi = z[,t] + dt*dz;
      vector[nn] gxi = SC_star * xi;
      vector[nn] dxi = 1.0 - xi.*xi.*xi - 2.0*xi.*xi - zi + I1;
      vector[nn] dzi = (1/tau0)*(4*(xi - eta) - zi - K*gxi);
      x[,t+1] = x[,t] + dt/2*(dx + dxi); 
      z[,t+1] = z[,t] + dt/2*(dz + dzi); 
    }
  } else if (ode_solver == 4) {
    // this is RK4
    for (t in 1:(nt - 1)) {
      // k1 = f(y)
      dx = 1.0 - x[,t].*x[,t].*x[,t] - 2.0*x[,t].*x[,t] - z[,t] + I1;
      dz = (1/tau0)*(4*(x[,t] - eta) - z[,t] - K*SC_star*x[,t]);
      // k2 = f(y + dt*k1/2)
      vector[nn] xi = x[,t] + dt*dx/2;
      vector[nn] zi = z[,t] + dt*dz/2;
      vector[nn] dxi = 1.0 - xi.*xi.*xi - 2.0*xi.*xi - zi + I1;
      vector[nn] dzi = (1/tau0)*(4*(xi - eta) - zi - K*SC_star*xi);
      // k3 = f(y + dt*k2/2)
      vector[nn] xii = x[,t] + dt*dxi/2;
      vector[nn] zii = z[,t] + dt*dzi/2;
      vector[nn] dxii = 1.0 - xii.*xii.*xii - 2.0*xii.*xii - zii + I1;
      vector[nn] dzii = (1/tau0)*(4*(xii - eta) - zii - K*SC_star*xii);
      // k4 = f(y + dt*k3)
      vector[nn] xiii = x[,t] + dt*dxii;
      vector[nn] ziii = z[,t] + dt*dzii;
      vector[nn] dxiii = 1.0 - xiii.*xiii.*xiii - 2.0*xiii.*xiii - ziii + I1;
      vector[nn] dziii = (1/tau0)*(4*(xiii - eta) - ziii - K*SC_star*xiii);
      x[,t+1] = x[,t] + dt/6*(dx + dxi*2 + dxii*2 + dxiii); 
      z[,t+1] = z[,t] + dt/6*(dz + dzi*2 + dzii*2 + dziii); 
    }
  } else if (ode_solver == 2) {
    
    vector[2*nn] ic = append_row(x_init, z_init);
    real it = 0.0;
    real ts[nt-1];
    for (t in 2:nt) ts[t-1] = t*dt;
    vector[2*nn] ode_sol[nt-1] = ode_ckrk_tol(
      ode_rhs, ic, it, ts, ode_rtol, ode_atol, ode_maxstep,
      SC_star, I1, tau0, K, eta);
    for (t in 2:nt) x[,t] = ode_sol[t-1][1:nn];
  } else if (ode_solver == 3) {
    vector[2*nn] ic = append_row(x_init, z_init);
    real it = 0.0;
    real ts[nt-1];
    for (t in 2:nt) ts[t-1] = t*dt;
    vector[2*nn] ode_sol[nt-1] = ode_adjoint_tol_ctl(
      ode_rhs, ic, it, ts, ode_rtol, ode_atol_fwd, ode_rtol, ode_atol_fwd,
      ode_rtol, ode_atol, ode_maxstep, 100, 1, 1, 1,
      SC_star, I1, tau0, K, eta);
    for (t in 2:nt) x[,t] = ode_sol[t-1][1:nn];
  }
  
  for(i in 1:nt){
    Seeg[i,]=to_row_vector(amplitude*(Gain*x[,i])+offset);
  }
  
  Seeg_vect=to_vector(Seeg);
  target += normal_lpdf(Obs_seeg_vect | Seeg_vect, eps);
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
        gx = gx + SC_star[i, j]*(x[j, t] - x[i, t]);
      dx = 1.0 - x[i, t]*x[i, t]*x[i, t] - 2.0*x[i, t]*x[i, t] - z[i, t] + I1;
      dz = (1/tau0)*(4*(x[i, t] - eta[i]) - z[i, t] - K*gx);
      x[i, t+1] = x[i, t] + dt*dx ; 
      z[i, t+1] = z[i, t] + dt*dz ; 
    }
  }  
  
  
  
  for(i in 1:nt){
    Seeg_qqc[i,]=to_row_vector(amplitude*(Gain*x[,i])+offset);
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
