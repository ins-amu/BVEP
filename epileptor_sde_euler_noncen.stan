data {
    int nt;
    real dt;
    real eta_true;
    real x_init;
    real z_init;
    real xlim[2];
    real zlim[2];
    real I1;
    real tau0;
    vector[nt] xs;
}

transformed data {
}

parameters {
    vector[nt-1] x_eta;
    vector[nt-1] z_eta;
    real x_init_star;
    real z_init_star;
    real eta;  
    real <lower=0.0> amplitude;
    real offset; 
    real<lower=0.0> eps;   
    real<lower=0.0> sig;
}

transformed parameters {
}

model {
    vector[nt] xhat;
    vector[nt] zhat;

    vector[nt] x;
    vector[nt] z;
             
    x_init_star ~ normal(0., 1.); 
    z_init_star ~ normal(0., 1.); 
        
    x[1] = x_init + x_init_star;
    z[1] = z_init + z_init_star;
                      
    to_vector(x_eta) ~ normal(0., 1.);
    to_vector(z_eta) ~ normal(0., 1.);
                         
   // eta ~ normal(eta_true, 1.); 

    amplitude ~ normal(1.,1.);
    offset ~ normal(0., 1.);
    eps ~ normal(0., 1.); 
    sig ~ normal(0., 1.);
 

    for (t in 1:(nt-1)) {
            real dx = 1.0 - x[t]*x[t]*x[t] - 2.0*x[t]*x[t] - z[t] + I1;
            real dz = (1.0/tau0)*(4*(x[t] - eta) - z[t] );
            x[t+1] =x[t] + dt*dx+ sqrt(dt)*x_eta[t]*sig;
            z[t+1] =z[t] + dt*dz+ sqrt(dt)*z_eta[t]*sig;  
    }  
    
    xhat=amplitude*x+offset;
    zhat=amplitude*z+offset;

    //xs ~ normal(xhat, eps); 

    target+=normal_lpdf(xs| xhat, eps);   

}


generated quantities {

    vector[nt] xhat_q;
    vector[nt] zhat_q;
    vector[nt] x_ppc;
    vector[nt] z_ppc;
    vector[nt] log_lik;

    vector[nt] x;
    vector[nt] z;

    x[1] = x_init + x_init_star;
    z[1] = z_init + z_init_star;

    for (t in 1:(nt-1)) {
            real dx = 1.0 - x[t]*x[t]*x[t] - 2.0*x[t]*x[t] - z[t] + I1;
            real dz = (1.0/tau0)*(4*(x[t] - eta) - z[t] );
            x[t+1] =x[t] + dt*dx+ sqrt(dt)*x_eta[t]*sig;
            z[t+1] =z[t] + dt*dz+ sqrt(dt)*z_eta[t]*sig;  
    }    

    xhat_q=amplitude*x+offset;
    zhat_q=amplitude*z+offset;

    for (i in 1:nt){
        x_ppc[i] = normal_rng(xhat_q[i], eps);
        z_ppc[i] = normal_rng(zhat_q[i], eps);

      }
      
    for (i in 1:nt){
        log_lik[i] = normal_lpdf(xs[i]| xhat_q[i], eps);
      }

      
}
