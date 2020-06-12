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
    vector<lower=xlim[1], upper=xlim[2]>[nt] x;
    vector<lower=zlim[1], upper=zlim[2]>[nt] z;
    real eta;  
    real <lower=0.0> amplitude;
    real offset; 
    real<lower=0.0> eps;   
    real<lower=0.0> sig;
    
}

model {
    vector[nt] xhat;
    vector[nt] zhat;

    x[1] ~ normal(x_init, 1.);
    z[1] ~ normal(z_init, 1.);
                         
    //eta ~ normal(eta_true, 1.); 

    amplitude ~ normal(1.,1.);
    offset ~ normal(0., 1.);
    eps ~ normal(0., 1.); 
    sig ~ normal(0., 1.);


    for (t in 1:(nt-1)) {
            real dx = 1.0 - x[t]*x[t]*x[t] - 2.0*x[t]*x[t] - z[t] + I1;
            real dz = (1.0/tau0)*(4*(x[t] - eta) - z[t] );
            x[t+1] ~ normal(x[t] + dt*dx, sqrt(dt)*sig); 
            z[t+1] ~ normal(z[t] + dt*dz, sqrt(dt)*sig);    
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