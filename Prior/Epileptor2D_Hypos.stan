data {
    int nt;
    real dt;
    vector[nt] xs;
    real xlim[2];
    real zlim[2];
    real I1;
    real tau0;
    real eta_sd; 
    real eta_mu;
}

transformed data {
}

parameters {
    vector<lower=xlim[1], upper=xlim[2]>[nt] x;
    vector<lower=zlim[1], upper=zlim[2]>[nt] z;
    real eta;  
    real <lower=0.0> amp;
    real offset; 
    real<lower=0.0> eps;   
    real<lower=0.0> sig;
}

model {
    vector[nt] xhat;

    x[1] ~ normal(-1.5, 1.);
    z[1] ~ normal(3.5, 1.);
                         
    eta ~ normal(eta_mu, eta_sd); 
    amp ~ normal(1,2);
    offset ~ normal(0, 2);
    eps ~ normal(0., 1.); 
    sig ~ normal(0., 1.);


    for (t in 1:(nt-1)) {
            real dx = 1.0 - x[t]*x[t]*x[t] - 2.0*x[t]*x[t] - z[t] + I1;
            real dz = (1.0/tau0)*(4*(x[t] - eta) - z[t] );
            x[t+1] ~ normal(x[t] + dt*dx, sqrt(dt)*sig); 
            z[t+1] ~ normal(z[t] + dt*dz, sqrt(dt)*sig);    
    }  

    xhat=amp*x+offset;

    //xs ~ normal(xhat, eps); 

    target+=normal_lpdf(xs| xhat, eps);    
}

generated quantities {

    vector[nt] xhat_q;
    vector[nt] x_ppc;
    vector[nt] log_lik;


    int num_params;
    int num_data;
    num_params=2*nt+5;
    num_data=nt;

    xhat_q=amp*x+offset;


    for (i in 1:nt){
        x_ppc[i] = normal_rng(xhat_q[i], eps);
      }
      
    for (i in 1:nt){
        log_lik[i] = normal_lpdf(xs[i]| xhat_q[i], eps);
      }

}