/* Written in INS Marseille*/
data {
    int nn; \\number of brain regions
    int nt;  \\number of data points per brain region
    real dt; \\time step in Euler integration
    real tau0; \\time scale in Eipleptor model
    real I1; \\input current in in Eipleptor model
    real Ks; \\global connectivity parameter
    matrix [nn, nn] SC; \\ brain structural connectivity
    matrix[nn, nt] Obs; \\fast activity variable at source-level
}

transformed data {
    vector[2] initial_val;
    vector[nn*nt] xs;
    xs=to_vector(Obs); \\vectorize observations
    initial_val[1]=-1.5; \\initial condition of fast activity variable in Epileoptor
    initial_val[2]=+3.5; \\initial condition of Epileoptor slow variable in Epileoptor
}

parameters {
    vector[nn] x_init;
    vector[nn] z_init;
    matrix[nn, nt-1] x_eta;
    matrix[nn, nt-1] z_eta;
    vector[nn] eta_star;
    real <lower=0.0> K_star;
    real <lower=0.0> amplitude_star;
    real offset; 
    real<lower=0.0> eps;
    real<lower=0.0> sig; 
}

transformed parameters {
    vector[nn] eta;
    real K;
    real amplitude;

    eta = -2.5 + eta_star;
    K = Ks+ K_star;
    amplitude = 1.0 + amplitude_star;
}

model {
    matrix[nn, nt] x;
    matrix[nn, nt] z;
    vector[nn*nt] xhat;
    real dx;
    real dz;
    real gx;

    /* priors*/

    eta_star ~ normal(0.0, 1.0);  
    amplitude_star ~ normal(0.0, 1.0);
    offset ~ normal(0.0, 1.0);
    K_star ~ normal(0.0, 1.0);  
    eps ~ normal(0.0, 1.0);   
    sig ~ normal(0.0, 1.0);   

    to_vector(x_eta) ~ normal(0.0, 1.0);
    to_vector(z_eta) ~ normal(0.0, 1.0);
    
    x_init ~ normal(0., 1.0);
    z_init ~ normal(0., 1.0);
    
    /* integrate & predict */

    for (i in 1:nn) {
      x[i, 1] = initial_val[1] + x_init[i];
      z[i, 1] = initial_val[2] + z_init[i];
    } 
    
    for (t in 1:(nt-1)) {
        for (i in 1:nn) {
            gx = 0;
            for (j in 1:nn)
                    gx = gx + SC[i, j]*(x[j, t] - x[i, t]);
            dx = 1.0 - x[i, t]*x[i, t]*x[i, t] - 2.0*x[i, t]*x[i, t] - z[i, t] + I1;
            dz = (1/tau0)*(4*(x[i, t] - eta[i]) - z[i, t] - K*gx);
            x[i, t+1] = x[i, t] + dt*dx + sqrt(dt)*x_eta[i, t]*sig; 
            z[i, t+1] = z[i, t] + dt*dz + sqrt(dt)*z_eta[i, t]*sig; 
                         }
                     }   

    xhat=amplitude*(to_vector(x) + offset);

    /* sampling*/

    target+=normal_lpdf(xs| xhat, eps);  

}    
          
generated quantities {
    matrix[nn, nt] x;
    matrix[nn, nt] z;
    vector[nn*nt] xhat_qqc;
    vector[nn*nt] x_ppc;
    real gx;
    real dx;
    real dz;

    for (i in 1:nn) {
      x[i, 1] = -1.5 + x_init[i];
      z[i, 1] = +3.5 + z_init[i];
    } 
    
    for (t in 1:(nt-1)) {
        for (i in 1:nn) {
            gx = 0;
            for (j in 1:nn)
                    gx = gx + SC[i, j]*(x[j, t] - x[i, t]);
            dx = 1.0 - x[i, t]*x[i, t]*x[i, t] - 2.0*x[i, t]*x[i, t] - z[i, t] + I1;
            dz = (1/tau0)*(4*(x[i, t] - eta[i]) - z[i, t] - K*gx);
            x[i, t+1] = x[i, t] + dt*dx + sqrt(dt)*x_eta[i, t]*sig; 
            z[i, t+1] = z[i, t] + dt*dz + sqrt(dt)*z_eta[i, t]*sig; 
                         }
                     }  


    xhat_qqc=amplitude*(to_vector(x) + offset);

    for (i in 1:(nn*nt)){
        x_ppc[i] = normal_rng(xhat_qqc[i], eps);
      }
}

