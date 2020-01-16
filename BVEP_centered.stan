/* Written in INS Marseille*/
data {
    int nn; \\number of brain regions
    int nt;  \\number of data points per time series
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
    matrix[nn, nt] x;
    matrix[nn, nt] z;
    vector[nn] eta;
    real <lower=0.0> K;
    real <lower=0.0> amplitude;
    real offset; 
    real<lower=0.0> eps;
    real<lower=0.0> sig;    
}

transformed parameters {
}

model {
    vector[nn*nt] xhat;
    real dx;
    real dz;
    real gx;

    /* priors*/

    eta ~ normal(-2.5, 1.0);       
    amplitude ~ normal(1.0 , 1.0);
    offset ~ normal(0.0, 1.0);
    K ~ normal(Ks, 1.0);  
    eps ~ normal(0.0, 1.0);   
    sig ~ normal(0.0, 1.0);   

    /* integrate & predict */
    
    for (i in 1:nn) {
     x[i, 1] ~ normal(initial_val[1], 1.0); 
     z[i, 1] ~ normal(initial_val[2], 1.0);
   } 
    
    for (t in 1:(nt-1)) {
        for (i in 1:nn) {
            gx = 0;
            for (j in 1:nn)
                    gx = gx + SC[i, j]*(x[j, t] - x[i, t]);
            dx = 1.0 - x[i, t]*x[i, t]*x[i, t] - 2.0*x[i, t]*x[i, t] - z[i, t] + I1;
            dz = (1/tau0)*(4*(x[i, t] - eta[i]) - z[i, t] - K*gx);
            x[i, t+1] ~ normal(x[i, t] + dt*dx, sqrt(dt)*sig); 
            z[i, t+1] ~ normal(z[i, t] + dt*dz, sqrt(dt)*sig); 
                         }
                     }
       
    xhat=amplitude*(to_vector(x) + offset);

    /* sampling*/

    target+=normal_lpdf(xs| xhat, eps); 

}

generated quantities {
    vector[nn*nt] xhat_qqc;
    vector[nn*nt] x_ppc;

    xhat_qqc=amplitude*(to_vector(x) + offset);

    for (i in 1:(nn*nt)){
        x_ppc[i] = normal_rng(xhat_qqc[i], eps);
      } 
}
