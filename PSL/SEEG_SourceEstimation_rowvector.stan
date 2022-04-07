data {
  int nt;
  int nn;
  int ns;

  row_vector[ns] Seeg[nt]; 
  matrix[ns, nn] Gain;
}

transformed data {
}

parameters {
  row_vector[nn] X_est[nt];
  real<lower=0.0> eps; 
}

transformed parameters {
}

model {
  row_vector[ns] Seeg_model[nt];

  for(i in 1:nn){
      X_est[:,i] ~ normal(-2.,1.);
  }

  eps ~ normal(0.,1.);

  for(t in 1:nt){
   Seeg_model[t] = (Gain * (X_est[t,]'))' ;
  }

  for (t in 1:nt) {
    target += normal_lpdf(Seeg[t] | Seeg_model[t], eps);
  }

}


generated quantities {      
}
