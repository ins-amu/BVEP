data {
  int nt;
  int nn;
  int ns;

  matrix[nt, ns] Seeg;
  matrix[ns, nn] Gain;
}

transformed data {
  vector[nt*ns] Seeg_obs_vect;
  Seeg_obs_vect=to_vector(Seeg);
}

parameters {
  matrix[nn, nt] X_est;
  real<lower=0.0> eps; 
}

transformed parameters {
}

model {
  matrix[nt, ns] Seeg_model;
  vector[nt*ns]  Seeg_model_vect;

  eps ~ normal(0.,1.);

  for(t in 1:nt){
    Seeg_model[t,]=to_row_vector((Gain*X_est[,t]));
  }

  Seeg_model_vect=to_vector(Seeg_model);

  for(i in 1:(nt*ns)){
    target+=  normal_lpdf(Seeg_obs_vect[i]| Seeg_model_vect[i], eps);  
  }

}


generated quantities {      
}
