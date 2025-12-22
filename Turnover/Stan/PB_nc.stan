data{
  int n;
  vector[n] PB_plant;
  vector[n] PB_detritus;
  array[n] int Community;
  int n_Community;
}

transformed data{
  // log10 transform
  vector[n] log_PB_plant = log10(PB_plant);
  vector[n] log_PB_detritus = log10(PB_detritus);
}

parameters{
  // Global parameters
  real alpha_mu; // intercept
  real beta_mu; // slope
  real<lower=0> alpha_sigma;
  real<lower=0> beta_sigma;
  real<lower=0> sigma; // likelihood standard deviation
  
  // Community parameters
  vector[n_Community] alpha_z; // z-scores
  vector[n_Community] beta_z;
}

transformed parameters{
  // Convert z-scores
  vector[n_Community] alpha = alpha_z * alpha_sigma + alpha_mu;
  vector[n_Community] beta = beta_z * beta_sigma + beta_mu;
}

model{
  // Priors
  /// Global parameters
  alpha_mu ~ normal( 0 , 1 );
  beta_mu ~ normal( 0 , 1 );
  alpha_sigma ~ normal( 0 , 1 ) T[0,]; // half-normal priors
  beta_sigma ~ normal( 0 , 1 ) T[0,];
  sigma ~ exponential( 1 );
  
  /// Community parameters
  alpha_z ~ normal( 0 , 1 );
  beta_z ~ normal( 0 , 1 );

  // Model
  vector[n] mu = alpha[Community] + beta[Community] .* log_PB_plant;
  log_PB_detritus ~ normal( mu , sigma );
}