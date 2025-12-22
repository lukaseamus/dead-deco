data{
  int n;
  vector[n] Turnover;
  array[n] int Reference;
  int n_Reference;
}

parameters{
  // Global parameters
  real alpha_mu; // intercept on log scale
  real<lower=0> alpha_sigma;
  real<lower=0> sigma; // likelihood standard deviation (log scale)
  
  // Group parameters
  vector[n_Reference] alpha_z; // z-scores
}

transformed parameters{
  // Convert z-scores
  vector[n_Reference] alpha = alpha_z * alpha_sigma + alpha_mu;
}

model{
  // Priors
  /// Global parameters
  alpha_mu ~ normal( log(7.54) , 1 );
  alpha_sigma ~ normal( 0 , 0.6 ) T[0,]; // half-normal prior
  sigma ~ exponential( 1 );
  
  /// Group parameters
  alpha_z ~ normal( 0 , 1 );

  // Model
  vector[n] mu = alpha[Reference];
  Turnover ~ lognormal( mu , sigma );
}