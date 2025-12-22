data{
  int n;
  vector[n] Turnover;
  array[n] int Species;
  int n_Species;
  array[n] int Level;
  int n_Level;
  array[n] int Reference;
  int n_Reference;
}

parameters{
  // Global parameters
  real alpha_mu; // intercept on log scale
  real<lower=0> alpha_sigma_s;
  real<lower=0> alpha_sigma_l;
  real<lower=0> alpha_sigma_r;
  real<lower=0> sigma; // likelihood standard deviation (log scale)
  
  // Group parameters
  vector[n_Species] alpha_s;
  vector[n_Level] alpha_l;
  vector[n_Reference] alpha_r;
}

model{
  // Priors
  /// Global parameters
  alpha_mu ~ normal( log(7.54) , 1 );
  alpha_sigma_s ~ normal( 0 , 0.3 ) T[0,]; // half-normal priors
  alpha_sigma_l ~ normal( 0 , 0.3 ) T[0,];
  alpha_sigma_r ~ normal( 0 , 0.3 ) T[0,];
  sigma ~ exponential( 1 );
  
  /// Group parameters
  alpha_s ~ normal( alpha_mu , alpha_sigma_s );
  alpha_l ~ normal( 0 , alpha_sigma_l );
  alpha_r ~ normal( 0 , alpha_sigma_r );

  // Model
  vector[n] mu = alpha_s[Species] + alpha_l[Level] + alpha_r[Reference];
  Turnover ~ lognormal( mu , sigma );
}