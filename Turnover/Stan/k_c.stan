data{
  int n;
  vector[n] k;
  array[n] int Species;
  int n_Species;
  array[n] int Review;
  int n_Review;
  array[n] int Reference;
  int n_Reference;
}

parameters{
  // Global parameters
  real alpha_mu;
  real<lower=0> alpha_sigma_s;
  real<lower=0> alpha_sigma_m;
  real<lower=0> alpha_sigma_r;
  real<lower=0> sigma; // likelihood standard deviation
  
  // Group parameters
  vector[n_Species] alpha_s;
  vector[n_Review] alpha_m;
  vector[n_Reference] alpha_r;
}

model{
  // Priors
  /// Global parameters
  alpha_mu ~ normal( 0.01 , 0.02 );
  alpha_sigma_s ~ normal( 0 , 0.02 ) T[0,]; // half-normal priors
  alpha_sigma_m ~ normal( 0 , 0.02 ) T[0,];
  alpha_sigma_r ~ normal( 0 , 0.02 ) T[0,];
  sigma ~ exponential( 1 );
  
  /// Group parameters
  alpha_s ~ normal( alpha_mu , alpha_sigma_s );
  alpha_m ~ normal( 0 , alpha_sigma_m );
  alpha_r ~ normal( 0 , alpha_sigma_r );

  // Model
  vector[n] mu = alpha_s[Species] + alpha_m[Review] + alpha_r[Reference];
  k ~ normal( mu , sigma );
}