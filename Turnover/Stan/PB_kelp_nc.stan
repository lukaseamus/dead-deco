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
  vector[n_Species] alpha_z_s; // z-scores
  vector[n_Level] alpha_z_l;
  vector[n_Reference] alpha_z_r;
}

transformed parameters{
  // Convert z-scores
  vector[n_Species] alpha_s = alpha_z_s * alpha_sigma_s + alpha_mu;
  vector[n_Level] alpha_l = alpha_z_l * alpha_sigma_l + 0;
  vector[n_Reference] alpha_r = alpha_z_r * alpha_sigma_r + 0;
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
  alpha_z_s ~ normal( 0 , 1 );
  alpha_z_l ~ normal( 0 , 1 );
  alpha_z_r ~ normal( 0 , 1 );

  // Model
  vector[n] mu = alpha_s[Species] + alpha_l[Level] + alpha_r[Reference];
  Turnover ~ lognormal( mu , sigma );
}