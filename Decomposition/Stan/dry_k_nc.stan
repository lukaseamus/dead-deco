data{
  int n;
  vector[n] Day;
  vector[n] Ratio_mean;
  array[n] int Species;
  int n_Species;
  array[n] int Treatment;
  int n_Treatment;
}

parameters{
  // Global parameters
  real log_k_mu;
  real log_sigma_mu;
  
  real<lower=0> log_k_sigma;
  real<lower=0> log_sigma_sigma;
  
  // Species/treatment parameters
  matrix[n_Species, n_Treatment] log_k_z; // z-scores
  matrix[n_Species, n_Treatment] log_sigma_z;
}

transformed parameters{
  matrix[n_Species, n_Treatment] log_k = log_k_z * log_k_sigma + log_k_mu;
  matrix[n_Species, n_Treatment] log_sigma = log_sigma_z * log_sigma_sigma + log_sigma_mu;
}

model{
  // Priors
  /// Global parameters
  log_k_mu ~ normal( log(0.06) , 0.6 );
  log_sigma_mu ~ normal( log(0.1) , 0.3 );
  
  log_k_sigma ~ normal( 0 , 0.6 ) T[0,]; // half-normal priors
  log_sigma_sigma ~ normal( 0 , 0.3 ) T[0,];
  
  /// Species/treatment parameters
  to_vector(log_k_z) ~ normal( 0 , 1 );
  to_vector(log_sigma_z) ~ normal( 0 , 1 );
  
  // Model
  /// Parameters
  vector[n] k;
  vector[n] sigma;

  for ( i in 1:n ) {
    k[i] = exp( log_k[ Species[i], Treatment[i] ] );
    sigma[i] = exp( log_sigma[ Species[i], Treatment[i] ] );
  }

  /// Function
  vector[n] r_mu = exp( -k .* Day );
  
  // Normal likelihood
  Ratio_mean ~ normal( r_mu , sigma );
}