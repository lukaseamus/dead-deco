data{
  int n;
  vector[n] Day;
  vector[n] Proportion_mean;
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
  matrix[n_Species, n_Treatment] log_k;
  matrix[n_Species, n_Treatment] log_sigma;
}

model{
  // Priors
  /// Global parameters
  log_k_mu ~ normal( log(0.06) , 0.5 );
  log_sigma_mu ~ normal( log(0.1) , 0.4 );
  
  log_k_sigma ~ normal( 0 , 0.6 ) T[0,]; // half-normal priors
  log_sigma_sigma ~ normal( 0 , 0.5 ) T[0,];
  
  /// Species/treatment parameters
  to_vector(log_k) ~ normal( log_k_mu , log_k_sigma );
  to_vector(log_sigma) ~ normal( log_sigma_mu , log_sigma_sigma );

  // Model
  /// Parameters
  vector[n] k;
  vector[n] sigma;

  for ( i in 1:n ) {
    k[i] = exp( log_k[ Species[i], Treatment[i] ] );
    sigma[i] = exp( log_sigma[ Species[i], Treatment[i] ] );
  }
  
  /// Function
  vector[n] p_mu = exp( -k .* Day );

  // Normal likelihood
  Proportion_mean ~ normal( p_mu , sigma );
}