data{
  int n;
  vector[n] Day;
  vector[n] Ratio;
  array[n] int Species;
  int n_Species;
  array[n] int Experiment;
  int n_Experiment;
}

parameters{
  // Global parameters
  real log_k_mu;
  real<lower=0> sigma;
  
  real<lower=0> log_k_sigma_s;
  real<lower=0> log_k_sigma_e;
  
  // Species and experiment parameters
  vector[n_Species] log_k_s;
  vector[n_Experiment] log_k_e;
}

model{
  // Priors
  /// Global parameters
  log_k_mu ~ normal( log(0.1) , 0.6 );
  sigma ~ exponential( 1 );
  
  log_k_sigma_s ~ normal( 0 , 0.6 ) T[0,]; // half-normal priors
  log_k_sigma_e ~ normal( 0 , 0.6 ) T[0,];
  
  /// Species and experiment parameters
  log_k_s ~ normal( log_k_mu , log_k_sigma_s );
  log_k_e ~ normal( 0 , log_k_sigma_e );

  // Model
  /// Parameters
  vector[n] k = exp( log_k_s[Species] + log_k_e[Experiment] );
  
  /// Function
  vector[n] r_mu = exp( -k .* Day );

  // Normal likelihood
  Ratio ~ normal( r_mu , sigma );
}