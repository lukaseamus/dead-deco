functions{
  // Beta prime log probability density function
  real betap_lpdf( real y , real alpha , real beta ) {
    return ( alpha - 1 ) * log( y )
    - ( alpha + beta ) * log1p( y ) -
    lbeta( alpha , beta );
  }
}

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
  // Parameters describing mean
  /// Global parameters
  real log_delta_mu; // delta = alpha + tau
  real log_mu_mu;
  real log_tau_mu;
  
  real<lower=0> log_delta_sigma_s;
  real<lower=0> log_mu_sigma_s;
  real<lower=0> log_tau_sigma_s;
  
  real<lower=0> log_delta_sigma_e;
  real<lower=0> log_mu_sigma_e;
  real<lower=0> log_tau_sigma_e;
  
  /// Species parameters
  vector[n_Species] log_delta_s; 
  vector[n_Species] log_mu_s;
  vector[n_Species] log_tau_s;
  
  /// Experiment parameters
  vector[n_Experiment] log_delta_e; 
  vector[n_Experiment] log_mu_e;
  vector[n_Experiment] log_tau_e;
  
  // Parameters describing precision
  real<lower=0> epsilon;
  real<lower=0> lambda;
  real<lower=0> theta;
}

model{
  // Priors
  /// Likelihood mean
  //// Global parameters
  log_delta_mu ~ normal( log(0.05) , 0.4 );
  log_mu_mu ~ normal( log(50) , 0.4 );
  log_tau_mu ~ normal( log(0.1) , 0.4 );
  
  log_delta_sigma_s ~ normal( 0 , 0.3 ) T[0,]; // half-normal priors
  log_mu_sigma_s ~ normal( 0 , 0.3 ) T[0,];
  log_tau_sigma_s ~ normal( 0 , 0.3 ) T[0,];
  
  log_delta_sigma_e ~ normal( 0 , 0.3 ) T[0,];
  log_mu_sigma_e ~ normal( 0 , 0.3 ) T[0,];
  log_tau_sigma_e ~ normal( 0 , 0.3 ) T[0,];
  
  //// Species parameters
  log_delta_s ~ normal( log_delta_mu , log_delta_sigma_s );
  log_mu_s ~ normal( log_mu_mu , log_mu_sigma_s );
  log_tau_s ~ normal( log_tau_mu , log_tau_sigma_s );
  
  //// Experiment parameters
  log_delta_e ~ normal( 0 , log_delta_sigma_e );
  log_mu_e ~ normal( 0 , log_mu_sigma_e );
  log_tau_e ~ normal( 0 , log_tau_sigma_e );
  
  /// Likelihood precision
  epsilon ~ gamma( square(4e4) / square(2e4) , 4e4 / square(2e4) );
  lambda ~ exponential( 1 );
  theta ~ gamma( square(500) / square(250) , 500 / square(250) );
  
  // Model
  /// Likelihood mean
  //// Parameters
  vector[n] delta = exp( log_delta_s[Species] + log_delta_e[Experiment] );
  vector[n] mu = exp( log_mu_s[Species] + log_mu_e[Experiment] );
  vector[n] tau = exp( log_tau_s[Species] + log_tau_e[Experiment] );
  vector[n] alpha = delta - tau;
  
  //// Function
  vector[n] r_mu = exp(
      Day .* alpha - ( alpha + tau ) .* 
      mu ./ 5 .* (
        log1p_exp( 5 ./ mu .* ( Day - mu ) ) -
        log1p_exp( -5 )
      )
    );
  
  /// Likelihood precision
  vector[n] nu = theta + exp(
    log( epsilon - theta ) - lambda .* Day
  );
  
  // Beta prime likelihood
  for ( i in 1:n ) { // loop because betap isn't vectorised
    Ratio[i] ~ betap( r_mu[i] * ( 1 + nu[i] ) , 2 + nu[i] );
  }
}