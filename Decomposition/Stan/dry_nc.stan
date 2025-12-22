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
  vector[n] Ratio_mean;
  vector[n] Ratio_sd;
  array[n] int Species;
  int n_Species;
  array[n] int Treatment;
  int n_Treatment;
}

transformed data{
  // Convert sd to nu because this is easier on the sampler
  vector[n] Ratio_nu = Ratio_mean .* ( 1 + Ratio_mean ) ./ Ratio_sd^2;
}

parameters{
  // Parameter describing true, unobserved ratio
  vector<lower=0>[n] r;
  
  // Parameters describing mean
  /// Global parameters
  real log_delta_mu; // delta = alpha + tau
  real log_mu_mu;
  real log_tau_mu;
  
  real<lower=0> log_delta_sigma;
  real<lower=0> log_mu_sigma;
  real<lower=0> log_tau_sigma;
  
  /// Species/treatment parameters
  matrix[n_Species, n_Treatment] log_delta_z; // z-scores
  matrix[n_Species, n_Treatment] log_mu_z;
  matrix[n_Species, n_Treatment] log_tau_z;
  
  // Parameters describing precision
  /// Global parameters
  real log_epsilon_mu;
  real log_lambda_mu;
  real log_theta_mu;

  real<lower=0> log_epsilon_sigma;
  real<lower=0> log_lambda_sigma;
  real<lower=0> log_theta_sigma;
  
  /// Species/treatment parameters
  matrix[n_Species, n_Treatment] log_epsilon_z;
  matrix[n_Species, n_Treatment] log_lambda_z;
  matrix[n_Species, n_Treatment] log_theta_z;
}

transformed parameters{
  matrix[n_Species, n_Treatment] log_delta = log_delta_z * log_delta_sigma + log_delta_mu;
  matrix[n_Species, n_Treatment] log_mu = log_mu_z * log_mu_sigma + log_mu_mu;
  matrix[n_Species, n_Treatment] log_tau = log_tau_z * log_tau_sigma + log_tau_mu;
  
  matrix[n_Species, n_Treatment] log_epsilon = log_epsilon_z * log_epsilon_sigma + log_epsilon_mu;
  matrix[n_Species, n_Treatment] log_lambda = log_lambda_z * log_lambda_sigma + log_lambda_mu;
  matrix[n_Species, n_Treatment] log_theta = log_theta_z * log_theta_sigma + log_theta_mu;
}

model{
  // Priors
  /// Likelihood mean
  //// Global parameters
  log_delta_mu ~ normal( log(0.03) , 0.3 );
  log_mu_mu ~ normal( log(40) , 0.3 );
  log_tau_mu ~ normal( log(0.06) , 0.3 );
  
  log_delta_sigma ~ normal( 0 , 0.3 ) T[0,]; // half-normal priors
  log_mu_sigma ~ normal( 0 , 0.3 ) T[0,];
  log_tau_sigma ~ normal( 0 , 0.3 ) T[0,];
  
  //// Species/treatment parameters
  to_vector(log_delta_z) ~ normal( 0 , 1 );
  to_vector(log_mu_z) ~ normal( 0 , 1 );
  to_vector(log_tau_z) ~ normal( 0 , 1 );
  
  /// Likelihood precision
  //// Global parameters
  log_epsilon_mu ~ normal( log(4e4) , 0.3 );
  log_lambda_mu ~ normal( log(0.1) , 0.3 );
  log_theta_mu ~ normal( log(500) , 0.3 );
  
  log_epsilon_sigma ~ normal( 0 , 0.3 ) T[0,];
  log_lambda_sigma ~ normal( 0 , 0.3 ) T[0,];
  log_theta_sigma ~ normal( 0 , 0.3 ) T[0,];
  
  //// Species/treatment parameters
  to_vector(log_epsilon_z) ~ normal( 0 , 1 );
  to_vector(log_lambda_z) ~ normal( 0 , 1 );
  to_vector(log_theta_z) ~ normal( 0 , 1 );
  
  // Model
  /// Likelihood mean
  //// Parameters
  vector[n] delta;
  vector[n] mu;
  vector[n] tau;
  vector[n] alpha;

  for ( i in 1:n ) {
    delta[i] = exp( log_delta[ Species[i], Treatment[i] ] );
    mu[i] = exp( log_mu[ Species[i], Treatment[i] ] );
    tau[i] = exp( log_tau[ Species[i], Treatment[i] ] );
    alpha[i] = delta[i] - tau[i];
  }

  //// Function
  vector[n] r_mu = exp(
      Day .* alpha - ( alpha + tau ) .* 
      mu ./ 5 .* (
        log1p_exp( 5 ./ mu .* ( Day - mu ) ) -
        log1p_exp( -5 )
      )
    );
  
  /// Likelihood precision
  //// Parameters
  vector[n] epsilon;
  vector[n] lambda;
  vector[n] theta;
  
  for ( i in 1:n ) {
    epsilon[i] = exp( log_epsilon[ Species[i], Treatment[i] ] );
    lambda[i] = exp( log_lambda[ Species[i], Treatment[i] ] );
    theta[i] = exp( log_theta[ Species[i], Treatment[i] ] );
  }
  
  //// Function
  vector[n] nu = theta + exp(
      log( epsilon - theta ) - lambda .* Day
    );
  
  // Beta prime likelihood
  for ( i in 1:n ) { // loop because betap isn't vectorised
    r[i] ~ betap( r_mu[i] * ( 1 + nu[i] ) , 2 + nu[i] );
  }
  
  // Beta prime measurement error model
  for ( i in 1:n ) {
   Ratio_mean[i] ~ betap(
      r[i] * ( 1 + Ratio_nu[i] ),
      2 + Ratio_nu[i]
    );
  }
}