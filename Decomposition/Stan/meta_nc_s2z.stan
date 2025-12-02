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
  array[n] int Condition;
  int n_Condition;
}

parameters{
  // Parameters describing mean
  /// Condition parameters
  vector[n_Condition] alpha_c;
  vector[n_Condition] log_mu_c;
  vector[n_Condition] log_tau_c;
  
  vector<lower=0>[n_Condition] alpha_sigma_s;
  vector<lower=0>[n_Condition] log_mu_sigma_s;
  vector<lower=0>[n_Condition] log_tau_sigma_s;
  
  vector<lower=0>[n_Condition] alpha_sigma_e;
  vector<lower=0>[n_Condition] log_mu_sigma_e;
  vector<lower=0>[n_Condition] log_tau_sigma_e;
  
  /// Species parameters
  array[n_Condition] sum_to_zero_vector[n_Species] alpha_z_s; // z-scores
  array[n_Condition] sum_to_zero_vector[n_Species] log_mu_z_s;
  array[n_Condition] sum_to_zero_vector[n_Species] log_tau_z_s;
  
  /// Experiment parameters
  array[n_Condition] sum_to_zero_vector[n_Experiment] alpha_z_e;
  array[n_Condition] sum_to_zero_vector[n_Experiment] log_mu_z_e;
  array[n_Condition] sum_to_zero_vector[n_Experiment] log_tau_z_e;
  
  // Parameters describing precision
  /// Condition parameters
  vector[n_Condition] log_epsilon_c;
  vector[n_Condition] log_lambda_c;
  vector[n_Condition] log_theta_c;

  vector<lower=0>[n_Condition] log_epsilon_sigma_s;
  vector<lower=0>[n_Condition] log_lambda_sigma_s;
  vector<lower=0>[n_Condition] log_theta_sigma_s;
  
  vector<lower=0>[n_Condition] log_epsilon_sigma_e;
  vector<lower=0>[n_Condition] log_lambda_sigma_e;
  vector<lower=0>[n_Condition] log_theta_sigma_e;
  
  /// Species parameters
  array[n_Condition] sum_to_zero_vector[n_Species] log_epsilon_z_s; // z-scores
  array[n_Condition] sum_to_zero_vector[n_Species] log_lambda_z_s;
  array[n_Condition] sum_to_zero_vector[n_Species] log_theta_z_s;
  
  /// Experiment parameters
  array[n_Condition] sum_to_zero_vector[n_Experiment] log_epsilon_z_e;
  array[n_Condition] sum_to_zero_vector[n_Experiment] log_lambda_z_e;
  array[n_Condition] sum_to_zero_vector[n_Experiment] log_theta_z_e;
}

model{
  // Priors
  /// Likelihood mean
  //// Condition parameters
  alpha_c ~ normal( 0 , 0.005 );
  log_mu_c ~ normal( log(100) , 0.1 );
  log_tau_c ~ normal( log(0.1) , 0.1 );
  
  alpha_sigma_s ~ normal( 0 , 0.005 ) T[0,]; // half-normal priors
  log_mu_sigma_s ~ normal( 0 , 0.1 ) T[0,];
  log_tau_sigma_s ~ normal( 0 , 0.1 ) T[0,];
  
  alpha_sigma_e ~ normal( 0 , 0.01 ) T[0,];
  log_mu_sigma_e ~ normal( 0 , 0.1 ) T[0,];
  log_tau_sigma_e ~ normal( 0 , 0.1 ) T[0,];
  
  //// Species and experiment parameters
  for (i in 1:n_Condition) {
    alpha_z_s[i][] ~ normal( 0 , 1 );
    log_mu_z_s[i][] ~ normal( 0 , 1 );
    log_tau_z_s[i][] ~ normal( 0 , 1 );
    alpha_z_e[i][] ~ normal( 0 , 1 );
    log_mu_z_e[i][] ~ normal( 0 , 1 );
    log_tau_z_e[i][] ~ normal( 0 , 1 );
  }
  
  //// Convert z-scores
  matrix[n_Condition, n_Species] alpha_s;
  matrix[n_Condition, n_Species] log_mu_s;
  matrix[n_Condition, n_Species] log_tau_s;
  matrix[n_Condition, n_Experiment] alpha_e;
  matrix[n_Condition, n_Experiment] log_mu_e;
  matrix[n_Condition, n_Experiment] log_tau_e;
      
  for (i in 1:n_Condition) {
    for (j in 1:n_Species) {
      alpha_s[i, j] = alpha_z_s[i][j] * sqrt( n_Species * inv( n_Species - 1 ) ) * 
                      alpha_sigma_s[i] + 0;
      log_mu_s[i, j] = log_mu_z_s[i][j] * sqrt( n_Species * inv( n_Species - 1 ) ) * 
                       log_mu_sigma_s[i] + 0;
      log_tau_s[i, j] = log_tau_z_s[i][j] * sqrt( n_Species * inv( n_Species - 1 ) ) * 
                        log_tau_sigma_s[i] + 0;
    }
    for (j in 1:n_Experiment) {
      alpha_e[i, j] = alpha_z_e[i][j] * sqrt( n_Experiment * inv( n_Experiment - 1 ) ) * 
                      alpha_sigma_e[i] + 0;
      log_mu_e[i, j] = log_mu_z_e[i][j] * sqrt( n_Experiment * inv( n_Experiment - 1 ) ) * 
                       log_mu_sigma_e[i] + 0;
      log_tau_e[i, j] = log_tau_z_e[i][j] * sqrt( n_Experiment * inv( n_Experiment - 1 ) ) *
                        log_tau_sigma_e[i] + 0;
    }
  }
  
  /// Likelihood precision
  //// Condition parameters
  log_epsilon_c ~ normal( log(4e4) , 0.1 );
  log_lambda_c ~ normal( log(0.1) , 0.1 );
  log_theta_c ~ normal( log(500) , 0.1 );
  
  log_epsilon_sigma_s ~ normal( 0 , 0.1 ) T[0,];
  log_lambda_sigma_s ~ normal( 0 , 0.1 ) T[0,];
  log_theta_sigma_s ~ normal( 0 , 0.1 ) T[0,];
  
  log_epsilon_sigma_e ~ normal( 0 , 0.1 ) T[0,];
  log_lambda_sigma_e ~ normal( 0 , 0.1 ) T[0,];
  log_theta_sigma_e ~ normal( 0 , 0.1 ) T[0,];
  
  //// Species and experiment parameters
  for (i in 1:n_Condition) {
    log_epsilon_z_s[i][] ~ normal( 0 , 1 );
    log_lambda_z_s[i][] ~ normal( 0 , 1 );
    log_theta_z_s[i][] ~ normal( 0 , 1 );
    log_epsilon_z_e[i][] ~ normal( 0 , 1 );
    log_lambda_z_e[i][] ~ normal( 0 , 1 );
    log_theta_z_e[i][] ~ normal( 0 , 1 );
  }
  
  //// Convert z-scores
  matrix[n_Condition, n_Species] log_epsilon_s;
  matrix[n_Condition, n_Species] log_lambda_s;
  matrix[n_Condition, n_Species] log_theta_s;
  matrix[n_Condition, n_Experiment] log_epsilon_e;
  matrix[n_Condition, n_Experiment] log_lambda_e;
  matrix[n_Condition, n_Experiment] log_theta_e;
  
  for (i in 1:n_Condition) {
    for (j in 1:n_Species) {
      log_epsilon_s[i, j] = log_epsilon_z_s[i][j] * sqrt( n_Species * inv( n_Species - 1 ) ) * 
                            log_epsilon_sigma_s[i] + 0;
      log_lambda_s[i, j] = log_lambda_z_s[i][j] * sqrt( n_Species * inv( n_Species - 1 ) ) * 
                           log_lambda_sigma_s[i] + 0;
      log_theta_s[i, j] = log_theta_z_s[i][j] * sqrt( n_Species * inv( n_Species - 1 ) ) * 
                          log_theta_sigma_s[i] + 0;
    }
    for (j in 1:n_Experiment) {
      log_epsilon_e[i, j] = log_epsilon_z_e[i][j] * sqrt( n_Experiment * inv( n_Experiment - 1 ) ) * 
                            log_epsilon_sigma_e[i] + 0;
      log_lambda_e[i, j] = log_lambda_z_e[i][j] * sqrt( n_Experiment * inv( n_Experiment - 1 ) ) * 
                           log_lambda_sigma_e[i] + 0;
      log_theta_e[i, j] = log_theta_z_e[i][j] * sqrt( n_Experiment * inv( n_Experiment - 1 ) ) *
                          log_theta_sigma_e[i] + 0;
    }
  }
  
  // Model
  /// Likelihood mean
  //// Parameters
  vector[n] alpha;
  vector[n] mu;
  vector[n] tau;

  for ( i in 1:n ) {
    alpha[i] = 
      alpha_c[ Condition[i] ] + 
      alpha_s[ Condition[i] , Species[i] ] + 
      alpha_e[ Condition[i] , Experiment[i] ];
    mu[i] = exp(
      log_mu_c[ Condition[i] ] + 
      log_mu_s[ Condition[i] , Species[i] ] + 
      log_mu_e[ Condition[i] , Experiment[i] ]
    );
    tau[i] = exp(
      log_tau_c[ Condition[i] ] + 
      log_tau_s[ Condition[i] , Species[i] ] + 
      log_tau_e[ Condition[i] , Experiment[i] ]
    );
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
    epsilon[i] = exp(
      log_epsilon_c[ Condition[i] ] + 
      log_epsilon_s[ Condition[i] , Species[i] ] + 
      log_epsilon_e[ Condition[i] , Experiment[i] ]
    );
    lambda[i] = exp(
      log_lambda_c[ Condition[i] ] + 
      log_lambda_s[ Condition[i] , Species[i] ] + 
      log_lambda_e[ Condition[i] , Experiment[i] ]
    );
    theta[i] = exp(
      log_theta_c[ Condition[i] ] + 
      log_theta_s[ Condition[i] , Species[i] ] + 
      log_theta_e[ Condition[i] , Experiment[i] ]
    );
  }
  
  //// Function (this structure is better for most multilevel models)
  vector[n] nu = theta + ( epsilon - theta ) .* exp( -lambda .* Day );
  
  // Beta prime likelihood
  for ( i in 1:n ) { // loop because betap isn't vectorised
    Ratio[i] ~ betap( r_mu[i] * ( 1 + nu[i] ) , 2 + nu[i] );
  }
}

generated quantities{
  // Save converted z-scores
  matrix[n_Condition, n_Species] alpha_s;
  matrix[n_Condition, n_Species] log_mu_s;
  matrix[n_Condition, n_Species] log_tau_s;
  matrix[n_Condition, n_Experiment] alpha_e;
  matrix[n_Condition, n_Experiment] log_mu_e;
  matrix[n_Condition, n_Experiment] log_tau_e;
      
  for (i in 1:n_Condition) {
    for (j in 1:n_Species) {
      alpha_s[i, j] = alpha_z_s[i][j] * sqrt( n_Species * inv( n_Species - 1 ) ) * 
                      alpha_sigma_s[i] + 0;
      log_mu_s[i, j] = log_mu_z_s[i][j] * sqrt( n_Species * inv( n_Species - 1 ) ) * 
                       log_mu_sigma_s[i] + 0;
      log_tau_s[i, j] = log_tau_z_s[i][j] * sqrt( n_Species * inv( n_Species - 1 ) ) * 
                        log_tau_sigma_s[i] + 0;
    }
    for (j in 1:n_Experiment) {
      alpha_e[i, j] = alpha_z_e[i][j] * sqrt( n_Experiment * inv( n_Experiment - 1 ) ) * 
                      alpha_sigma_e[i] + 0;
      log_mu_e[i, j] = log_mu_z_e[i][j] * sqrt( n_Experiment * inv( n_Experiment - 1 ) ) * 
                       log_mu_sigma_e[i] + 0;
      log_tau_e[i, j] = log_tau_z_e[i][j] * sqrt( n_Experiment * inv( n_Experiment - 1 ) ) *
                        log_tau_sigma_e[i] + 0;
    }
  }
  
  matrix[n_Condition, n_Species] log_epsilon_s;
  matrix[n_Condition, n_Species] log_lambda_s;
  matrix[n_Condition, n_Species] log_theta_s;
  matrix[n_Condition, n_Experiment] log_epsilon_e;
  matrix[n_Condition, n_Experiment] log_lambda_e;
  matrix[n_Condition, n_Experiment] log_theta_e;
  
  for (i in 1:n_Condition) {
    for (j in 1:n_Species) {
      log_epsilon_s[i, j] = log_epsilon_z_s[i][j] * sqrt( n_Species * inv( n_Species - 1 ) ) * 
                            log_epsilon_sigma_s[i] + 0;
      log_lambda_s[i, j] = log_lambda_z_s[i][j] * sqrt( n_Species * inv( n_Species - 1 ) ) * 
                           log_lambda_sigma_s[i] + 0;
      log_theta_s[i, j] = log_theta_z_s[i][j] * sqrt( n_Species * inv( n_Species - 1 ) ) * 
                          log_theta_sigma_s[i] + 0;
    }
    for (j in 1:n_Experiment) {
      log_epsilon_e[i, j] = log_epsilon_z_e[i][j] * sqrt( n_Experiment * inv( n_Experiment - 1 ) ) * 
                            log_epsilon_sigma_e[i] + 0;
      log_lambda_e[i, j] = log_lambda_z_e[i][j] * sqrt( n_Experiment * inv( n_Experiment - 1 ) ) * 
                           log_lambda_sigma_e[i] + 0;
      log_theta_e[i, j] = log_theta_z_e[i][j] * sqrt( n_Experiment * inv( n_Experiment - 1 ) ) *
                          log_theta_sigma_e[i] + 0;
    }
  }
}