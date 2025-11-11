data{
  int n;
  vector<lower=0>[n] Fresh;
  vector<lower=0>[n] Dry;
  array[n] int Species;
  int n_Species;
}

parameters{
  // Hyperparameters
  real<lower=0, upper=1> beta_mu;
  real<lower=0> beta_nu;
  real<lower=0> cv_mu;
  real<lower=0> cv_theta;
  
  // Species parameters
  vector<lower=0, upper=1>[n_Species] beta;
  vector<lower=0>[n_Species] cv;
}

model{
  // Hyperpriors
  beta_mu ~ beta( 0.29 * 30 , (1 - 0.29) * 30 );
  beta_nu ~ gamma( square(20) / square(10) , 20 / square(10) );
  cv_mu ~ exponential( 1 );
  cv_theta ~ exponential( 1 );
  
  // Species priors
  beta ~ beta( beta_mu * beta_nu , (1 - beta_mu) * beta_nu );
  cv ~ gamma( cv_mu / cv_theta , 1 / cv_theta );
  
  // Model
  vector[n] mu = beta[Species] .* Fresh;
  vector[n] sigma = cv[Species] .* mu;
  
  // Truncated Gaussian likelihood
  Dry ~ normal( mu , sigma ) T[0,];
}