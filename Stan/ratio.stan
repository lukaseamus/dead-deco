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
  real<lower=0> sigma_mu;
  real<lower=0> sigma_theta;
  
  // Species parameters
  vector<lower=0, upper=1>[n_Species] beta;
  vector<lower=0>[n_Species] sigma;
}

model{
  // Hyperpriors
  beta_mu ~ beta( 0.29 * 30 , (1 - 0.29) * 30 );
  beta_nu ~ gamma( square(30) / square(20) , 30 / square(20) );
  sigma_mu ~ exponential( 1 );
  sigma_theta ~ exponential( 1 );
  
  // Species priors
  beta ~ beta( beta_mu * beta_nu , (1 - beta_mu) * beta_nu );
  sigma ~ gamma( sigma_mu / sigma_theta , 1 / sigma_theta );
  
  // Model
  vector[n] mu = beta[Species] .* Fresh;
  
  // Truncated Gaussian likelihood
  Dry ~ normal( mu , sigma[Species] ) T[0,];
}