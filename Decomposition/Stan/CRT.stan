data{
  int n;
  vector[n] CRT;
}

parameters{
  real mu;
  real<lower=0> sigma;
}

model{
  // Priors
  mu ~ normal( log(75) , 1 );
  sigma ~ exponential( 1 );

  // Lognormal likelihood
  CRT ~ lognormal( mu , sigma );
}