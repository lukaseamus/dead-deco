data{
  int n;
  vector[n] Proportion;
}

parameters{
  real<lower=0, upper=1> mu;
  real<lower=0> nu;
}

model{
  // Priors
  mu ~ beta( 0.71 * 20 , (1 - 0.71) * 20 );
  nu ~ gamma( square(20) / square(10) , 20 / square(10) );

  // Beta likelihood
  Proportion ~ beta( mu * nu , (1 - mu) * nu );
}