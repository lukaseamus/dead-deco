#### Dead deco: kelp decomposition meta-analysis ####
#### Luka Seamus Wright                          ####

# 1. Prepare data ####
# 1.1 Load data and simulate observations ####
require(tidyverse)
require(magrittr)
require(here)
require(extraDistr)
set.seed(100)
deco <- here("Decomposition", "Decomposition_Meta.csv") %>% 
  read_csv(col_types = list("f", "c", "f", "f", "f", "f", "f")) %>%
  # Replace 0 with small constant within measurement error
  # because the model is undefined for y = 0
  mutate(Ratio_mean = if_else(Ratio_mean == 0, 1e-5, Ratio_mean)) %>%
  rowwise() %>%
  mutate(Ratio = if( !is.na(Ratio_sd) & !is.na(Ratio_n) ) {
    list(
      rbetapr(
        Ratio_n ,
        Ratio_mean * ( 1 + Ratio_mean * (1 + Ratio_mean) / Ratio_sd^2 ) , 
        2 + Ratio_mean * (1 + Ratio_mean) / Ratio_sd^2
      )
    )
  } else {
    list( Ratio_mean )
  }) %>%
  unnest(Ratio) %T>%
  print()

# 1.2 Check data ####
deco %>% nrow() # 2305 observations

deco %>%
  distinct(Reference, Location, Species, Treatment, 
           Condition) %>%
  nrow() # 186 experiments

deco %>%
  distinct(Reference, Location, Species, Treatment, 
           Condition, Experiment) %>%
  print(n = 186)

deco %>%
  distinct(Reference, Location, Species, Treatment, Condition, Experiment) %>%
  count(Experiment) %>%
  filter(n > 1)
# No duplicated experiment numbers 
           
deco %$% nlevels(Species) # 14 species
deco %>% 
  distinct(Species) %>% 
  mutate(Species = Species %>% as.character()) %>%
  arrange(Species)

deco %$% nlevels(Reference) # 25 references
deco %>% 
  distinct(Reference) %>%
  mutate(Reference = Reference %>% as.character()) %>%
  arrange(Reference) %>%
  print(n = 25)

# All seems to be in order

# 2. Model data ####
# 2.1 Optimal model ####
# 2.1.1 Visualise data ####
deco %>%
  ggplot() +
    geom_point(aes(Day, Ratio, colour = Condition),
               shape = 16, alpha = 0.5) +
    facet_wrap(~ Species) +
    theme_minimal()

require(ggh4x)
deco %>%
  ggplot() +
    geom_point(aes(Day, Ratio, colour = Condition),
               shape = 16, alpha = 0.5) +
    facet_nested_wrap(~ Species + Experiment,
                      nest_line = TRUE) +
    theme_minimal()

# The best strategy seems to be to use partial pooling for
# species and experiments, but within condition. This is
# a hierarchical model that does not share variance between
# conditions.

# 2.1.2 Prior simulation ####
# Here I use a model I developed specifically for macroalgal 
# decomposition (github.com/lukaseamus/limbodeco).

# R doesn't have a built-in log1p_exp function
log1p_exp <- function(x) {
  ifelse(
    x > 0, 
    x + log1p(exp(-x)),
    log1p(exp(x))
  )
}

# Simulate based on known averages (see Decomposition.R in
# this R project and Examples.R in github.com/lukaseamus/limbodeco)
tibble(n = 1:1e3,
       alpha_mu = rnorm( 1e3 , 0 , 0.005 ),
       log_mu_mu = rnorm( 1e3 , log(100) , 0.1 ),
       log_tau_mu = rnorm( 1e3 , log(0.1) , 0.1 ),
       alpha_sigma_s = rtnorm( 1e3 , 0 , 0.005 , 0 ), # half-normal priors
       log_mu_sigma_s = rtnorm( 1e3 , 0 , 0.1 , 0 ),
       log_tau_sigma_s = rtnorm( 1e3 , 0 , 0.1 , 0 ),
       alpha_sigma_e = rtnorm( 1e3 , 0 , 0.005 , 0 ),
       log_mu_sigma_e = rtnorm( 1e3 , 0 , 0.1 , 0 ),
       log_tau_sigma_e = rtnorm( 1e3 , 0 , 0.1 , 0 ),
       log_epsilon_mu = rnorm( 1e3 , log(4e4) , 0.1 ),
       log_lambda_mu = rnorm( 1e3 , log(0.1) , 0.1 ),
       log_theta_mu = rnorm( 1e3 , log(500) , 0.1 ),
       log_epsilon_sigma_s = rtnorm( 1e3 , 0 , 0.1 , 0 ),
       log_lambda_sigma_s = rtnorm( 1e3 , 0 , 0.1 , 0 ),
       log_theta_sigma_s = rtnorm( 1e3 , 0 , 0.1 , 0 ),
       log_epsilon_sigma_e = rtnorm( 1e3 , 0 , 0.1 , 0 ),
       log_lambda_sigma_e = rtnorm( 1e3 , 0 , 0.1 , 0 ),
       log_theta_sigma_e = rtnorm( 1e3 , 0 , 0.1 , 0 ),
       alpha = 
         rnorm( 1e3 , alpha_mu , alpha_sigma_s ) + 
         rnorm( 1e3 , 0 , alpha_sigma_e ),
       mu = exp(
         rnorm( 1e3 , log_mu_mu , log_mu_sigma_s ) + 
         rnorm( 1e3 , 0 , log_mu_sigma_e )
       ),
       tau = exp(
         rnorm( 1e3 , log_tau_mu , log_tau_sigma_s ) + 
         rnorm( 1e3 , 0 , log_tau_sigma_e )
       ),
       epsilon = exp(
         rnorm( 1e3 , log_epsilon_mu , log_epsilon_sigma_s ) + 
         rnorm( 1e3 , 0 , log_epsilon_sigma_e )
       ),
       lambda = exp(
         rnorm( 1e3 , log_lambda_mu , log_lambda_sigma_s ) + 
         rnorm( 1e3 , 0 , log_lambda_sigma_e )
       ),
       theta = exp(
         rnorm( 1e3 , log_theta_mu , log_theta_sigma_s ) + 
         rnorm( 1e3 , 0 , log_theta_sigma_e )
       )) %>%
  expand_grid(Day = deco %$% seq(min(Day), max(Day), length.out = 100)) %>%
  mutate(
    r_mu = exp(
      Day * alpha - ( alpha + tau ) * mu / 5 * (
        log1p_exp( 5 / mu * ( Day - mu ) ) - log1p_exp( -5 )
      )
    ),
    k = ( alpha + tau ) / ( 1 + exp( 5 / mu * ( Day - mu ) ) ) - tau,
    nu = theta + exp( log(epsilon - theta) - lambda * Day ),
    r = rbetapr( n() , r_mu * ( 1 + nu ) , 2 + nu )
  ) %>%
  pivot_longer(cols = c(r_mu, k, nu, r),
               names_to = "parameter") %>%
  ggplot(aes(Day, value, group = n)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    facet_wrap(~parameter, scale = "free", nrow = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# This allows for too much variability in alpha 
# but it should work given the amount of data.

# 2.1.3 Stan models ####
require(cmdstanr)
meta_c_model <- here("Decomposition", "Stan", "meta_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

meta_nc_model <- here("Decomposition", "Stan", "meta_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

meta_c_s2z_model <- here("Decomposition", "Stan", "meta_c_s2z.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

meta_nc_s2z_model <- here("Decomposition", "Stan", "meta_nc_s2z.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
meta_c_samples <- meta_c_model$sample(
          data = deco %>%
            select(Day, Ratio, Species, 
                   Experiment, Condition) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

meta_nc_samples <- meta_nc_model$sample(
          data = deco %>%
            select(Day, Ratio, Species, 
                   Experiment, Condition) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

meta_c_s2z_samples <- meta_c_s2z_model$sample(
          data = deco %>%
            select(Day, Ratio, Species, 
                   Experiment, Condition) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

meta_nc_s2z_samples <- meta_nc_s2z_model$sample(
          data = deco %>%
            select(Day, Ratio, Species, 
                   Experiment, Condition) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

# 2.1.4 Model checks ####
# Rhat
meta_c_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 93% of rhat above 1.001. rhat = 1.00 ± 0.00123.

dry_nc_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 99% of rhat above 1.001. rhat = 1.00 ± 0.0104.

# Chains
dry_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay()

dry_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Chains aren't great but somewhat better for the centred model

# Pairs
dry_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[1,1]", "log_mu[1,1]", "log_tau[1,1]"))
dry_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[1,2]", "log_mu[1,2]", "log_tau[1,2]"))

dry_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[2,1]", "log_mu[2,1]", "log_tau[2,1]"))
dry_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[2,2]", "log_mu[2,2]", "log_tau[2,2]"))

dry_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[1,1]", "log_mu[1,1]", "log_tau[1,1]"))
dry_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[1,2]", "log_mu[1,2]", "log_tau[1,2]"))

dry_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[2,1]", "log_mu[2,1]", "log_tau[2,1]"))
dry_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[2,2]", "log_mu[2,2]", "log_tau[2,2]"))
# Correlations look pretty similar for centred and non-centred models

# 2.1.5 Prior-posterior comparison ####
dry_c_prior <- prior_samples(
  model = dry_c_model,
  data = deco_dry_summary %>%
    select(Day, Proportion_mean, Proportion_sd,
           Species, Treatment) %>%
    compose_data()
  )
# Too many divergences, cannot effectively sample centred hierarchical priors

dry_nc_prior <- prior_samples(
  model = dry_nc_model,
  data = deco_dry_summary %>%
    select(Day, Proportion_mean, Proportion_sd,
           Species, Treatment) %>%
    compose_data()
  )
# Works smoothly -> use only non-centred priors, because they are the same

dry_nc_prior %>% 
  prior_posterior_draws(
    posterior_samples = dry_c_samples,
    group = deco_dry_summary %>%
      select(Species, Treatment),
    parameters = c("alpha_mu", "alpha_sigma",
                   "alpha[Species, Treatment]", 
                   "log_mu_mu", "log_mu_sigma",
                   "log_mu[Species, Treatment]", 
                   "log_tau_mu", "log_tau_sigma",
                   "log_tau[Species, Treatment]",
                   "log_epsilon_mu", "log_epsilon_sigma",
                   "log_epsilon[Species, Treatment]", 
                   "log_lambda_mu", "log_lambda_sigma",
                   "log_lambda[Species, Treatment]",
                   "log_theta_mu", "log_theta_sigma",
                   "log_theta[Species, Treatment]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Species", 
    second_group_name = "Treatment"
  )

dry_nc_prior %>% 
  prior_posterior_draws(
    posterior_samples = dry_nc_samples,
    group = deco_dry_summary %>%
      select(Species, Treatment),
    parameters = c("alpha_mu", "alpha_sigma",
                   "alpha[Species, Treatment]", 
                   "log_mu_mu", "log_mu_sigma",
                   "log_mu[Species, Treatment]", 
                   "log_tau_mu", "log_tau_sigma",
                   "log_tau[Species, Treatment]",
                   "log_epsilon_mu", "log_epsilon_sigma",
                   "log_epsilon[Species, Treatment]", 
                   "log_lambda_mu", "log_lambda_sigma",
                   "log_lambda[Species, Treatment]",
                   "log_theta_mu", "log_theta_sigma",
                   "log_theta[Species, Treatment]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Species", 
    second_group_name = "Treatment"
  )
# Posteriors are somewhat smoother for the centred model
# Since both models were suboptimal and the non-centred
# model has a better chance of improving with slower
# sampling, I will re-run the non-centred model

# 2.1.6 Rerun optimal model ####
dry_nc_samples <- dry_nc_model$sample(
          data = deco_dry_summary %>%
            select(Day, Proportion_mean, Proportion_sd,
                   Species, Treatment) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4,
          adapt_delta = 0.99, # force slower sampling to reduce divergences
          max_treedepth = 15
        ) %T>%
  print(max_rows = 200)
# Patience!

dry_nc_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 1% of rhat above 1.001. rhat = 1.00 ± 0.000159.

dry_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Chains are very good

# 2.1.7 Prediction ####
# Global parameters
dry_prior_posterior_global <- dry_nc_prior %>% 
  prior_posterior_draws(
    posterior_samples = dry_nc_samples,
    parameters = c("alpha_mu", "alpha_sigma",
                   "log_mu_mu", "log_mu_sigma",
                   "log_tau_mu", "log_tau_sigma",
                   "log_epsilon_mu", "log_epsilon_sigma",
                   "log_lambda_mu", "log_lambda_sigma",
                   "log_theta_mu", "log_theta_sigma"),
    format = "short"
  ) %>%
  mutate( # Calculate parameters for unobserved species/treatments
    alpha = rnorm( n() , alpha_mu , alpha_sigma ),
    mu = rnorm( n() , log_mu_mu , log_mu_sigma ) %>% exp(),
    tau = rnorm( n() , log_tau_mu , log_tau_sigma ) %>% exp(),
    epsilon = rnorm( n() , log_epsilon_mu , log_epsilon_sigma ) %>% exp(),
    lambda = rnorm( n() , log_lambda_mu , log_lambda_sigma ) %>% exp(),
    theta = rnorm( n() , log_theta_mu , log_theta_sigma ) %>% exp()
  ) %>%
  select(starts_with("."), distribution, 
         alpha, mu, tau, epsilon, lambda, theta) %T>%
  print()

dry_prior_posterior_global %>%
  pivot_longer(cols = -c(starts_with("."), distribution),
               names_to = "parameter") %>%
  group_by(distribution, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n())

# Species/treatment parameters
dry_prior_posterior <- dry_nc_prior %>% 
  prior_posterior_draws(
    posterior_samples = dry_nc_samples,
    group = deco_dry_summary %>%
      select(Species, Treatment),
    parameters = c("alpha[Species, Treatment]", 
                   "log_mu[Species, Treatment]", 
                   "log_tau[Species, Treatment]", 
                   "log_epsilon[Species, Treatment]", 
                   "log_lambda[Species, Treatment]", 
                   "log_theta[Species, Treatment]"),
    format = "short"
  ) %>% 
  mutate(across( # Exponentiate all logged parameters
    starts_with("log"), ~ exp(.x), .names = "{sub('^log_', '', .col)}"
  )) %>% # Remove redundant priors within species
  filter(!(Treatment == "Dead" & distribution == "prior")) %>%
  mutate( # Embed prior in treatment
    Treatment = if_else(
      distribution == "prior", "Prior", Treatment
    ) %>% fct()
  ) %>%
  select(starts_with("."), Species, Treatment, 
         alpha, mu, tau, epsilon, lambda, theta) %T>%
  print()

dry_prior_posterior %>%
  pivot_longer(cols = -c(starts_with("."), Species, Treatment),
               names_to = "parameter") %>%
  group_by(Species, Treatment, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n()) %>%
  print(n = 36)
# For dead E. radiata initial decay (0.147 d^-1) is faster than 
# final decay (0.0186 d^-1), suggesting that the logistic function
# describes an decrease in decay with time. The direction of the 
# logistic slope is not enforced in the model but in the normal 
# case where detrital photosynthesis is present, there is always 
# an increase in decay with time.

# Calculate rounded values for supplementary table
dry_parameters <- dry_prior_posterior %>%
  select(!starts_with(".")) %>%
  filter(Treatment != "Prior") %>%
  group_by(Species, Treatment) %>%
  summarise(
    across( everything(), list(mean = mean, sd = sd) ),
    n = n()
  ) %>%
  ungroup() %>%
  mutate( # Note I am converting dimensionless rates to %
    alpha = glue("{signif(alpha_mean*100, 2)} ± {signif(alpha_sd*100, 2)}"),
    mu = glue("{signif(mu_mean, 2)} ± {signif(mu_sd, 2)}"),
    tau = glue("{signif(tau_mean*100, 2)} ± {signif(tau_sd*100, 2)}"),
    epsilon = glue("{signif(epsilon_mean, 5)} ± {signif(epsilon_sd, 5)}"),
    lambda = glue("{signif(lambda_mean, 2)} ± {signif(lambda_sd, 2)}"),
    theta_mean_rounded = case_when(
      theta_mean < 100 ~ signif(theta_mean, 2), 
      theta_mean < 1e3 ~ signif(theta_mean, 3),
      theta_mean < 1e4 ~ signif(theta_mean, 4)
    ),
    theta_sd_rounded = case_when(
      theta_sd < 100 ~ signif(theta_sd, 2), 
      theta_sd < 1e3 ~ signif(theta_sd, 3),
      theta_sd < 1e4 ~ signif(theta_sd, 4)
    ),
    theta = glue("{theta_mean_rounded} ± {theta_sd_rounded}")
  ) %>%
  select(!(contains("mean") | contains("sd"))) %T>%
  print()
# Species               Treatment     n alpha       mu       tau         epsilon       lambda      theta      
# Ecklonia radiata      Live      80000 -1.7 ± 0.95 35 ± 37  3.2 ± 2.8   40732 ± 19290 0.42 ± 0.1  0.49 ± 0.29
# Ecklonia radiata      Dead      80000 -15 ± 1.9   31 ± 4.8 1.9 ± 1.1   39689 ± 18998 1.3 ± 0.79  8.5 ± 4.6  
# Amphibolis griffithii Live      80000 -1.1 ± 0.38 44 ± 56  1.2 ± 1.3   41924 ± 20599 0.18 ± 0.24 1772 ± 3279
# Amphibolis griffithii Dead      80000 -2.8 ± 1.1  20 ± 9.6 0.77 ± 0.22 42161 ± 21956 0.19 ± 0.24 1804 ± 3334

# Predict across predictor range
dry_prediction <- dry_prior_posterior %>%
  spread_continuous(data = deco_dry_summary %>%
                      # Ensure predictor range starts at 0
                      bind_rows(
                        expand_grid(
                          Day = 0, 
                          Proportion_mean = 1,
                          Species = c("Ecklonia radiata", "Amphibolis griffithii"),
                          Treatment = c("Live", "Dead")
                        ) 
                      ), 
                    # all groups have the same predictor range
                    predictor_name = "Day", length = 200) %>%
  mutate(
    p_mu = exp(
      Day * alpha - ( alpha + tau ) * mu / 5 * (
        log1p_exp( 5 / mu * ( Day - mu ) ) - log1p_exp( -5 )
      )
    ),
    k = ( alpha + tau ) / ( 1 + exp( 5 / mu * ( Day - mu ) ) ) - tau,
    nu = ( epsilon - theta ) * exp( -lambda * Day ) + theta,
    p = rbetapr( n() , p_mu * ( 1 + nu ) , 2 + nu )
  ) %>% # Summarise predictions
  group_by(Day, Species, Treatment) %>%
  median_qi(p_mu, k, nu, p, .width = c(.5, .8, .9)) %T>%
  print()

# Save progress and clean up
dry_prior_posterior %>%
  write_rds(here("Decomposition", "RDS", "dry_prior_posterior.rds"))
dry_prediction %>%
  write_rds(here("Decomposition", "RDS", "dry_prediction.rds"))

rm(list = ls(pattern = "ratio"), dry_c_model, dry_c_prior, dry_c_samples,
   dry_nc_model, dry_nc_prior, dry_nc_samples)










# 2.2 Naive model ####

