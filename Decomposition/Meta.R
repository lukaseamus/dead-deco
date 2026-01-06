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
deco %>% filter(Day != 0) %>% nrow() # 2093 observations

deco %>%
  distinct(Reference, Location, Species, Treatment, 
           Condition) %>%
  nrow() # 185 experiments

deco %>%
  count(Reference, Location, Species, Treatment, 
        Condition, Experiment) %>%
  print(n = 185)

deco %>%
  distinct(Reference, Location, Species, Treatment, 
           Condition, Experiment) %>%
  count(Experiment) %>%
  filter(n > 1)
# No duplicated experiment numbers 
           
deco %$% nlevels(Species) # 14 species
deco %>% 
  distinct(Species) %>% 
  mutate(Species = Species %>% fct_relevel(sort)) %>%
  arrange(Species)

deco %$% nlevels(Reference) # 25 references
deco %>% 
  distinct(Reference) %>%
  mutate(Reference = Reference %>% fct_relevel(sort)) %>%
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
# Quite sparse spread across time for some experiments. 
# The best strategy seems to be to use partial pooling for
# species and experiments, but within condition. This is
# a hierarchical model that does not share variance between
# conditions. If that doesn't work only partial pool across
# species.

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
# this R project and Examples.R in github.com/lukaseamus/limbodeco).
# I need to enforce an increase in exponential decay over time
# due to the known decline in physiology. Therefore I am re-
# parameterising the model in terms of delta, mu and tau.
tibble(n = 1:1e3,
       log_delta_mu = rnorm( 1e3 , log(0.05) , 0.4 ), # delta = alpha + tau
       log_mu_mu = rnorm( 1e3 , log(50) , 0.4 ),
       log_tau_mu = rnorm( 1e3 , log(0.1) , 0.4 ),
       log_delta_sigma_s = rtnorm( 1e3 , 0 , 0.3 , 0 ), # half-normal priors
       log_mu_sigma_s = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       log_tau_sigma_s = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       log_delta_sigma_e = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       log_mu_sigma_e = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       log_tau_sigma_e = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       epsilon = rgamma( 1e3 , 4e4^2 / 2e4^2 , 4e4 / 2e4^2 ),
       lambda = rexp( 1e3 , 1 ),
       theta = rgamma( 1e3 , 500^2 / 250^2 , 500 / 250^2 ),
       delta = exp(
         rnorm( 1e3 , log_delta_mu , log_delta_sigma_s ) +
           rnorm( 1e3 , 0 , log_delta_sigma_e )
       ),
       mu = exp(
         rnorm( 1e3 , log_mu_mu , log_mu_sigma_s ) +
           rnorm( 1e3 , 0 , log_mu_sigma_e )
       ),
       tau = exp(
         rnorm( 1e3 , log_tau_mu , log_tau_sigma_s ) +
           rnorm( 1e3 , 0 , log_tau_sigma_e )
       ),
       alpha = delta - tau) %>%
  expand_grid(Day = deco %$% seq(min(Day), max(Day), length.out = 100)) %>%
  mutate(
    r_mu = exp(
      Day * alpha - ( alpha + tau ) * mu / 5 * (
        log1p_exp( 5 / mu * ( Day - mu ) ) - log1p_exp( -5 )
      )
    ),
    k = ( alpha + tau ) / ( 1 + exp( 5 / mu * ( Day - mu ) ) ) - tau,
    nu = theta + (epsilon - theta) * exp( -lambda * Day ),
    r = rbetapr( n() , r_mu * ( 1 + nu ) , 2 + nu )
  ) %>%
  pivot_longer(cols = c(r_mu, k, nu, r),
               names_to = "parameter") %>%
  ggplot(aes(Day, value, group = n)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    facet_wrap(~parameter, scale = "free", nrow = 1) +
    facetted_pos_scales(
      y = list(
        parameter == "r" ~ scale_y_continuous(limits = c(0, 2)),
        parameter == "r_mu" ~ scale_y_continuous(limits = c(0, 2))
      )
    ) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# With the exception of a bit of runaway exponential growth, this
# looks like reasonable prior uncertainty.

# 2.1.3 Stan models ####
# I previously tried to code a hierarchical model design for a
# simple intercept model and it struggled unless I enforced
# a sum-to-zero constraint on the partially pooled variables. 
# Here, I have opted to run separate models per condition with
# the same priors. This is mathematically the same because 
# information about interspecific and inter-experiment variability 
# is not shared between conditions in either. The benefit of 
# separate models is that for such a complex case as this, 
# I can partition computational load, make the model more readable 
# and avoid using sum-to-zero constraints.

require(cmdstanr)
meta_c_model <- here("Decomposition", "Stan", "meta_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

meta_nc_model <- here("Decomposition", "Stan", "meta_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
meta_live_c_samples <- meta_c_model$sample(
          data = deco %>%
            filter(Condition == "Live" & # Note that I filter for live
                     Day != 0) %>% # t0 = 1 is pre-determined
            droplevels() %>% # This makes sure no unused levels are floating about
            select(Day, Ratio, Species, Experiment) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)
# All chains fail because of large prior uncertainty in epsilon and theta
# paired with the log( epsilon - theta ) structure lead to the log probability
# evaluating to -Inf. However, the alternative formulation leads to overflow in 
# this case, which in turn causes Inf Rhat and NA effective sample size.
# Let's run the non-centred version.

meta_live_nc_samples <- meta_nc_model$sample(
          data = deco %>%
            filter(Condition == "Live" & Day != 0) %>%
            droplevels() %>%
            select(Day, Ratio, Species, Experiment) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

meta_dead_c_samples <- meta_c_model$sample(
          data = deco %>%
            filter(Condition == "Dead" & Day != 0) %>% # Note that I filter for dead
            droplevels() %>%
            select(Day, Ratio, Species, Experiment) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

meta_dead_nc_samples <- meta_nc_model$sample(
          data = deco %>%
            filter(Condition == "Dead" & Day != 0) %>%
            droplevels() %>%
            select(Day, Ratio, Species, Experiment) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

# 2.1.4 Model checks ####
# Rhat
meta_live_c_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 65% of rhat above 1.001. rhat = 1.02 ± 0.0236.

meta_live_nc_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 15% of rhat above 1.001. rhat = 1.00 ± 0.00107.

meta_dead_c_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 100% of rhat above 1.001. rhat = 1.01 ± 0.00740.

meta_dead_nc_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 7% of rhat above 1.001. rhat = 1.00 ± 0.000318.
# The non-centred models seem better.

# Chains
require(bayesplot)
meta_live_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Two chains were lost

meta_live_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# No chains lost, chains look ok

meta_dead_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Bad chains

meta_dead_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Good look ok

# Pairs
meta_live_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta_mu", "log_mu_mu", "log_tau_mu"))
meta_live_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta_s[1]", "log_delta_e[1]"))
meta_live_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_mu_s[1]", "log_mu_e[1]"))
meta_live_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_tau_s[1]", "log_tau_e[1]"))

meta_live_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta_mu", "log_mu_mu", "log_tau_mu"))
meta_live_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta_s[1]", "log_delta_e[1]"))
meta_live_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_mu_s[1]", "log_mu_e[1]"))
meta_live_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_tau_s[1]", "log_tau_e[1]"))

meta_dead_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta_mu", "log_mu_mu", "log_tau_mu"))
meta_dead_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta_s[1]", "log_delta_e[1]"))
meta_dead_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_mu_s[1]", "log_mu_e[1]"))
meta_dead_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_tau_s[1]", "log_tau_e[1]"))

meta_dead_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta_mu", "log_mu_mu", "log_tau_mu"))
meta_dead_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta_s[1]", "log_delta_e[1]"))
meta_dead_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_mu_s[1]", "log_mu_e[1]"))
meta_dead_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_tau_s[1]", "log_tau_e[1]"))
# No or weak correlations in all cases

# 2.1.5 Prior-posterior comparison ####
source("functions.R")
meta_live_prior <- prior_samples(
  model = meta_nc_model,
  data = deco %>%
    filter(Condition == "Live" & Day != 0) %>%
    droplevels() %>%
    select(Day, Ratio, Species, Experiment) %>%
    compose_data()
  )

meta_dead_prior <- prior_samples(
  model = meta_nc_model,
  data = deco %>%
    filter(Condition == "Dead" & Day != 0) %>%
    droplevels() %>%
    select(Day, Ratio, Species, Experiment) %>%
    compose_data()
)

# Live centred species
meta_live_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_live_c_samples,
    group = deco %>%
            filter(Condition == "Live" & Day != 0) %>%
            droplevels() %>%
            select(Species),
    parameters = c("log_delta_mu", "log_delta_sigma_s", "log_delta_s[Species]",
                   "log_mu_mu", "log_mu_sigma_s", "log_mu_s[Species]",
                   "log_tau_mu", "log_tau_sigma_s", "log_tau_s[Species]",
                   "epsilon", "lambda", "theta"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Species")

# Live centred experiments
meta_live_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_live_c_samples,
    group = deco %>%
            filter(Condition == "Live" & Day != 0) %>%
            droplevels() %>%
            select(Experiment),
    parameters = c("log_delta_mu", "log_delta_sigma_e", "log_delta_e[Experiment]",
                   "log_mu_mu", "log_mu_sigma_e", "log_mu_e[Experiment]",
                   "log_tau_mu", "log_tau_sigma_e", "log_tau_e[Experiment]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Experiment")

# Live non-centred species
meta_live_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_live_nc_samples,
    group = deco %>%
            filter(Condition == "Live" & Day != 0) %>%
            droplevels() %>%
            select(Species),
    parameters = c("log_delta_mu", "log_delta_sigma_s", "log_delta_s[Species]",
                   "log_mu_mu", "log_mu_sigma_s", "log_mu_s[Species]",
                   "log_tau_mu", "log_tau_sigma_s", "log_tau_s[Species]",
                   "epsilon", "lambda", "theta"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Species")

# Live non-centred experiments
meta_live_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_live_nc_samples,
    group = deco %>%
      filter(Condition == "Live" & Day != 0) %>%
      droplevels() %>%
      select(Experiment),
    parameters = c("log_delta_mu", "log_delta_sigma_e", "log_delta_e[Experiment]",
                   "log_mu_mu", "log_mu_sigma_e", "log_mu_e[Experiment]",
                   "log_tau_mu", "log_tau_sigma_e", "log_tau_e[Experiment]"),
    format = "long"
  ) %>%
  prior_posterior_plot(group_name = "Experiment")
# Non-centred posteriors look better

# Dead centred species and experiments
meta_dead_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_dead_c_samples,
    group = deco %>%
            filter(Condition == "Dead" & Day != 0) %>%
            droplevels() %>%
            select(Species, Experiment),
    parameters = c("log_delta_mu", 
                   "log_delta_sigma_s", "log_delta_sigma_e",
                   "log_delta_s[Species]", "log_delta_e[Experiment]",
                   "log_mu_mu", 
                   "log_mu_sigma_s", "log_mu_sigma_e",
                   "log_mu_s[Species]", "log_mu_e[Experiment]",
                   "log_tau_mu", 
                   "log_tau_sigma_s", "log_tau_sigma_e",
                   "log_tau_s[Species]", "log_tau_e[Experiment]",
                   "epsilon", "lambda", "theta"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Species",
    second_group_name = "Experiment"
  )

# Dead non-centred species and experiments
meta_dead_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_dead_nc_samples,
    group = deco %>%
            filter(Condition == "Dead" & Day != 0) %>%
            droplevels() %>%
            select(Species, Experiment),
    parameters = c("log_delta_mu", 
                   "log_delta_sigma_s", "log_delta_sigma_e",
                   "log_delta_s[Species]", "log_delta_e[Experiment]",
                   "log_mu_mu", 
                   "log_mu_sigma_s", "log_mu_sigma_e",
                   "log_mu_s[Species]", "log_mu_e[Experiment]",
                   "log_tau_mu", 
                   "log_tau_sigma_s", "log_tau_sigma_e",
                   "log_tau_s[Species]", "log_tau_e[Experiment]",
                   "epsilon", "lambda", "theta"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Species",
    second_group_name = "Experiment"
  )
# Non-centred posteriors look better

# Choose non-cetred model as optimal in both cases

# 2.1.6 Prediction ####
# Global parameters
meta_live_prior_posterior_global <- meta_live_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_live_nc_samples,
    parameters = c("log_delta_mu", "log_delta_sigma_s", "log_delta_sigma_e",
                   "log_mu_mu", "log_mu_sigma_s", "log_mu_sigma_e",
                   "log_tau_mu", "log_tau_sigma_s", "log_tau_sigma_e",
                   "epsilon", "lambda", "theta"),
    format = "short"
  ) %>%
  mutate( # Calculate parameters for new species and experiments
    delta = exp(
      rnorm( n() , log_delta_mu , log_delta_sigma_s ) +
        rnorm( n() , 0 , log_delta_sigma_e )
    ),
    mu = exp(
      rnorm( n() , log_mu_mu , log_mu_sigma_s ) +
        rnorm( n() , 0 , log_mu_sigma_e )
    ),
    tau = exp(
      rnorm( n() , log_tau_mu , log_tau_sigma_s ) +
        rnorm( n() , 0 , log_tau_sigma_e )
    ),
    alpha = delta - tau
  ) %>%
  rowwise() %>% # Simulate a distribution of alpha per row and 
  mutate( # calculate the species and experiment standard deviations.
    alpha_sigma_s = sd(
      exp( rnorm( 5e3 , log_delta_mu , log_delta_sigma_s ) ) - 
        exp( rnorm( 5e3 , log_tau_mu , log_tau_sigma_s ) )
    ),
    alpha_sigma_e = sd(
      exp( rnorm( 5e3 , log_delta_mu , log_delta_sigma_e ) ) -
        exp( rnorm( 5e3 , log_tau_mu , log_tau_sigma_e ) )
    )
  ) %>%
  ungroup() %>%
  select(starts_with("."), distribution, 
         alpha, mu, tau, epsilon, lambda, theta,
         alpha_sigma_s, alpha_sigma_e,
         log_mu_sigma_s, log_mu_sigma_e,
         log_tau_sigma_s, log_tau_sigma_e) %T>%
  print()

meta_dead_prior_posterior_global <- meta_dead_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_dead_nc_samples,
    parameters = c("log_delta_mu", "log_delta_sigma_s", "log_delta_sigma_e",
                   "log_mu_mu", "log_mu_sigma_s", "log_mu_sigma_e",
                   "log_tau_mu", "log_tau_sigma_s", "log_tau_sigma_e",
                   "epsilon", "lambda", "theta"),
    format = "short"
  ) %>%
  mutate(
    delta = exp(
      rnorm( n() , log_delta_mu , log_delta_sigma_s ) +
        rnorm( n() , 0 , log_delta_sigma_e )
    ),
    mu = exp(
      rnorm( n() , log_mu_mu , log_mu_sigma_s ) +
        rnorm( n() , 0 , log_mu_sigma_e )
    ),
    tau = exp(
      rnorm( n() , log_tau_mu , log_tau_sigma_s ) +
        rnorm( n() , 0 , log_tau_sigma_e )
    ),
    alpha = delta - tau
  ) %>%
  rowwise() %>% 
  mutate(
    alpha_sigma_s = sd(
      exp( rnorm( 5e3 , log_delta_mu , log_delta_sigma_s ) ) - 
        exp( rnorm( 5e3 , log_tau_mu , log_tau_sigma_s ) )
    ),
    alpha_sigma_e = sd(
      exp( rnorm( 5e3 , log_delta_mu , log_delta_sigma_e ) ) -
        exp( rnorm( 5e3 , log_tau_mu , log_tau_sigma_e ) )
    )
  ) %>%
  ungroup() %>%
  select(starts_with("."), distribution, 
         alpha, mu, tau, epsilon, lambda, theta,
         alpha_sigma_s, alpha_sigma_e,
         log_mu_sigma_s, log_mu_sigma_e,
         log_tau_sigma_s, log_tau_sigma_e) %T>%
  print()

meta_prior_posterior_global <- meta_live_prior_posterior_global %>%
  mutate(Condition = "Live" %>% fct()) %>%
  bind_rows(
    meta_dead_prior_posterior_global %>%
      mutate(Condition = "Dead" %>% fct())
  ) %>% # Priors for both conditions are identical
  filter(!(Condition == "Dead" & distribution == "prior")) %>%
  mutate(
    # Embed prior in condition
    Condition = if_else(
      distribution == "prior", "Prior", Condition
    ) %>% fct()
  ) %>%
  select(-distribution) %T>%
  print()

meta_prior_posterior_global %>%
  pivot_longer(cols = -c(starts_with("."), Condition),
               names_to = "parameter") %>%
  group_by(Condition, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n()) %>%
  print(n = 36)

# Species parameters
meta_live_prior_posterior_species <- meta_live_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_live_nc_samples,
    group = deco %>%
            filter(Condition == "Live" & Day != 0) %>%
            droplevels() %>%
            select(Species),
    parameters = c("log_delta_s[Species]", "log_delta_sigma_e",
                   "log_mu_s[Species]", "log_mu_sigma_e",
                   "log_tau_s[Species]", "log_tau_sigma_e"),
    format = "short"
  ) %>%
  # mutate(across( # Exponentiate all log parameters
  #   starts_with("log"), ~ exp(.x), .names = "{sub('^log_', '', .col)}"
  # )) %>%
  mutate(
    delta = exp(
      rnorm( n() , log_delta_s , log_delta_sigma_e )
    ),
    mu = exp(
      rnorm( n() , log_mu_s , log_mu_sigma_e )
    ),
    tau = exp(
      rnorm( n() , log_tau_s , log_tau_sigma_e )
    ),
    alpha = delta - tau
  ) %>%
  filter(Species == "Macrocystis pyrifera" & distribution == "prior" |
           distribution == "posterior") %>% # Remove redundant priors
  mutate( # Embed prior in species
    Species = if_else(
      distribution == "prior", "Prior", Species
    ) %>% fct()
  ) %>%
  select(starts_with("."), Species, 
         alpha, mu, tau) %T>%
  print()

meta_dead_prior_posterior_species <- meta_dead_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_dead_nc_samples,
    group = deco %>%
            filter(Condition == "Dead" & Day != 0) %>%
            droplevels() %>%
            select(Species),
    parameters = c("log_delta_s[Species]", "log_delta_sigma_e",
                   "log_mu_s[Species]", "log_mu_sigma_e",
                   "log_tau_s[Species]", "log_tau_sigma_e"),
    format = "short"
  ) %>%
  mutate(
    delta = exp(
      rnorm( n() , log_delta_s , log_delta_sigma_e )
    ),
    mu = exp(
      rnorm( n() , log_mu_s , log_mu_sigma_e )
    ),
    tau = exp(
      rnorm( n() , log_tau_s , log_tau_sigma_e )
    ),
    alpha = delta - tau
  ) %>%
  filter(Species == "Nereocystis luetkeana" & distribution == "prior" |
           distribution == "posterior") %>%
  mutate(
    Species = if_else(
      distribution == "prior", "Prior", Species
    ) %>% fct()
  ) %>%
  select(starts_with("."), Species, 
         alpha, mu, tau) %T>%
  print()

meta_prior_posterior_species <- meta_live_prior_posterior_species %>%
  mutate(Condition = "Live" %>% fct()) %>%
  bind_rows(
    meta_dead_prior_posterior_species %>%
      mutate(Condition = "Dead" %>% fct())
  ) %>%
  filter(!(Condition == "Dead" & Species == "Prior")) %>%
  mutate(
    Condition = if_else(
      Species == "Prior", "Prior", Condition
    ) %>% fct()
  ) %T>%
  print()

meta_prior_posterior_species %>%
  pivot_longer(cols = -c(starts_with("."), Species, Condition),
               names_to = "parameter") %>%
  group_by(Species, Condition, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n()) %>%
  print(n = 57)

# Experiment parameters
meta_live_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_live_nc_samples,
    group = deco %>%
            filter(Condition == "Live" & Day != 0) %>%
            droplevels() %>%
            select(Species, Experiment),
    parameters = c("log_delta_s[Species]", "log_delta_e[Experiment]",
                   "log_mu_s[Species]", "log_mu_e[Experiment]",
                   "log_tau_s[Species]", "log_tau_e[Experiment]"),
    format = "long"
  ) %>%
  distinct(Species, Experiment) %>%
  print(n = 200)
# This doesn't match experiments to their species,
# so a workaround is required.

meta_live_prior_posterior_experiment <- deco %>%
  filter(Condition == "Live" & Day != 0) %>%
  droplevels() %>%
  distinct(Species, Experiment) %>% # Get species-experiment pairs from data
  left_join( # Join the species samples by species
    meta_live_prior %>% 
      prior_posterior_draws(
        posterior_samples = meta_live_nc_samples,
        group = deco %>%
          filter(Condition == "Live" & Day != 0) %>%
          droplevels() %>%
          select(Species),
        parameters = c("log_delta_s[Species]",
                       "log_mu_s[Species]",
                       "log_tau_s[Species]"),
        format = "short"
      ),
    by = "Species",
    relationship = "many-to-many"
  ) %>%
  left_join( # Join the experiment samples by experiment and sample ID
    meta_live_prior %>% 
      prior_posterior_draws(
        posterior_samples = meta_live_nc_samples,
        group = deco %>%
          filter(Condition == "Live" & Day != 0) %>%
          droplevels() %>%
          select(Experiment),
        parameters = c("log_delta_e[Experiment]",
                       "log_mu_e[Experiment]",
                       "log_tau_e[Experiment]"),
        format = "short"
      ),
    by = c("Experiment", ".chain", ".iteration", ".draw", "distribution"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    delta = exp( log_delta_s + log_delta_e ),
    mu = exp( log_mu_s + log_mu_e ),
    tau = exp( log_tau_s + log_tau_e ),
    alpha = delta - tau
  ) %>%
  filter(Species == "Macrocystis pyrifera" & Experiment == 1 & # Remove redundant priors
           distribution == "prior" | distribution == "posterior") %>% 
  mutate( # Embed prior in species and experiment
    Species = if_else(
      distribution == "prior", "Prior", Species
    ) %>% fct(),
    Experiment = if_else(
      distribution == "prior", "Prior", Experiment
    ) %>% fct()
  ) %>%
  select(starts_with("."), Species, Experiment,
         alpha, mu, tau) %T>%
  print()
  
meta_dead_prior_posterior_experiment <- deco %>%
  filter(Condition == "Dead" & Day != 0) %>%
  droplevels() %>%
  distinct(Species, Experiment) %>%
  left_join(
    meta_dead_prior %>% 
      prior_posterior_draws(
        posterior_samples = meta_dead_nc_samples,
        group = deco %>%
          filter(Condition == "Dead" & Day != 0) %>%
          droplevels() %>%
          select(Species),
        parameters = c("log_delta_s[Species]",
                       "log_mu_s[Species]",
                       "log_tau_s[Species]"),
        format = "short"
      ),
    by = "Species",
    relationship = "many-to-many"
  ) %>%
  left_join(
    meta_dead_prior %>% 
      prior_posterior_draws(
        posterior_samples = meta_dead_nc_samples,
        group = deco %>%
          filter(Condition == "Dead" & Day != 0) %>%
          droplevels() %>%
          select(Experiment),
        parameters = c("log_delta_e[Experiment]",
                       "log_mu_e[Experiment]",
                       "log_tau_e[Experiment]"),
        format = "short"
      ),
    by = c("Experiment", ".chain", ".iteration", ".draw", "distribution"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    delta = exp( log_delta_s + log_delta_e ),
    mu = exp( log_mu_s + log_mu_e ),
    tau = exp( log_tau_s + log_tau_e ),
    alpha = delta - tau
  ) %>%
  filter(Species == "Nereocystis luetkeana" & Experiment == 14 &
           distribution == "prior" | distribution == "posterior") %>% 
  mutate(
    Species = if_else(
      distribution == "prior", "Prior", Species
    ) %>% fct(),
    Experiment = if_else(
      distribution == "prior", "Prior", Experiment
    ) %>% fct()
  ) %>%
  select(starts_with("."), Species, Experiment,
         alpha, mu, tau) %T>%
  print()

meta_prior_posterior_experiment <- meta_live_prior_posterior_experiment %>%
  mutate(Condition = "Live" %>% fct()) %>%
  bind_rows(
    meta_dead_prior_posterior_experiment %>%
      mutate(Condition = "Dead" %>% fct())
  ) %>%
  filter(!(Condition == "Dead" & Species == "Prior")) %>%
  mutate(
    Condition = if_else(
      Species == "Prior", "Prior", Condition
    ) %>% fct()
  ) %T>%
  print()

meta_prior_posterior_experiment %>%
  pivot_longer(cols = -c(starts_with("."), Species, Experiment, Condition),
               names_to = "parameter") %>%
  group_by(Species, Experiment, Condition, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n()) %>%
  print()

# Save progress and clean up
meta_prior_posterior_global %>%
  write_rds(here("Decomposition", "RDS", "meta_prior_posterior_global.rds"))
meta_prior_posterior_species %>%
  write_rds(here("Decomposition", "RDS", "meta_prior_posterior_species.rds"))
meta_prior_posterior_experiment %>%
  write_rds(here("Decomposition", "RDS", "meta_prior_posterior_experiment.rds"))

rm(meta_c_model, meta_nc_model, meta_live_prior, meta_dead_prior,
   meta_live_c_samples, meta_live_nc_samples, meta_dead_c_samples, meta_dead_nc_samples,
   meta_live_prior_posterior_global, meta_dead_prior_posterior_global,
   meta_live_prior_posterior_species, meta_dead_prior_posterior_species,
   meta_live_prior_posterior_experiment, meta_dead_prior_posterior_experiment)

# Calculate global rounded values for text
require(glue)
meta_parameters_global <- meta_prior_posterior_global %>%
  filter(Condition != "Prior") %>%
  select(!starts_with(".")) %>%
  group_by(Condition) %>%
  summarise(
    across( everything(), list(mean = mean, sd = sd, median = median) ),
    n = n()
  ) %>%
  ungroup() %>%
  mutate( # Note I am converting dimensionless rates to %
    alpha = glue("{signif(alpha_mean*100, 2)} ± {signif(alpha_sd*100, 2)}"),
    mu_mean_rounded = case_when(
      mu_mean < 100 ~ signif(mu_mean, 2), 
      mu_mean < 1e3 ~ signif(mu_mean, 3),
      TRUE ~ signif(mu_mean, 4)
    ),
    mu_sd_rounded = case_when(
      mu_sd < 100 ~ signif(mu_sd, 2), 
      mu_sd < 1e3 ~ signif(mu_sd, 3),
      TRUE ~ signif(mu_sd, 4)
    ),
    mu_median_rounded = case_when(
      mu_median < 100 ~ signif(mu_median, 2), 
      mu_median < 1e3 ~ signif(mu_median, 3),
      TRUE ~ signif(mu_median, 4)
    ),
    mu = glue("{mu_mean_rounded} ± {mu_sd_rounded} ({mu_median_rounded})"),
    tau = glue("{signif(tau_mean*100, 2)} ± {signif(tau_sd*100, 2)}"),
    # Note standard deviations are not converted
    alpha_sigma_s = glue("{signif(alpha_sigma_s_mean, 2)} ± {signif(alpha_sigma_s_sd, 2)}"),
    alpha_sigma_e = glue("{signif(alpha_sigma_e_mean, 2)} ± {signif(alpha_sigma_e_sd, 2)}"),
    log_mu_sigma_s = glue("{signif(log_mu_sigma_s_mean, 2)} ± {signif(log_mu_sigma_s_sd, 2)}"),
    log_mu_sigma_e = glue("{signif(log_mu_sigma_e_mean, 2)} ± {signif(log_mu_sigma_e_sd, 2)}"),
    log_tau_sigma_s = glue("{signif(log_tau_sigma_s_mean, 2)} ± {signif(log_tau_sigma_s_sd, 2)}"),
    log_tau_sigma_e = glue("{signif(log_tau_sigma_e_mean, 2)} ± {signif(log_tau_sigma_e_sd, 2)}"),
    epsilon_mean_rounded = case_when(
      epsilon_mean < 100 ~ signif(epsilon_mean, 2), 
      epsilon_mean < 1e3 ~ signif(epsilon_mean, 3),
      epsilon_mean < 1e4 ~ signif(epsilon_mean, 4),
      epsilon_mean < 1e5 ~ signif(epsilon_mean, 5),
      TRUE ~ signif(epsilon_mean, 6)
    ),
    epsilon_sd_rounded = case_when(
      epsilon_sd < 100 ~ signif(epsilon_sd, 2), 
      epsilon_sd < 1e3 ~ signif(epsilon_sd, 3),
      epsilon_sd < 1e4 ~ signif(epsilon_sd, 4),
      epsilon_sd < 1e5 ~ signif(epsilon_sd, 5),
      TRUE ~ signif(epsilon_sd, 6)
    ),
    epsilon = glue("{epsilon_mean_rounded} ± {epsilon_sd_rounded}"),
    lambda = glue("{signif(lambda_mean, 2)} ± {signif(lambda_sd, 2)}"),
    theta_mean_rounded = case_when(
      theta_mean < 100 ~ signif(theta_mean, 2), 
      theta_mean < 1e3 ~ signif(theta_mean, 3),
      theta_mean < 1e4 ~ signif(theta_mean, 4),
      theta_mean < 1e5 ~ signif(theta_mean, 5),
      TRUE ~ signif(theta_mean, 6)
    ),
    theta_sd_rounded = case_when(
      theta_sd < 100 ~ signif(theta_sd, 2), 
      theta_sd < 1e3 ~ signif(theta_sd, 3),
      theta_sd < 1e4 ~ signif(theta_sd, 4),
      theta_sd < 1e5 ~ signif(theta_sd, 5),
      TRUE ~ signif(theta_sd, 6)
    ),
    theta = glue("{theta_mean_rounded} ± {theta_sd_rounded}")
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# Calculate species rounded values for supplementary table
meta_parameters_species <- meta_prior_posterior_species %>%
  select(!starts_with(".")) %>%
  filter(Species != "Prior") %>%
  group_by(Species, Condition) %>%
  summarise(
    across( everything(), list(mean = mean, sd = sd, median = median) ),
    n = n()
  ) %>%
  ungroup() %>%
  mutate( # Note I am converting dimensionless rates to %
    alpha = glue("{signif(alpha_mean*100, 2)} ± {signif(alpha_sd*100, 2)}"),
    mu_mean_rounded = case_when(
      mu_mean < 100 ~ signif(mu_mean, 2), 
      mu_mean < 1e3 ~ signif(mu_mean, 3),
      TRUE ~ signif(mu_mean, 4)
    ),
    mu_sd_rounded = case_when(
      mu_sd < 100 ~ signif(mu_sd, 2), 
      mu_sd < 1e3 ~ signif(mu_sd, 3),
      TRUE ~ signif(mu_sd, 4)
    ),
    mu_median_rounded = case_when(
      mu_median < 100 ~ signif(mu_median, 2), 
      mu_median < 1e3 ~ signif(mu_median, 3),
      TRUE ~ signif(mu_median, 4)
    ),
    mu = glue("{mu_mean_rounded} ± {mu_sd_rounded} ({mu_median_rounded})"),
    tau = glue("{signif(tau_mean*100, 2)} ± {signif(tau_sd*100, 2)}")
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# Calculate experiment rounded values for text
meta_parameters_experiment <- meta_prior_posterior_experiment %>%
  select(!starts_with(".")) %>%
  filter(Species != "Prior") %>%
  group_by(Species, Experiment, Condition) %>%
  summarise(
    across( everything(), list(mean = mean, sd = sd, median = median) ),
    n = n()
  ) %>%
  ungroup() %>%
  mutate( # Note I am converting dimensionless rates to %
    alpha = glue("{signif(alpha_mean*100, 2)} ± {signif(alpha_sd*100, 2)}"),
    mu_mean_rounded = case_when(
      mu_mean < 100 ~ signif(mu_mean, 2), 
      mu_mean < 1e3 ~ signif(mu_mean, 3),
      TRUE ~ signif(mu_mean, 4)
    ),
    mu_sd_rounded = case_when(
      mu_sd < 100 ~ signif(mu_sd, 2), 
      mu_sd < 1e3 ~ signif(mu_sd, 3),
      TRUE ~ signif(mu_sd, 4)
    ),
    mu_median_rounded = case_when(
      mu_median < 100 ~ signif(mu_median, 2), 
      mu_median < 1e3 ~ signif(mu_median, 3),
      TRUE ~ signif(mu_median, 4)
    ),
    mu = glue("{mu_mean_rounded} ± {mu_sd_rounded} ({mu_median_rounded})"),
    tau = glue("{signif(tau_mean*100, 2)} ± {signif(tau_sd*100, 2)}")
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# Calculate sigma contrasts
meta_sigma <- meta_prior_posterior_global %>%
  filter(Condition != "Prior") %>%
  select(Condition, starts_with("."), contains("sigma")) %>%
  pivot_longer(cols = contains("sigma"),
               names_to = "Parameter") %>%
  pivot_wider(names_from = Condition, values_from = value) %>%
  mutate(difference = Live - Dead,
         ratio = Live / Dead,
         log_ratio = log10(ratio)) %T>%
  print()

# Summarise contrasts for main table
meta_sigma_summary <- meta_sigma %>%
  select(!starts_with(".")) %>%
  group_by(Parameter) %>%
  summarise(
    across( everything(), list(mean = mean, sd = sd) ),
    n = n(),
    P = mean( difference > 0 ) %>% signif(2)
  ) %>%
  ungroup() %>%
  mutate(
    Live = glue("{signif(Live_mean, 2)} ± {signif(Live_sd, 2)}"),
    Dead = glue("{signif(Dead_mean, 2)} ± {signif(Dead_sd, 2)}"),
    difference = glue("{signif(difference_mean, 2)} ± {signif(difference_sd, 2)}"),
    log_ratio = glue("{signif(log_ratio_mean, 2)} ± {signif(log_ratio_sd, 2)}")
  ) %>%
  select(!(contains("mean") | contains("sd"))) %T>%
  print()

# Predict across predictor range
# Predict for new experiments on new species
meta_prediction_global <- meta_prior_posterior_global %>%
  spread_continuous(data = deco, 
                    predictor_name = "Day", 
                    group_name = "Condition") %>% 
  mutate(
    r_mu = exp(
      Day * alpha - ( alpha + tau ) * mu / 5 * (
        log1p_exp( 5 / mu * ( Day - mu ) ) - log1p_exp( -5 )
      )
    ),
    k = ( alpha + tau ) / ( 1 + exp( 5 / mu * ( Day - mu ) ) ) - tau,
    nu = ( epsilon - theta ) * exp( -lambda * Day ) + theta,
    r = rbetapr( n() , r_mu * ( 1 + nu ) , 2 + nu )
  ) %>%
  drop_na() %T>%
  print()
# There was some underflow to r_mu = 0 because of tau > 4, 
# which results in r = NA. Dropping these NAs before summarising
# avoids the summary becoming NA due to a few extreme samples.
# Warning still shows but is safe to ignore since NAs were filtered out.

meta_prediction_global %>%
  group_by(Day, Condition) %>%
  count() %>%
  filter(n < 8e4) %>%
  print(n = 30)
# Only 23 samples were removed. All groups have n >= 79999.
# Only Prior is affected.

# Summarise predictions
meta_prediction_global_summary <- meta_prediction_global %>%
  group_by(Day, Condition) %>%
  median_qi(r_mu, k, nu, r, .width = c(.5, .8, .9)) %T>%
  print()

rm(meta_prediction_global)

# Predict for new experiments on existing species
meta_prediction_species <- meta_prior_posterior_species %>%
  spread_continuous(data = deco, 
                    predictor_name = "Day",
                    group_name = "Species",
                    second_group_name = "Condition") %>%
  mutate(
    r_mu = exp(
      Day * alpha - ( alpha + tau ) * mu / 5 * (
        log1p_exp( 5 / mu * ( Day - mu ) ) - log1p_exp( -5 )
      )
    ),
    k = ( alpha + tau ) / ( 1 + exp( 5 / mu * ( Day - mu ) ) ) - tau
  ) %T>%
  print()

# Summarise predictions

# meta_prediction_species_summary <- meta_prediction_species %>%
#   group_by(Day, Species, Condition) %>%
#   median_qi(r_mu, k, .width = c(.5, .8, .9)) %T>%
#   print()

# median_qi as above crashes R. Use custom function instead.

meta_prediction_species_summary <- meta_prediction_species %>%
  group_by(Day, Species, Condition) %>%
  summarise(
    across(
      c(r_mu, k),
      list(median = median, 
           lower_0.5 = ~ qi(.x, .width = .5)[1],
           upper_0.5 = ~ qi(.x, .width = .5)[2],
           lower_0.8 = ~ qi(.x, .width = .8)[1],
           upper_0.8 = ~ qi(.x, .width = .8)[2],
           lower_0.9 = ~ qi(.x, .width = .9)[1],
           upper_0.9 = ~ qi(.x, .width = .9)[2]),
      .names = "{.col}.{.fn}"
    )
  ) %>%
  ungroup() %>%
  rename(r_mu = r_mu.median, k = k.median) %>%
  pivot_longer(cols = contains("lower") | contains("upper")) %>%
  separate(col = name, into = c("name", ".width"), sep = "_(?=[^_]*$)") %>%
  mutate(.width = .width %>% as.numeric()) %>%
  pivot_wider(names_from = name, values_from = value) %T>%
  print()

rm(meta_prediction_species)

# Predict for existing experiments
meta_prediction_experiment <- meta_prior_posterior_experiment %>% 
  group_by(Experiment, Condition, .chain) %>%
  slice_sample(n = 1500) %>% # Reduce the number of samples per group and chain
  ungroup() %>%
  spread_continuous(data = deco, 
                    predictor_name = "Day",
                    group_name = "Experiment",
                    second_group_name = "Condition") %>%
  mutate(
    r_mu = exp(
      Day * alpha - ( alpha + tau ) * mu / 5 * (
        log1p_exp( 5 / mu * ( Day - mu ) ) - log1p_exp( -5 )
      )
    ),
    k = ( alpha + tau ) / ( 1 + exp( 5 / mu * ( Day - mu ) ) ) - tau
  ) %T>%
  print()

# Summarise predictions
meta_prediction_experiment_summary <- meta_prediction_experiment %>%
  group_by(Day, Species, Experiment, Condition) %>%
  summarise(
    across(
      c(r_mu, k),
      list(median = median, 
           lower_0.5 = ~ qi(.x, .width = .5)[1],
           upper_0.5 = ~ qi(.x, .width = .5)[2],
           lower_0.8 = ~ qi(.x, .width = .8)[1],
           upper_0.8 = ~ qi(.x, .width = .8)[2],
           lower_0.9 = ~ qi(.x, .width = .9)[1],
           upper_0.9 = ~ qi(.x, .width = .9)[2]),
      .names = "{.col}.{.fn}"
    )
  ) %>%
  ungroup() %>%
  rename(r_mu = r_mu.median, k = k.median) %>%
  pivot_longer(cols = contains("lower") | contains("upper")) %>%
  separate(col = name, into = c("name", ".width"), sep = "_(?=[^_]*$)") %>%
  mutate(.width = .width %>% as.numeric()) %>%
  pivot_wider(names_from = name, values_from = value) %T>%
  print()

rm(meta_prediction_experiment)

# Save progress
meta_prediction_global_summary %>%
  write_rds(here("Decomposition", "RDS", "meta_prediction_global.rds"))
meta_prediction_species_summary %>%
  write_rds(here("Decomposition", "RDS", "meta_prediction_species.rds"))
meta_prediction_experiment_summary %>%
  write_rds(here("Decomposition", "RDS", "meta_prediction_experiment.rds"))

# 2.2 Naive model ####
# 2.2.1 Prior simulation ####
# The classic exponential decay model for proportions is e^-k*t.
# I am again taking 0.1 d^-1 as my prior.
tibble(n = 1:1e3, # Much higher prior uncertainty because there's only one parameter
       log_k_mu = rnorm( 1e3 , log(0.1) , 0.6 ), 
       log_k_sigma_s = rtnorm( 1e3 , 0 , 0.6 , 0 ), # half-normal prior
       log_k_sigma_e = rtnorm( 1e3 , 0 , 0.6 , 0 ),
       k = exp(
         rnorm( 1e3 , log_k_mu , log_k_sigma_s ) +
           rnorm( 1e3 , 0 , log_k_sigma_e )
       ),
       sigma = rexp( 1e3 , 1 )) %>%
  expand_grid(Day = deco %$% seq(min(Day), max(Day), length.out = 100)) %>%
  mutate(
    r_mu = exp( -k * Day ),
    r = rnorm( n() , r_mu , sigma )
  ) %>%
  pivot_longer(cols = c(r_mu, r),
               names_to = "parameter") %>%
  ggplot(aes(Day, value, group = n)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    facet_wrap(~parameter, scale = "free", nrow = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 2.2.2 Stan models ####
meta_k_c_model <- here("Decomposition", "Stan", "meta_k_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

meta_k_nc_model <- here("Decomposition", "Stan", "meta_k_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

meta_live_k_c_samples <- meta_k_c_model$sample(
          data = deco %>%
            filter(Condition == "Live" & Day != 0) %>%
            droplevels() %>%
            select(Day, Ratio, Species, Experiment) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

meta_live_k_nc_samples <- meta_k_nc_model$sample(
          data = deco %>%
            filter(Condition == "Live" & Day != 0) %>%
            droplevels() %>%
            select(Day, Ratio, Species, Experiment) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

meta_dead_k_c_samples <- meta_k_c_model$sample(
          data = deco %>%
            filter(Condition == "Dead" & Day != 0) %>%
            droplevels() %>%
            select(Day, Ratio, Species, Experiment) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

meta_dead_k_nc_samples <- meta_k_nc_model$sample(
          data = deco %>%
            filter(Condition == "Dead" & Day != 0) %>%
            droplevels() %>%
            select(Day, Ratio, Species, Experiment) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

# 2.2.3 Model checks ####
# Rhat
meta_live_k_c_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.000175.

meta_live_k_nc_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.000112.

meta_dead_k_c_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.000190.

meta_dead_k_nc_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000661.

# Models are similar in both cases.

# Chains
meta_live_k_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
meta_live_k_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Little difference between chains

meta_dead_k_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
meta_dead_k_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Non-centred chains are better

# Pairs
meta_live_k_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k_mu", "log_k_s[1]", "log_k_e[1]"))
meta_live_k_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("sigma", "log_k_sigma_s", "log_k_sigma_e"))

meta_live_k_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k_mu", "log_k_s[1]", "log_k_e[1]"))
meta_live_k_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("sigma", "log_k_sigma_s", "log_k_sigma_e"))

meta_dead_k_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k_mu", "log_k_s[1]", "log_k_e[1]"))
meta_dead_k_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("sigma", "log_k_sigma_s", "log_k_sigma_e"))

meta_dead_k_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k_mu", "log_k_s[1]", "log_k_e[1]"))
meta_dead_k_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("sigma", "log_k_sigma_s", "log_k_sigma_e"))

# Some non-identifiability between species and experiment.
# Generally not much difference between centred and non-centred.

# 2.2.4 Prior-posterior comparison ####
meta_live_prior <- prior_samples(
  model = meta_k_nc_model,
  data = deco %>%
    filter(Condition == "Live" & Day != 0) %>%
    droplevels() %>%
    select(Day, Ratio, Species, Experiment) %>%
    compose_data()
  )

meta_dead_prior <- prior_samples(
  model = meta_k_nc_model,
  data = deco %>%
    filter(Condition == "Dead" & Day != 0) %>%
    droplevels() %>%
    select(Day, Ratio, Species, Experiment) %>%
    compose_data()
)

# Live centred species
meta_live_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_live_k_c_samples,
    group = deco %>%
            filter(Condition == "Live" & Day != 0) %>%
            droplevels() %>%
            select(Species),
    parameters = c("log_k_mu", "log_k_sigma_s", 
                   "log_k_s[Species]", "sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Species")

# Live centred experiments
meta_live_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_live_k_c_samples,
    group = deco %>%
            filter(Condition == "Live" & Day != 0) %>%
            droplevels() %>%
            select(Experiment),
    parameters = c("log_k_mu", "log_k_sigma_e", "log_k_e[Experiment]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Experiment")

# Live non-centred species
meta_live_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_live_k_nc_samples,
    group = deco %>%
            filter(Condition == "Live" & Day != 0) %>%
            droplevels() %>%
            select(Species),
    parameters = c("log_k_mu", "log_k_sigma_s", 
                   "log_k_s[Species]", "sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Species")

# Live non-centred experiments
meta_live_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_live_k_nc_samples,
    group = deco %>%
      filter(Condition == "Live" & Day != 0) %>%
      droplevels() %>%
      select(Experiment),
    parameters = c("log_k_mu", "log_k_sigma_e", "log_k_e[Experiment]"),
    format = "long"
  ) %>%
  prior_posterior_plot(group_name = "Experiment")
# Posteriors look similar

# Dead centred species and experiments
meta_dead_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_dead_k_c_samples,
    group = deco %>%
            filter(Condition == "Dead" & Day != 0) %>%
            droplevels() %>%
            select(Species, Experiment),
    parameters = c("log_k_mu", "log_k_sigma_s", "log_k_s[Species]", 
                   "log_k_sigma_e", "log_k_e[Experiment]", "sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Species",
    second_group_name = "Experiment"
  )

# Dead non-centred species and experiments
meta_dead_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_dead_k_nc_samples,
    group = deco %>%
            filter(Condition == "Dead" & Day != 0) %>%
            droplevels() %>%
            select(Species, Experiment),
    parameters = c("log_k_mu", "log_k_sigma_s", "log_k_s[Species]", 
                   "log_k_sigma_e", "log_k_e[Experiment]", "sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Species",
    second_group_name = "Experiment"
  )
# Posteriors look similar

# Models are very similar in all respects. Choose non-centred model 
# as optimal in both cases due to marginally better rhat and chains.

# 2.2.5 Prediction ####
# Global parameters
meta_live_k_prior_posterior_global <- meta_live_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_live_k_nc_samples,
    parameters = c("log_k_mu", "log_k_sigma_s",
                   "log_k_sigma_e", "sigma"),
    format = "short"
  ) %>%
  mutate( # Calculate k for new species and experiments
    k = exp(
      rnorm( n() , log_k_mu , log_k_sigma_s ) +
        rnorm( n() , 0 , log_k_sigma_e )
    )
  ) %>%
  select(starts_with("."), distribution,
         k, log_k_sigma_s, log_k_sigma_e, sigma) %T>%
  print()

meta_dead_k_prior_posterior_global <- meta_dead_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_dead_k_nc_samples,
    parameters = c("log_k_mu", "log_k_sigma_s",
                   "log_k_sigma_e", "sigma"),
    format = "short"
  ) %>%
  mutate(
    k = exp(
      rnorm( n() , log_k_mu , log_k_sigma_s ) +
        rnorm( n() , 0 , log_k_sigma_e )
    )
  ) %>%
  select(starts_with("."), distribution, 
         k, log_k_sigma_s, log_k_sigma_e, sigma) %T>%
  print()

meta_k_prior_posterior_global <- meta_live_k_prior_posterior_global %>%
  mutate(Condition = "Live" %>% fct()) %>%
  bind_rows(
    meta_dead_k_prior_posterior_global %>%
      mutate(Condition = "Dead" %>% fct())
  ) %>% # Priors for both conditions are identical
  filter(!(Condition == "Dead" & distribution == "prior")) %>%
  mutate(
    # Embed prior in condition
    Condition = if_else(
      distribution == "prior", "Prior", Condition
    ) %>% fct()
  ) %>%
  select(-distribution) %T>%
  print()

meta_k_prior_posterior_global %>%
  pivot_longer(cols = -c(starts_with("."), Condition),
               names_to = "parameter") %>%
  group_by(Condition, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n()) %>%
  print()

# Species parameters
meta_live_k_prior_posterior_species <- meta_live_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_live_k_nc_samples,
    group = deco %>%
            filter(Condition == "Live" & Day != 0) %>%
            droplevels() %>%
            select(Species),
    parameters = c("log_k_s[Species]", "log_k_sigma_e"),
    format = "short"
  ) %>%
  mutate( # Calculate k for new experiments on observed species
    k = exp(
      rnorm( n() , log_k_s , log_k_sigma_e )
    )
  ) %>%
  filter(Species == "Macrocystis pyrifera" & distribution == "prior" |
           distribution == "posterior") %>% # Remove redundant priors
  mutate( # Embed prior in species
    Species = if_else(
      distribution == "prior", "Prior", Species
    ) %>% fct()
  ) %>%
  select(starts_with("."), Species, k) %T>%
  print()

meta_dead_k_prior_posterior_species <- meta_dead_prior %>% 
  prior_posterior_draws(
    posterior_samples = meta_dead_k_nc_samples,
    group = deco %>%
            filter(Condition == "Dead" & Day != 0) %>%
            droplevels() %>%
            select(Species),
    parameters = c("log_k_s[Species]", "log_k_sigma_e"),
    format = "short"
  ) %>%
  mutate(
    k = exp(
      rnorm( n() , log_k_s , log_k_sigma_e )
    )
  ) %>%
  filter(Species == "Nereocystis luetkeana" & distribution == "prior" |
           distribution == "posterior") %>%
  mutate(
    Species = if_else(
      distribution == "prior", "Prior", Species
    ) %>% fct()
  ) %>%
  select(starts_with("."), Species, k) %T>%
  print()

meta_k_prior_posterior_species <- meta_live_k_prior_posterior_species %>%
  mutate(Condition = "Live" %>% fct()) %>%
  bind_rows(
    meta_dead_k_prior_posterior_species %>%
      mutate(Condition = "Dead" %>% fct())
  ) %>%
  filter(!(Condition == "Dead" & Species == "Prior")) %>%
  mutate(
    Condition = if_else(
      Species == "Prior", "Prior", Condition
    ) %>% fct()
  ) %T>%
  print()

meta_k_prior_posterior_species %>%
  pivot_longer(cols = -c(starts_with("."), Species, Condition),
               names_to = "parameter") %>%
  group_by(Species, Condition, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n()) %>%
  print(n = 100)

# Experiment parameters
meta_live_k_prior_posterior_experiment <- deco %>%
  filter(Condition == "Live" & Day != 0) %>%
  droplevels() %>%
  distinct(Species, Experiment) %>% # Get species-experiment pairs from data
  left_join( # Join the species samples by species
    meta_live_prior %>% 
      prior_posterior_draws(
        posterior_samples = meta_live_k_nc_samples,
        group = deco %>%
          filter(Condition == "Live" & Day != 0) %>%
          droplevels() %>%
          select(Species),
        parameters = c("log_k_s[Species]"),
        format = "short"
      ),
    by = "Species",
    relationship = "many-to-many"
  ) %>%
  left_join( # Join the experiment samples by experiment and sample ID
    meta_live_prior %>% 
      prior_posterior_draws(
        posterior_samples = meta_live_k_nc_samples,
        group = deco %>%
          filter(Condition == "Live" & Day != 0) %>%
          droplevels() %>%
          select(Experiment),
        parameters = c("log_k_e[Experiment]"),
        format = "short"
      ),
    by = c("Experiment", ".chain", ".iteration", ".draw", "distribution"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    k = exp( log_k_s + log_k_e )
  ) %>%
  filter(Species == "Macrocystis pyrifera" & Experiment == 1 & # Remove redundant priors
           distribution == "prior" | distribution == "posterior") %>% 
  mutate( # Embed prior in species and experiment
    Species = if_else(
      distribution == "prior", "Prior", Species
    ) %>% fct(),
    Experiment = if_else(
      distribution == "prior", "Prior", Experiment
    ) %>% fct()
  ) %>%
  select(starts_with("."), Species, Experiment, k) %T>%
  print()
  
meta_dead_k_prior_posterior_experiment <- deco %>%
  filter(Condition == "Dead" & Day != 0) %>%
  droplevels() %>%
  distinct(Species, Experiment) %>%
  left_join(
    meta_dead_prior %>% 
      prior_posterior_draws(
        posterior_samples = meta_dead_k_nc_samples,
        group = deco %>%
          filter(Condition == "Dead" & Day != 0) %>%
          droplevels() %>%
          select(Species),
        parameters = c("log_k_s[Species]"),
        format = "short"
      ),
    by = "Species",
    relationship = "many-to-many"
  ) %>%
  left_join(
    meta_dead_prior %>% 
      prior_posterior_draws(
        posterior_samples = meta_dead_k_nc_samples,
        group = deco %>%
          filter(Condition == "Dead" & Day != 0) %>%
          droplevels() %>%
          select(Experiment),
        parameters = c("log_k_e[Experiment]"),
        format = "short"
      ),
    by = c("Experiment", ".chain", ".iteration", ".draw", "distribution"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    k = exp( log_k_s + log_k_e )
  ) %>%
  filter(Species == "Nereocystis luetkeana" & Experiment == 14 &
           distribution == "prior" | distribution == "posterior") %>% 
  mutate(
    Species = if_else(
      distribution == "prior", "Prior", Species
    ) %>% fct(),
    Experiment = if_else(
      distribution == "prior", "Prior", Experiment
    ) %>% fct()
  ) %>%
  select(starts_with("."), Species, Experiment, k) %T>%
  print()

meta_k_prior_posterior_experiment <- meta_live_k_prior_posterior_experiment %>%
  mutate(Condition = "Live" %>% fct()) %>%
  bind_rows(
    meta_dead_k_prior_posterior_experiment %>%
      mutate(Condition = "Dead" %>% fct())
  ) %>%
  filter(!(Condition == "Dead" & Species == "Prior")) %>%
  mutate(
    Condition = if_else(
      Species == "Prior", "Prior", Condition
    ) %>% fct()
  ) %T>%
  print()

meta_k_prior_posterior_experiment %>%
  pivot_longer(cols = -c(starts_with("."), Species, Experiment, Condition),
               names_to = "parameter") %>%
  group_by(Species, Experiment, Condition, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n()) %>%
  print(n = 200)

# Save progress and clean up
meta_k_prior_posterior_global %>%
  write_rds(here("Decomposition", "RDS", "meta_k_prior_posterior_global.rds"))
meta_k_prior_posterior_species %>%
  write_rds(here("Decomposition", "RDS", "meta_k_prior_posterior_species.rds"))
meta_k_prior_posterior_experiment %>%
  write_rds(here("Decomposition", "RDS", "meta_k_prior_posterior_experiment.rds"))

rm(meta_k_c_model, meta_k_nc_model, meta_live_prior, meta_dead_prior,
   meta_live_k_c_samples, meta_live_k_nc_samples, meta_dead_k_c_samples, meta_dead_k_nc_samples,
   meta_live_k_prior_posterior_global, meta_dead_k_prior_posterior_global,
   meta_live_k_prior_posterior_species, meta_dead_k_prior_posterior_species,
   meta_live_k_prior_posterior_experiment, meta_dead_k_prior_posterior_experiment)

# Calculate global rounded values
meta_k_parameters_global <- meta_k_prior_posterior_global %>%
  filter(Condition != "Prior") %>%
  select(!starts_with(".")) %>%
  group_by(Condition) %>%
  summarise(
    across( everything(), list(mean = mean, sd = sd, median = median) ),
    n = n()
  ) %>%
  ungroup() %>%
  mutate( # Note I am converting dimensionless rates to %
    k_mean_rounded = if_else(
      k_mean < 1, signif(k_mean*100, 2), signif(k_mean*100, 3)
    ),
    k_sd_rounded = if_else(
      k_sd < 1, signif(k_sd*100, 2), signif(k_sd*100, 3)
    ),
    k_median_rounded = if_else(
      k_median < 1, signif(k_median*100, 2), signif(k_median*100, 3)
    ),
    k = glue("{k_mean_rounded} ± {k_sd_rounded} ({k_median_rounded})"),
    sigma = glue("{signif(sigma_mean, 2)} ± {signif(sigma_sd, 2)}")
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# Calculate species rounded values for supplementary table
meta_k_parameters_species <- meta_k_prior_posterior_species %>%
  select(!starts_with(".")) %>%
  filter(Species != "Prior") %>%
  group_by(Species, Condition) %>%
  summarise(
    across( everything(), list(mean = mean, sd = sd, median = median) ),
    n = n()
  ) %>%
  ungroup() %>%
  mutate( # Note I am converting dimensionless rates to %
    k_mean_rounded = if_else(
      k_mean < 1, signif(k_mean*100, 2), signif(k_mean*100, 3)
    ),
    k_sd_rounded = if_else(
      k_sd < 1, signif(k_sd*100, 2), signif(k_sd*100, 3)
    ),
    k_median_rounded = if_else(
      k_median < 1, signif(k_median*100, 2), signif(k_median*100, 3)
    ),
    k = glue("{k_mean_rounded} ± {k_sd_rounded} ({k_median_rounded})")
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# Calculate experiment rounded values for text
meta_k_parameters_experiment <- meta_k_prior_posterior_experiment %>%
  select(!starts_with(".")) %>%
  filter(Species != "Prior") %>%
  group_by(Species, Experiment, Condition) %>%
  summarise(
    across( everything(), list(mean = mean, sd = sd, median = median) ),
    n = n()
  ) %>%
  ungroup() %>%
  mutate( # Note I am converting dimensionless rates to %
    k_mean_rounded = if_else(
      k_mean < 1, signif(k_mean*100, 2), signif(k_mean*100, 3)
    ),
    k_sd_rounded = if_else(
      k_sd < 1, signif(k_sd*100, 2), signif(k_sd*100, 3)
    ),
    k_median_rounded = if_else(
      k_median < 1, signif(k_median*100, 2), signif(k_median*100, 3)
    ),
    k = glue("{k_mean_rounded} ± {k_sd_rounded} ({k_median_rounded})")
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print(n = 185)

# Calculate contrasts
meta_contrast <- meta_k_prior_posterior_global %>%
  filter(Condition != "Prior") %>%
  select(!contains("sigma")) %>%
  pivot_wider(names_from = Condition, values_from = k) %>%
  mutate(difference = Dead - Live,
         ratio = Dead / Live,
         log_ratio = log10(ratio)) %T>%
  print()

meta_k_sigma <- meta_k_prior_posterior_global %>%
  filter(Condition != "Prior") %>%
  select(Condition, starts_with("."), contains("k_sigma")) %>%
  pivot_longer(cols = contains("sigma"),
               names_to = "Parameter") %>%
  pivot_wider(names_from = Condition, values_from = value) %>%
  mutate(difference = Live - Dead,
         ratio = Live / Dead,
         log_ratio = log10(ratio)) %T>%
  print()

# Summarise
meta_contrast_summary <- meta_contrast %>%
  select(!starts_with(".")) %>%
  summarise(
    across( everything(), list(mean = mean, sd = sd, median = median) ),
    n = n(),
    P = mean( difference > 0 ) %>% signif(2)
  ) %>%
  mutate( # Note I am converting dimensionless rates to %
    Live_mean_rounded = if_else(
      Live_mean < 1, signif(Live_mean*100, 2), signif(Live_mean*100, 3)
    ),
    Live_sd_rounded = if_else(
      Live_sd < 1, signif(Live_sd*100, 2), signif(Live_sd*100, 3)
    ),
    Live_median_rounded = if_else(
      Live_median < 1, signif(Live_median*100, 2), signif(Live_median*100, 3)
    ),
    Live = glue("{Live_mean_rounded} ± {Live_sd_rounded} ({Live_median_rounded})"),
    Dead_mean_rounded = if_else(
      Dead_mean < 1, signif(Dead_mean*100, 2), signif(Dead_mean*100, 3)
    ),
    Dead_sd_rounded = if_else(
      Dead_sd < 1, signif(Dead_sd*100, 2), signif(Dead_sd*100, 3)
    ),
    Dead_median_rounded = if_else(
      Dead_median < 1, signif(Dead_median*100, 2), signif(Dead_median*100, 3)
    ),
    Dead = glue("{Dead_mean_rounded} ± {Dead_sd_rounded} ({Dead_median_rounded})"),
    difference_mean_rounded = if_else(
      difference_mean < 1, signif(difference_mean*100, 2), signif(difference_mean*100, 3)
    ),
    difference_sd_rounded = if_else(
      difference_sd < 1, signif(difference_sd*100, 2), signif(difference_sd*100, 3)
    ),
    difference_median_rounded = if_else(
      difference_median < 1, signif(difference_median*100, 2), signif(difference_median*100, 3)
    ),
    difference = glue("{difference_mean_rounded} ± {difference_sd_rounded} ({difference_median_rounded})"),
    ratio_mean_rounded = if_else(
      ratio_mean < 100, signif(ratio_mean, 2), signif(ratio_mean, 3)
    ),
    ratio_sd_rounded = if_else(
      ratio_sd < 100, signif(ratio_sd, 2), signif(ratio_sd, 3)
    ),
    ratio_median_rounded = if_else(
      ratio_median < 100, signif(ratio_median, 2), signif(ratio_median, 3)
    ),
    ratio = glue("{ratio_mean_rounded} ± {ratio_sd_rounded} ({ratio_median_rounded})"),
    log_ratio = glue("{signif(log_ratio_mean, 2)} ± {signif(log_ratio_sd, 2)}")
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

meta_k_sigma_summary <- meta_k_sigma %>%
  select(!starts_with(".")) %>%
  group_by(Parameter) %>%
  summarise(
    across( everything(), list(mean = mean, sd = sd) ),
    n = n(),
    P = mean( difference > 0 ) %>% signif(2)
  ) %>%
  ungroup() %>%
  mutate(
    Live = glue("{signif(Live_mean, 2)} ± {signif(Live_sd, 2)}"),
    Dead = glue("{signif(Dead_mean, 2)} ± {signif(Dead_sd, 2)}"),
    difference = glue("{signif(difference_mean, 2)} ± {signif(difference_sd, 2)}"),
    log_ratio = glue("{signif(log_ratio_mean, 2)} ± {signif(log_ratio_sd, 2)}")
  ) %>%
  select(!(contains("mean") | contains("sd"))) %T>%
  print()

# Predict across predictor range
# Predict for new experiments on new species
meta_k_prediction_global <- meta_k_prior_posterior_global %>%
  spread_continuous(data = deco, 
                    predictor_name = "Day", 
                    group_name = "Condition") %>% 
  mutate(
    r_mu = exp( -k * Day ),
    r = rnorm( n() , r_mu , sigma )
  ) %T>%
  print()

# Summarise predictions
meta_k_prediction_global_summary <- meta_k_prediction_global %>%
  group_by(Day, Condition) %>%
  median_qi(r_mu, r, .width = c(.5, .8, .9)) %T>%
  print()

rm(meta_k_prediction_global)

# Predict for new experiments on existing species
meta_k_prediction_species <- meta_k_prior_posterior_species %>%
  spread_continuous(data = deco, 
                    predictor_name = "Day",
                    group_name = "Species",
                    second_group_name = "Condition") %>%
  mutate(
    r_mu = exp( -k * Day )
  ) %T>%
  print()

# Summarise predictions
meta_k_prediction_species_summary <- meta_k_prediction_species %>%
  group_by(Day, Species, Condition) %>%
  summarise(
    across(
      r_mu,
      list(median = median, 
           lower_0.5 = ~ qi(.x, .width = .5)[1],
           upper_0.5 = ~ qi(.x, .width = .5)[2],
           lower_0.8 = ~ qi(.x, .width = .8)[1],
           upper_0.8 = ~ qi(.x, .width = .8)[2],
           lower_0.9 = ~ qi(.x, .width = .9)[1],
           upper_0.9 = ~ qi(.x, .width = .9)[2]),
      .names = "{.col}.{.fn}"
    )
  ) %>%
  ungroup() %>%
  rename(r_mu = r_mu.median) %>%
  pivot_longer(cols = contains("lower") | contains("upper")) %>%
  separate(col = name, into = c("name", ".width"), sep = "_(?=[^_]*$)") %>%
  mutate(.width = .width %>% as.numeric()) %>%
  pivot_wider(names_from = name, values_from = value) %T>%
  print()

rm(meta_k_prediction_species)

# Predict for existing experiments
meta_k_prediction_experiment <- meta_k_prior_posterior_experiment %>% 
  group_by(Experiment, Condition, .chain) %>%
  slice_sample(n = 1500) %>% # Reduce the number of samples per group and chain
  ungroup() %>%
  spread_continuous(data = deco, 
                    predictor_name = "Day",
                    group_name = "Experiment",
                    second_group_name = "Condition") %>%
  mutate(
    r_mu = exp( -k * Day )
  ) %T>%
  print()

# Summarise predictions
meta_k_prediction_experiment_summary <- meta_k_prediction_experiment %>%
  group_by(Day, Species, Experiment, Condition) %>%
  summarise(
    across(
      r_mu,
      list(median = median, 
           lower_0.5 = ~ qi(.x, .width = .5)[1],
           upper_0.5 = ~ qi(.x, .width = .5)[2],
           lower_0.8 = ~ qi(.x, .width = .8)[1],
           upper_0.8 = ~ qi(.x, .width = .8)[2],
           lower_0.9 = ~ qi(.x, .width = .9)[1],
           upper_0.9 = ~ qi(.x, .width = .9)[2]),
      .names = "{.col}.{.fn}"
    )
  ) %>%
  ungroup() %>%
  rename(r_mu = r_mu.median) %>%
  pivot_longer(cols = contains("lower") | contains("upper")) %>%
  separate(col = name, into = c("name", ".width"), sep = "_(?=[^_]*$)") %>%
  mutate(.width = .width %>% as.numeric()) %>%
  pivot_wider(names_from = name, values_from = value) %T>%
  print()

rm(meta_k_prediction_experiment)

# Save progress
meta_k_prediction_global_summary %>%
  write_rds(here("Decomposition", "RDS", "meta_k_prediction_global.rds"))
meta_k_prediction_species_summary %>%
  write_rds(here("Decomposition", "RDS", "meta_k_prediction_species.rds"))
meta_k_prediction_experiment_summary %>%
  write_rds(here("Decomposition", "RDS", "meta_k_prediction_experiment.rds"))

# 3. Figures ####
# 3.1 Combine predictions ####
prediction <- bind_rows(
    meta_prediction_global_summary %>%
      select(starts_with("."), Day, Condition, 
             starts_with("r_mu"), starts_with("k")) %>%
      mutate(Predictor = "Species" %>% fct(),
             Group = "Global" %>% fct()),
    meta_prediction_global_summary %>%
      select(starts_with("."), Day, Condition, 
             starts_with("r_mu"), starts_with("k")) %>%
      mutate(Predictor = "Experiments" %>% fct(),
             Group = "Global" %>% fct()),
    meta_prediction_species_summary %>%
      mutate(Predictor = "Species" %>% fct()) %>%
      rename(Group = Species),
    meta_prediction_experiment_summary %>%
      mutate(Predictor = "Experiments" %>% fct()) %>%
      rename(Group = Experiment) %>%
      select(!Species)
  ) %>%
  mutate( # Make sure Species is the first panel
    Predictor = Predictor %>% fct_relevel("Species"),
    # Global needs to be highlighted in the plot
    Highlight = Group == "Global",
    # Make sure Global is plotted last (on top)
    Group = Group %>% fct_reorder(Highlight)
  ) %T>%
  print()

k <- bind_rows(
    meta_k_prior_posterior_global %>%
      mutate(Predictor = "Species" %>% fct(),
             Group = "Global" %>% fct()) %>%
      select(!contains("sigma")),
    meta_k_prior_posterior_global %>%
      mutate(Predictor = "Experiments" %>% fct(),
             Group = "Global" %>% fct()) %>%
      select(!contains("sigma")),
    meta_k_prior_posterior_species %>%
      mutate(Predictor = "Species" %>% fct()) %>%
      rename(Group = Species),
    meta_k_prior_posterior_experiment %>%
      mutate(Predictor = "Experiments" %>% fct()) %>%
      rename(Group = Experiment) %>%
      select(!Species)
  ) %T>%
  print()

# Calculate densities manually
k_dens <- k %>%
  group_by(Predictor, Condition, Group) %>%
  reframe(x = c(0, density(k, n = 2^10, from = 0, to = 0.6, bw = 0.6 * 0.02)$x, 0.6),
          y = c(0, density(k, n = 2^10, from = 0, to = 0.6, bw = 0.6 * 0.02)$y, 0)) %>%
  group_by(Predictor, Condition, Group) %>% # Standardise area with Riemann sum (avoid manually added x[1]).
  mutate(y = y * 0.1 / ( sum(y) * ( x[3] - x[2] ) )) %>%
  ungroup() %>%
  mutate( # Make sure Species is the first panel
    Predictor = Predictor %>% fct_relevel("Species"),
    # Global needs to be highlighted in the plot
    Highlight = Group == "Global"
  ) %>%
  group_by(Predictor) %>%
  arrange(Highlight, .by_group = TRUE) %>%
  ungroup() %>%
  mutate( # A nested group is required for densities
    Group_nested = str_c(Predictor, Condition, Group, sep = "_") %>%
      fct_inorder() # This relevels according to specified arrange order
  ) %T>%
  print()

# 3.2 Figure 3 ####
# 3.2.1 Figure 3a ####
# Define custom theme
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = margin(0.2, 0.5, 0.2, 0.2, unit = "cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 12, hjust = 0),
                 axis.text = element_text(size = 10, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black", lineend = "square"),
                 legend.key = element_blank(),
                 legend.key.width = unit(.25, "cm"),
                 legend.key.height = unit(.45, "cm"),
                 legend.key.spacing.x = unit(.5, "cm"),
                 legend.key.spacing.y = unit(.05, "cm"),
                 legend.background = element_blank(),
                 legend.position = "top",
                 legend.justification = 0,
                 legend.text = element_text(size = 12, hjust = 0),
                 legend.title = element_blank(),
                 legend.margin = margin(0, 0, 0, 0, unit = "cm"),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 12, hjust = 0),
                 panel.spacing.x = unit(1, "cm"),
                 panel.spacing.y = unit(0.6, "cm"),
                 text = element_text(family = "Futura"))

Fig_3a <- prediction %>%
  filter(Condition != "Prior") %>%
  ggplot() +
    geom_hline(yintercept = c(0, 1)) +
    geom_point(data = deco %>% filter(Day != 0), 
               aes(Day, Ratio),
               colour = "#a29400", shape = 16, 
               alpha = 0.05, size = 1.2) +
    geom_line(aes(Day, r_mu, group = Group, 
                  colour = Highlight, alpha = Highlight)) +
    scale_colour_manual(values = c("#a29400", "black"), guide = "none") +
    scale_alpha_manual(values = c(0.3, 1), guide = "none") +
    scale_y_continuous(breaks = seq(0, 2, 1)) +
    labs(x = "Detrital age (days)",
         y = expression("Relative detrital mass ("*italic(m)*"/"*italic(m)[0]*")")) +
    coord_cartesian(xlim = c(0, 200), ylim = c(0, 2),
                    expand = F, clip = "off") +
    facet_grid2(Condition ~ Predictor,
                switch = "y",
                strip = strip_nested(text_y = element_text(angle = 0, hjust = 0, vjust = 1))) +
    mytheme +
    theme(plot.margin = margin(0, 0.5, 0, 0.2, unit = "cm"))

Fig_3a

# 3.3.2 Figure 3b ####
Fig_3b <- k_dens %>%
  filter(Condition != "Prior") %>%
  ggplot() +
    geom_line(aes(x, if_else(Condition == "Live", y + 3.3, y),
                  group = Group_nested, alpha = Highlight,
                  colour = Highlight)) +
    geom_text(
      data = tibble(
        Predictor = c("Species", "Experiments") %>% rep(each = 2) %>% fct(),
        y = c(3.3, 0) %>% rep(2),
        label = c("Live", "Dead", NA %>% rep(2))
      ),
      aes(x = -0.157, y = y, label = label),
      family = "Futura", size = 12, size.unit = "pt",
      hjust = 0, vjust = 0
    ) +
    scale_colour_manual(values = c("#a29400", "black"), guide = "none") +
    scale_alpha_manual(values = c(0.3, 1), guide = "none") +
    scale_x_continuous(breaks = seq(0, 0.6, by = 0.2),
                       labels = scales::label_number(accuracy = c(1, 0.1 %>% rep(3)))) +
    facet_grid(~ Predictor) +
    labs(x = expression("Exponential decay ("*italic(k)*", day"^-1*")")) +
    coord_cartesian(xlim = c(0, 0.6), expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          strip.text = element_blank())

Fig_3b
# Safely ignore warning, which is due to intentional NAs in geom_text.

# 3.3.3 Figure 3c ####
# Add labels to meta_contrast
meta_contrast %<>%
  mutate(label_Live = meta_contrast_summary %$% 
           ( P * 100 ) %>% str_c("%"),
         label_Dead = meta_contrast_summary %$%
           ( (1 - P) * 100 ) %>% str_c("%")) %T>%
  print()

require(ggridges)
require(geomtextpath)
Fig_3c <- meta_contrast %>%
  ggplot() +
    stat_density_ridges(aes(ratio, 0), fill = "#a29400",
                        colour = NA, n = 2^10, from = -2, to = 4,
                        bandwidth = 6*0.02, scale = 1) +
    geom_textdensity(aes(x = ratio, y = after_stat(density),
                         label = label_Live),
                     colour = "#a29400", family = "Futura",
                     size = 3.5, hjust = 0.6, vjust = 0,
                     n = 2^10, bw = 6*0.02, text_only = T) +
    geom_textdensity(aes(x = ratio, y = after_stat(density),
                         label = label_Dead),
                     colour = "#a29400", family = "Futura",
                     size = 3.5, hjust = 0.25, vjust = 0,
                     n = 2^10, bw = 6*0.02, text_only = T) +
    geom_vline(aes(xintercept = 10^0)) +
    scale_x_log10(limits = c(10^-2, 10^4), 
                  breaks = 10^(-2:4),
                  labels = scales::label_log(),
                  oob = scales::oob_keep) +
    labs(x = expression("Relative exponential decay ("*italic(k)["Dead"]*"/"*italic(k)["Live"]*")")) +
    coord_cartesian(expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          strip.text = element_blank())

Fig_3c

# 3.3.4 Figure 3d ####
Fig_3d <- meta_k_prior_posterior_global %>%
  filter(Condition != "Prior") %>%
  pivot_longer(cols = starts_with("log"),
               names_to = "Predictor",
               values_to = "log_k_sigma") %>%
  mutate(
    Predictor = if_else(
      Predictor %>% str_detect("s$"), 
      "Species", "Experiments"
    ) %>% fct_relevel("Species"),
    Condition = Condition %>% fct_relevel("Dead"),
    log_k_variance = log_k_sigma^2
  ) %>%
  ggplot() +
    stat_density_ridges(aes(log_k_sigma, y = Condition), fill = "#a29400",
                        colour = NA, n = 2^10, from = 0, to = 2.4,
                        bandwidth = 2.4*0.02, scale = 2, alpha = 0.7) +
    geom_text(
      data = tibble(
        Predictor = c("Species", "Experiments") %>% rep(each = 2) %>% fct(),
        Condition = c("Live", "Dead") %>% rep(2) %>% fct(),
        label = c("Live", "Dead", NA %>% rep(2))
      ),
      aes(x = -0.628, y = Condition, label = label),
      family = "Futura", size = 12, size.unit = "pt",
      hjust = 0, vjust = 0
    ) +
    scale_x_continuous(limits = c(0, 2.4), 
                       breaks = seq(0, 2.4, by = 0.8),
                       labels = scales::label_number(accuracy = c(1, 0.1 %>% rep(3))),
                       oob = scales::oob_keep) +
    facet_grid(~ Predictor) +
    labs(x = expression("Variation in exponential decay ("*italic(σ)["ln"*italic(k)]*")")) +
    coord_cartesian(expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          strip.text = element_blank())

Fig_3d

# 3.3.4 Combine panels ####
require(patchwork)
Fig_3 <- ( Fig_3a / Fig_3b / Fig_3c / Fig_3d ) +
  plot_layout(heights = c(1, 0.555, 0.1, 0.2))

Fig_3 %>%
  ggsave(filename = "Fig_3.pdf", path = "Figures",
         device = cairo_pdf, height = 20, width = 14, units = "cm")

# 4. Tables ####
# 4.1 Table 1 ####
Table_1.2 <- meta_contrast_summary %>%
  mutate(Species = "Meta-analysis" %>% fct()) %>%
  select(Species, Live, Dead, difference, log_ratio, P) %T>%
  print()

Table_1.2 %>%
  write_csv(here("Tables", "Table_1.2.csv"))

require(officer)
read_docx() %>%
  body_add_table(value = Table_1.2) %>%
  print(target = here("Tables", "Table_1.2.docx"))

# 4.2 Table S2 ####
Table_S2.2 <- meta_k_parameters_global %>% 
  select(Condition, k, sigma) %>%
  left_join(
    meta_parameters_global %>%
      select(Condition, alpha, mu, tau, epsilon, lambda, theta),
    by = "Condition"
  ) %T>%
  print()

Table_S2.2 %>%
  write_csv(here("Tables", "Table_S2.2.csv"))

read_docx() %>%
  body_add_table(value = Table_S2.2) %>%
  print(target = here("Tables", "Table_S2.2.docx"))

Table_S2.2_reduced <- meta_k_parameters_global %>% 
  select(Condition, k) %>%
  left_join(
    meta_parameters_global %>%
      select(Condition, alpha, mu, tau),
    by = "Condition"
  ) %T>%
  print()

Table_S2.2_reduced %>%
  write_csv(here("Tables", "Table_S2.2_reduced.csv"))

read_docx() %>%
  body_add_table(value = Table_S2.2_reduced) %>%
  print(target = here("Tables", "Table_S2.2_reduced.docx"))

# 4.3 Table S4 ####
Table_S4 <- meta_k_parameters_species %>% 
  select(Species, Condition, k) %>%
  left_join(
    meta_parameters_species %>%
      select(Species, Condition, alpha, mu, tau),
    by = c("Species", "Condition")
  ) %>%
  mutate(Species = Species %>% fct_relevel(sort)) %>%
  arrange(Species) %T>%
  print()

Table_S4 %>%
  write_csv(here("Tables", "Table_S4.csv"))

read_docx() %>%
  body_add_table(value = Table_S4) %>%
  print(target = here("Tables", "Table_S4.docx"))

# 4.4 Table S5 ####
Table_S5 <- meta_k_sigma_summary %>% 
  bind_rows(meta_sigma_summary) %>%
  mutate(
    Predictor = if_else(
      Parameter %>% str_detect("s$"), 
      "Species", "Experiments"
    ) %>% fct_relevel("Species"),
    Parameter = Parameter %>% 
      str_remove("_sigma.*") %>%
      fct()
  ) %>%
  select(Parameter, Predictor, Live, Dead, difference, log_ratio, P) %>%
  arrange(Parameter, Predictor) %T>%
  print()

Table_S5 %>%
  write_csv(here("Tables", "Table_S5.csv"))

read_docx() %>%
  body_add_table(value = Table_S5) %>%
  print(target = here("Tables", "Table_S5.docx"))

# 4.5 Extra data for text ####
meta_global <- meta_parameters_global %>%
  select(-n) %>%
  pivot_longer(cols = -Condition,
               names_to = "Parameter") %>%
  pivot_wider(names_from = Condition,
              values_from = value) %T>%
  print()

meta_global %>%
  write_csv(here("Tables", "meta_global.csv"))

read_docx() %>%
  body_add_table(value = meta_global) %>%
  print(target = here("Tables", "meta_global.docx"))

meta_experiment <- meta_k_parameters_experiment %>% 
  select(Experiment, Species, Condition, k) %>%
  left_join(
    meta_parameters_experiment %>%
      select(Experiment, Species, Condition, alpha, mu, tau),
    by = c("Species", "Experiment", "Condition")
  ) %>%
  mutate(Species = Species %>% fct_relevel(sort)) %>%
  arrange(Experiment, Species) %T>%
  print()

meta_experiment %>%
  write_csv(here("Tables", "meta_experiment.csv"))

read_docx() %>%
  body_add_table(value = meta_experiment) %>%
  print(target = here("Tables", "meta_experiment.docx"))

# 5. Carbon sequestration ####
# Filbee-Dexter et al. 2024 (doi: 10.1038/s41561-024-01449-7) use 
# decomposition to calculate carbon sequestration by retention
# below 200 m depth as a proportion of net primary production.
# The equation is 0.71*exp(-k*CRT), 0.71 is the proportion of 
# that is exported, k is the exponential decay constant per day 
# and CRT is the coastal residence time in days. The 0.25, 0.5 
# and 0.75 quantiles of CRT are given as 10, 75 and 441 days.
# It's best to model the exported proportion and CRT which is
# possible since the necessary data are provided in Supplementary 
# Data 1 (https://www.nature.com/articles/s41561-024-01449-7#Sec13).

# 5.1 Model export ####
# 5.1.1 Prepare data ####
export <- here("Decomposition", "Export.csv") %>% 
  read_csv() %>% # NB provided units are incorrect but calculation is correct
  mutate(Proportion = Total_detritus_g_C_m2_yr / Avg_ann_prod_kg_C_m2_y) %T>% 
  print()

export %$% mean(Proportion)

export %>%
  ggplot(aes(Proportion)) +
    geom_density() +
    theme_minimal()
# Proportions of 1 or more are impossible.

export %<>%
  mutate(Proportion = if_else(Proportion < 1, Proportion, 0.9999)) %T>%
  print()

export %$% mean(Proportion)

export %>%
  ggplot(aes(Proportion)) +
    geom_density() +
    theme_minimal()

# 5.1.2 Stan model ####
# I am using the reported mean proportion of 0.71 as my prior.
export_model <- here("Decomposition", "Stan", "export.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

export_samples <- export_model$sample(
          data = export %>%
            select(Proportion) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

# 5.1.3 Model checks ####
# Rhat
export_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000393.

# Chains
export_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Chains are good.

# 5.1.4 Prior-posterior comparison ####
export_prior <- prior_samples(
  model = export_model,
  data = export %>%
    select(Proportion) %>%
    compose_data()
  )

export_prior %>% 
  prior_posterior_draws(
    posterior_samples = export_samples,
    parameters = c("mu", "nu"),
    format = "long"
    ) %>%
  prior_posterior_plot()

# 5.1.5 Prediction ####
export_prior_posterior <- export_prior %>% 
  prior_posterior_draws(
    posterior_samples = export_samples,
    parameters = c("mu", "nu"),
    format = "short"
  ) %>%
  mutate(
    Proportion = rbeta( n() , mu * nu , (1 - mu) * nu )
  ) %>%
  select(starts_with("."), distribution, Proportion) %T>%
  print()

# 5.2 Model coastal residence time ####
# 5.2.1 Prepare data ####
CRT <- here("Decomposition", "CRT.csv") %>% 
  read_csv() %T>% 
  print()

CRT %>%
  ggplot(aes(CRT)) +
    geom_density() +
    scale_x_log10() +
    theme_minimal()
# Fine but some NAs. Take care to remove.

# 5.2.2 Stan model ####
# I am using the reported median export time of 75 days as my prior.
CRT_model <- here("Decomposition", "Stan", "CRT.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

CRT_samples <- CRT_model$sample(
          data = CRT %>%
            select(CRT) %>%
            drop_na() %>% # NB NAs have to be removed
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

# 5.2.3 Model checks ####
# Rhat
CRT_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000212.

# Chains
CRT_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Chains are good.

# 5.2.4 Prior-posterior comparison ####
CRT_prior <- prior_samples(
  model = CRT_model,
  data = CRT %>%
    select(CRT) %>%
    drop_na() %>%
    compose_data()
  )

CRT_prior %>% 
  prior_posterior_draws(
    posterior_samples = CRT_samples,
    parameters = c("mu", "sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot()

# 5.2.5 Prediction ####
CRT_prior_posterior <- CRT_prior %>% 
  prior_posterior_draws(
    posterior_samples = CRT_samples,
    parameters = c("mu", "sigma"),
    format = "short"
  ) %>%
  mutate(
    CRT = rlnorm( n() , mu , sigma )
  ) %>%
  select(starts_with("."), distribution, CRT) %T>%
  print()

# 5.3 Carbon sequestration potential ####
# 5.3.1 Merge posteriors ####
CSP <- meta_k_prior_posterior_global %>%
  filter(Condition != "Prior") %>%
  select(starts_with("."), Condition, k) %>%
  full_join(
    export_prior_posterior %>%
      filter(distribution == "posterior") %>%
      select(-distribution),
    by = c(".chain", ".iteration", ".draw")
  ) %>%
  full_join(
    CRT_prior_posterior %>%
      filter(distribution == "posterior") %>%
      select(-distribution),
    by = c(".chain", ".iteration", ".draw")
  ) %T>%
  print()

# 5.3.2 Prediction ####
CSP %<>%
  mutate(CSP = Proportion * exp( -k * CRT ),
         # I am also calculating in log space because of underflow to 0
         log_CSP = log(Proportion) -k * CRT) %>%
  select(-c(k, Proportion, CRT)) %T>%
  print()

# Contrast
CSP_contrast <- CSP %>%
  select(-log_CSP) %>%
  pivot_wider(names_from = Condition, values_from = CSP) %>%
  mutate(difference = Live - Dead,
         log_ratio = log10(Live / Dead)) %T>%
  print()

CSP_contrast %>%
  count(is.finite(log_ratio))
# Quite a lot of Infs or NaNs introduced. Calculate
# log_ratio as difference of logs instead.

CSP_contrast %<>%
  select(-log_ratio) %>%
  full_join(
    CSP %>%
      select(-CSP) %>%
      pivot_wider(names_from = Condition, values_from = log_CSP) %>%
      # Ratio as difference of natural logs, then converted to log10
      mutate(log_ratio = (Live - Dead) / log(10)) %>%
      select(starts_with("."), log_ratio),
    by = c(".chain", ".iteration", ".draw")
  ) %T>%
  print()

CSP_contrast %>%
  count(is.finite(log_ratio))
# All values are finite.

# Summarise
CSP_contrast_summary <- CSP_contrast %>%
  select(!starts_with(".")) %>%
  summarise(
    across(
      everything(), 
      list(
        mean = mean, sd = sd, median = median
      )
    ),
    n = n(), # Calculate based on log_ratio rather than difference because 
    # this avoids values that underflowed to 0
    P = mean( log_ratio > 0 ) %>% signif(2)
  ) %>%
  mutate(
    Live = glue(
      "{signif(Live_mean, 2)} ± {signif(Live_sd, 2)} ({signif(Live_median, 2)})"
    ),
    Dead = glue(
      "{signif(Dead_mean, 2)} ± {signif(Dead_sd, 2)} ({signif(Dead_median, 2)})"
    ),
    difference = glue(
      "{signif(difference_mean, 2)} ± {signif(difference_sd, 2)} ({signif(difference_median, 2)})"
    ),
    log_ratio = glue(
      "{signif(log_ratio_mean, 2)} ± {signif(log_ratio_sd, 4)} ({signif(log_ratio_median, 2)})"
    )
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()
# Filbee-Dexter et al. 2024 report average CSP as 15 ± 2 (3 to 38) %.
# This is much closer to th median live estimate of 23% than the median
# dead estimate of 0.00058%. The live estimate is over an order of magnitude
# greater and there is an 86% chance that it is greater.

# 5.4 Figure S4 ####
# 5.4.1 Figure S4a ####
Fig_S4a <- CSP %>%
  mutate(Condition = Condition %>% fct_relevel("Dead")) %>%
  ggplot() +
    stat_density_ridges(aes(CSP, Condition), fill = "#a29400",
                        colour = NA, n = 2^10, from = 0, to = 1,
                        bandwidth = 1*0.02, scale = 0.9) +
    scale_x_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, 0.2),
                       labels = scales::label_number(accuracy = c(1, rep(0.1, 4), 1)),
                       oob = scales::oob_keep) +
    labs(x = "Proportion exported below 200 metres") +
    coord_cartesian(expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = 12, vjust = 0, hjust = 0,
                                     margin = margin(r = 10)),
          axis.line.y = element_blank(),
          plot.margin = margin(0.2, 0.5, 0.2, 0, unit = "cm"))

Fig_S4a

# 5.4.2 Figure S4b ####
# Add labels to CSP_contrast
CSP_contrast %<>%
  mutate(label_Live = CSP_contrast_summary %$% 
           ( P * 100 ) %>% str_c("%"),
         label_Dead = CSP_contrast_summary %$%
           ( (1 - P) * 100 ) %>% str_c("%")) %T>%
  print()

Fig_S4b <- CSP_contrast %>%
  ggplot() +
    stat_density_ridges(aes(difference, 0), fill = "#a29400",
                        colour = NA, n = 2^10, from = -1, to = 1,
                        bandwidth = 2*0.02, scale = 1) +
    geom_textdensity(aes(x = difference, y = after_stat(density),
                         label = label_Live),
                     colour = "#a29400", family = "Futura",
                     size = 3.5, hjust = 0.8, vjust = 0,
                     n = 2^10, bw = 2*0.02, text_only = T) +
    geom_textdensity(aes(x = difference, y = after_stat(density),
                         label = label_Dead),
                     colour = "#a29400", family = "Futura",
                     size = 3.5, hjust = 0.4, vjust = 0,
                     n = 2^10, bw = 2*0.02, text_only = T) +
    geom_vline(aes(xintercept = 0)) +
    scale_x_continuous(limits = c(-1, 1),
                       breaks = seq(-1, 1, 0.5),
                       labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1, 1),
                                                     style_negative = "minus"),
                       oob = scales::oob_keep) +
    labs(x = "Difference (Live − Dead)") +
    coord_cartesian(expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          plot.margin = margin(0.2, 0.5, 0.2, 0, unit = "cm"))

Fig_S4b

# 5.4.3 Figure S4c ####
Fig_S4c <- CSP_contrast %>%
  ggplot() +
    stat_density_ridges(aes(log_ratio, 0), fill = "#a29400",
                        colour = NA, n = 2^10, from = -10, to = 10,
                        bandwidth = 20*0.02, scale = 1) +
    geom_textdensity(aes(x = log_ratio, y = after_stat(density),
                         label = label_Live),
                     colour = "#a29400", family = "Futura",
                     size = 3.5, hjust = 0.8, vjust = 0,
                     n = 2^10, bw = 20*0.02, text_only = T) +
    geom_textdensity(aes(x = log_ratio, y = after_stat(density),
                         label = label_Dead),
                     colour = "#a29400", family = "Futura",
                     size = 3.5, hjust = 0.4, vjust = 0,
                     n = 2^10, bw = 20*0.02, text_only = T) +
    geom_vline(aes(xintercept = 0)) +
    scale_x_continuous(limits = c(-10, 10),
                       breaks = seq(-10, 10, 5),
                       labels = scales::label_math(10^.x), # have to fix minus
                       oob = scales::oob_keep) +
    labs(x = "Ratio (Live / Dead)") +
    coord_cartesian(expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          plot.margin = margin(0.2, 0.5, 0.2, 0, unit = "cm"))

Fig_S4c

# 5.4.4 Combine panels ####
Fig_S4 <- ( Fig_S4a / Fig_S4b / Fig_S4c ) +
  plot_layout(heights = c(1, 0.5, 0.5))

Fig_S4 %>%
  ggsave(filename = "Fig_S4.pdf", path = "Figures",
         device = cairo_pdf, height = 15, width = 10, units = "cm")