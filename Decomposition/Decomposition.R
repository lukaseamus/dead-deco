#### Dead deco: kelp decomposition without physiology ####
#### Luka Seamus Wright                               ####

# 1. Prepare data ####
# 1.1 Load data ####
require(tidyverse)
require(magrittr)
require(here)

decouk <- here("Decomposition", "Decomposition_UK.csv") %>% read_csv(col_types = list("f", "f", "f")) %>%
  mutate(Day = Hour / 24,
         # Replace 0 with small constant within balance measurement error
         # because the model is undefined for y = 0
         Final = if_else(Final == 0, 1e-5, Final)) %T>%
  print()

decoau <- here("Decomposition", "Decomposition_AU.csv") %>% read_csv(col_types = list("f", "f", "f", "f")) %>%
  mutate(Deployment = Deployment %>% dmy(),
         Retrieval = Retrieval %>% dmy(),
         Day = Deployment %--% Retrieval / ddays(),
         Final = if_else(Final == 0, 1e-5, Final),
         Dry = if_else(Dry == 0, 1e-6, Dry)) %T>%
  print()

ratio <- here("Decomposition", "Ratio.csv") %>% read_csv(col_types = list("f", "f")) %>%
  mutate(Ratio = Dry / Fresh) %T>%
  print()

# 1.2 Mass ratio ####
# 1.2.1 Summary ####
ratio %>%
  group_by(Species) %>%
  summarise(Ratio_mean = mean(Ratio),
            Ratio_sd = sd(Ratio),
            n = n())
# Better use a predictive model

# 1.2.2 Prior simulation ####
ratio %>%
  ggplot(aes(Fresh, Dry)) +
    geom_point() +
    facet_wrap(~ Species, scales = "free") +
    theme_minimal()

# For Amphibolis griffithii the dry-wet mass ratio is typically around
# 0.27 (de los Santos et al. 2012, doi 10.3354/meps09757) to 
# 0.31 (Borum et al. 2016, doi 10.1111/pce.12658). I am not aware of
# accessible dry-wet mass ratios for Ecklonia radiata. The mass
# ratio is beta distributed because dry mass cannot exceed wet mass.

tibble(n = 1:1e3,
       beta_mu = rbeta( 1e3 , 0.29 * 30 , (1 - 0.29) * 30 ), 
       beta_nu = rgamma( 1e3 , 20^2 / 10^2 , 20 / 10^2 ), # low nu = large variation
       beta = rbeta( 1e3 , beta_mu * beta_nu , (1 - beta_mu) * beta_nu )) %>% # %$% hist(beta)
  expand_grid(Fresh = ratio %$% 
                seq(min(Fresh), max(Fresh), length.out = 50)) %>%
  mutate(Dry = beta * Fresh) %>%
  ggplot(aes(Fresh, Dry, group = n)) +
    geom_abline(slope = 1) +
    geom_hline(yintercept = ratio %$% 
                 range(Dry)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Likelihood should be a truncated normal to avoid negative predictions.

# 1.2.3 Stan model ####
require(cmdstanr)
ratio_model <- here("Decomposition", "Stan", "ratio.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
ratio_samples <- ratio_model$sample(
          data = ratio %>% 
            select(Species, Fresh, Dry) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4,
          adapt_delta = 0.99,
          max_treedepth = 15 # force sampler to slow down
        )

# 1.2.4 Model checks ####
# Rhat
ratio_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000416.

# Chains
require(bayesplot)
ratio_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Chains are good.

# Pairs
ratio_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("beta[1]", "beta_mu", "beta_nu", "cv[1]"))
ratio_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("beta[2]", "beta_mu", "beta_nu", "cv[2]"))
# No correlation.

# 1.2.5 Prior-posterior comparison ####
source("functions.R")
ratio_prior <- prior_samples(
  model = ratio_model,
  data = ratio %>% 
    select(Species, Fresh, Dry) %>%
    compose_data(),
  adapt_delta = 0.99,
  max_treedepth = 15
  # force sampler to slow down for smoother priors
  )

ratio_prior %>% 
  prior_posterior_draws(
    posterior_samples = ratio_samples,
    group = ratio %>% select(Species),
    parameters = c("beta_mu", "beta_nu", "cv_mu", "cv_theta",
                   "beta[Species]", "cv[Species]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Species", ridges = FALSE)
# Posteriors look good.

# 1.2.6 Prediction ####
# Priors and posteriors for hyperparameters
ratio_prior_posterior_hyper <- ratio_prior %>% 
  prior_posterior_draws(
    posterior_samples = ratio_samples,
    parameters = c("beta_mu", "beta_nu", "cv_mu", "cv_theta"),
    format = "short"
  ) %>% # Calculate predictions for new species, i.e. unobserved.
  mutate(beta = rbeta( n() , beta_mu * beta_nu , (1 - beta_mu) * beta_nu ),
         cv = rgamma( n() , cv_mu / cv_theta , 1 / cv_theta )) %T>%
  print()

# Priors and posteriors for species parameters
ratio_prior_posterior_species <- ratio_prior %>% 
  prior_posterior_draws(
    posterior_samples = ratio_samples,
    parameters = c("beta[Species]", "cv[Species]"),
    format = "short"
  ) %>% 
  # Since I want only one grouping variable, there is redundancy in distribution.
  filter(!(Species == "Ecklonia radiata" & 
             distribution == "prior")) %>% # Remove one redundant prior.
  mutate(Species = if_else(distribution == "prior", # Add Prior to Species.
                           "Prior", Species) %>% fct()) %>%
  select(-distribution) %T>%
  print()

# Combine species- and hyper-parameters
ratio_prior_posterior <- ratio_prior_posterior_hyper %>%
  # Prior is the same for species and hyper
  filter(distribution == "posterior") %>%
  mutate(Species = "Unobserved" %>% fct()) %>%
  select(-c(beta_mu, beta_nu, cv_mu, cv_theta, distribution)) %>%
  bind_rows(ratio_prior_posterior_species) %T>%
  print()

# Summarise parameters
ratio_parameters <- ratio_prior_posterior %>%
  group_by(Species) %>%
  summarise(beta_mean = mean(beta) %>% signif(2),
            beta_sd = sd(beta) %>% signif(2),
            cv_mean = mean(cv) %>% signif(2),
            cv_sd = sd(cv) %>% signif(2),
            n = n()) %T>%
  print()

# Predict across predictor range
require(extraDistr) # R doesn't have a native truncated normal
ratio_prediction <- ratio_prior_posterior %>%
  spread_continuous(data = ratio, 
                    predictor_name = "Fresh",
                    group_name = "Species") %>%
  mutate(mu = beta * Fresh,
         sigma = cv * mu,
         Dry = rtnorm( n() , mean = mu , sd = sigma , a = 0 )) %>%
  group_by(Species, Fresh) %>% # Summarise prediction distributions
  mean_qi(mu, Dry, .width = c(.5, .8, .9)) %T>%
  print()

# Save progress
ratio_prior_posterior %>%
  write_rds(here("Decomposition", "RDS", "ratio_prior_posterior.rds"))
ratio_prediction %>%
  write_rds(here("Decomposition", "RDS", "ratio_prediction.rds"))

# 1.2.7 Figure S1 ####
# Summarise parameters for annotation
require(glue)
ratio_annotation <- ratio_parameters %>%
  mutate(
    label_mu = glue(
      "italic(μ)*' = {beta_mean} ± {beta_sd} × '*italic(x)"
    ),
    label_sigma = glue(
      "italic(σ)*' = {cv_mean} ± {cv_sd} × '*italic(μ)"
    )
  ) %>%
  pivot_longer(cols = contains("label"),
               names_to = "parameter",
               values_to = "label",
               names_prefix = "label_") %T>%
  print()

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

# Plot
require(geomtextpath)
require(ggh4x)
Fig_S1 <- ggplot() + 
    # manually limit 1:1 line to species-specific range
    geom_textline(data = tibble(x = c(0, 18, 0, 2.4), y = c(0, 18, 0, 2.4),
                                Species = c("Ecklonia radiata", "Amphibolis griffithii") %>%
                                  rep(each = 2) %>% fct()), aes(x, y),
                  label = "1:1", family = "Futura", size = 3.5, hjust = 1) +
    geom_point(data = ratio,
               aes(Fresh, Dry, colour = Species),
               size = 2.4, alpha = 0.5, shape = 16) +
    # geom_ribbon(data = ratio_prediction %>%
    #               filter(Species == "Prior" &
    #                        .width == 0.9) %>%
    #               select(-Species),
    #             aes(Fresh, ymin = Dry.lower, ymax = Dry.upper),
    #             alpha = 0.1) +
    geom_line(data = ratio_prediction %>%
                filter(!Species %in% c("Unobserved", "Prior")),
              aes(Fresh, mu, colour = Species)) +
    # geom_ribbon(data = ratio_prediction %>%
    #               filter(!Species %in% c("Unobserved", "Prior")),
    #             aes(Fresh, ymin = mu.lower, ymax = mu.upper,
    #                 fill = Species, alpha = factor(.width))) +
    geom_ribbon(data = ratio_prediction %>%
                  filter(!Species %in% c("Unobserved", "Prior")),
                aes(Fresh, ymin = Dry.lower, ymax = Dry.upper,
                    fill = Species, alpha = factor(.width))) +
    geom_text(data = ratio_annotation %>%
                filter(!Species %in% c("Unobserved", "Prior")),
              aes(x = c(33, 33, 4.95, 4.95), 
                  y = c(2.5, 1, 0.339, 0.1335), 
                  label = label),
              family = "Futura", size = 10, size.unit = "pt", 
              hjust = 0, parse = TRUE, lineheight = 0.8) +
    scale_fill_manual(values = c("#c3b300", "#4a7518"), guide = "none") +
    scale_colour_manual(values = c("#c3b300", "#4a7518"), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    facet_wrap(~ Species, scales = "free") +
    facetted_pos_scales(
      x = list(
        Species == "Ecklonia radiata" ~ 
          scale_x_continuous(limits = c(0, 60), 
                             breaks = seq(0, 60, by = 20)),
        Species == "Amphibolis griffithii" ~ 
          scale_x_continuous(limits = c(0, 9),
                             breaks = seq(0, 9, by = 3))
        ),
      y = list(
        Species == "Ecklonia radiata" ~ 
          scale_y_continuous(limits = c(0, 18), 
                             breaks = seq(0, 18, by = 6)),
        Species == "Amphibolis griffithii" ~ 
          scale_y_continuous(limits = c(0, 2.4),
                             breaks = seq(0, 2.4, by = 0.6),
                             labels = scales::label_number(
                               accuracy = c(1, rep(0.1, 4))
                             ))
        )
      ) +
    labs(x = "Fresh mass (g)",
         y = "Dry mass (g)") +
    coord_cartesian(expand = FALSE, clip = "off") +
    mytheme +
    theme(strip.text = element_text(face = "italic"))

Fig_S1

Fig_S1 %>%
  ggsave(filename = "Fig_S1.pdf", path = "Figures",
         device = cairo_pdf, height = 10, width = 20, 
         units = "cm")
# For some reason annotation font turns bold, so needs to be edited.

# 1.3 Remaining mass ####
# 1.3.1 Dry mass ####
# Calculate with ratio model parameters
deco_dry <- decoau %>%
  left_join(
    ratio_prior_posterior %>%
      select(Species, beta, cv),
    by = "Species",
    relationship = "many-to-many"
  ) %>%
  droplevels() %>% # ratio_prior_posterior added levels
  mutate(
    mu = beta * Initial,
    sigma = cv * mu,
    Initial_dry = rtnorm( n() , mean = mu , sd = sigma , a = 0 ),
    Ratio = Dry / Initial_dry
  ) %T>%
  print()

# Summarise
deco_dry_summary <- deco_dry %>%
  group_by(ID, Species, Treatment, Tank, Deployment, Initial, 
           Retrieval, Final, Dry, Initials, Notes, Day) %>%
  summarise(
    beta_mean = mean(beta),
    beta_sd = sd(beta),
    Initial_dry_mean = mean(Initial_dry),
    Initial_dry_sd = sd(Initial_dry),
    Ratio_mean = mean(Ratio),
    Ratio_median = median(Ratio),
    Ratio_sd = sd(Ratio),
    Ratio_lwr = qi(Ratio, .width = 0.99)[1],
    Ratio_upr = qi(Ratio, .width = 0.99)[2],
    n = n()
  ) %>%
  ungroup() %T>%
  print()

# 1.3.2 Fresh mass ####
# Combine and calculate for both datasets 
deco_fresh <- decoau %>%
  select(ID, Species, Treatment, Initial, 
         Final, Day, Initials) %>%
  bind_rows(decouk %>% select(-Hour)) %>%
  mutate(Ratio = Final / Initial) %T>%
  print()

# 2. Model data ####
# 2.1 Dry mass macroalgal model ####
# 2.1.1 Visualisation ####
ggplot() +
  geom_violin(data = deco_dry,
              aes(Day, Ratio, group = ID)) +
  facet_grid(Treatment ~ Species) +
  theme_minimal()
# This won't do for visualisation.

ggplot() +
  geom_pointrange(data = deco_dry_summary,
                  aes(Day, Ratio_mean,
                      ymin = Ratio_mean - Ratio_sd,
                      ymax = Ratio_mean + Ratio_sd)) +
  facet_grid(Treatment ~ Species) +
  theme_minimal()
# But ± s.d. clearly doesn't represent the skew.

ggplot() +
  geom_pointrange(data = deco_dry_summary,
                  aes(Day, Ratio_mean,
                      ymin = Ratio_lwr,
                      ymax = Ratio_upr)) +
  facet_grid(Treatment ~ Species) +
  theme_minimal()
# Better

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

# I am taking the mean k value for Ecklonia radiata from 
# Simpkins et al. 2025 (doi: 10.1002/lno.70006) as my prior
# for tau: 0.06 d^-1. Amphibolis griffithii likely decomposes
# much slower so I'll make the variability large enough to
# include the smallest positive values. My prior for mu is
# based on estimates of the longevity of detrital photosynthesis
# for E. radiata from Wright et al. 2024 (doi: 10.1093/aob/mcad167)
# as well as observations on Amphibolis antarctica (unpublished).

# I need to enforce an increase in exponential decay over time
# due to the known decline in physiology. Therefore I am re-
# parameterising the model in terms of delta, mu and tau.
# delta = alpha + tau since alpha is the initial exponential  
# growth rate and tau is the final exponential decay rate.
# delta is chosen to allow for initial growth and decay, but
# favour decay, since half the experiments concern dead detritus.

tibble(n = 1:1e3,
       log_delta_mu = rnorm( 1e3 , log(0.03) , 0.3 ), # delta = alpha + tau
       log_mu_mu = rnorm( 1e3 , log(40) , 0.3 ),
       log_tau_mu = rnorm( 1e3 , log(0.06) , 0.3 ),
       log_delta_sigma = rtnorm( 1e3 , 0 , 0.3 , 0 ), # half-normal priors
       log_mu_sigma = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       log_tau_sigma = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       log_epsilon_mu = rnorm( 1e3 , log(4e4) , 0.3 ),
       log_lambda_mu = rnorm( 1e3 , log(0.1) , 0.3 ),
       log_theta_mu = rnorm( 1e3 , log(500) , 0.3 ),
       log_epsilon_sigma = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       log_lambda_sigma = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       log_theta_sigma = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       delta = exp( rnorm( 1e3 , log_delta_mu , log_delta_sigma ) ),
       mu = exp( rnorm( 1e3 , log_mu_mu , log_mu_sigma ) ),
       tau = exp( rnorm( 1e3 , log_tau_mu , log_tau_sigma ) ),
       epsilon = exp( rnorm( 1e3 , log_epsilon_mu , log_epsilon_sigma ) ),
       lambda = exp( rnorm( 1e3 , log_lambda_mu , log_lambda_sigma ) ),
       theta = exp( rnorm( 1e3 , log_theta_mu , log_theta_sigma ) ),
       alpha = delta - tau) %>%
  expand_grid(Day = deco_dry_summary %$% 
                seq(min(Day), max(Day), length.out = 100)) %>%
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
    facetted_pos_scales(
      y = list(
        parameter == "r" ~ scale_y_continuous(limits = c(0, 2)),
        parameter == "r_mu" ~ scale_y_continuous(limits = c(0, 2))
      )
    ) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 2.1.3 Stan models ####
dry_c_model <- here("Decomposition", "Stan", "dry_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

dry_nc_model <- here("Decomposition", "Stan", "dry_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

dry_c_samples <- dry_c_model$sample(
          data = deco_dry_summary %>%
            select(Day, Ratio_mean, Ratio_sd,
                   Species, Treatment) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

dry_nc_samples <- dry_nc_model$sample(
          data = deco_dry_summary %>%
            select(Day, Ratio_mean, Ratio_sd,
                   Species, Treatment) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

# 2.1.4 Model checks ####
# Rhat
dry_c_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 47% of rhat above 1.001. rhat = 1.00 ± 0.00108.

dry_nc_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 100% of rhat above 1.001. rhat = 1.07 ± 0.0410.

# Chains
dry_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Chains are ok 

dry_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Chains are bad

# Pairs
dry_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[1,1]", "log_mu[1,1]", "log_tau[1,1]"))
dry_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[1,2]", "log_mu[1,2]", "log_tau[1,2]"))

dry_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[2,1]", "log_mu[2,1]", "log_tau[2,1]"))
dry_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[2,2]", "log_mu[2,2]", "log_tau[2,2]"))

dry_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[1,1]", "log_mu[1,1]", "log_tau[1,1]"))
dry_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[1,2]", "log_mu[1,2]", "log_tau[1,2]"))

dry_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[2,1]", "log_mu[2,1]", "log_tau[2,1]"))
dry_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[2,2]", "log_mu[2,2]", "log_tau[2,2]"))
# Correlations are pretty bad for both models, but centred posteriors 
# look better

# 2.1.5 Prior-posterior comparison ####
dry_prior <- prior_samples(
  model = dry_c_model,
  data = deco_dry_summary %>%
    select(Day, Ratio_mean, Ratio_sd,
           Species, Treatment) %>%
    compose_data()
  )
# Too many divergences, cannot effectively sample centred hierarchical priors

dry_prior <- prior_samples(
  model = dry_nc_model,
  data = deco_dry_summary %>%
    select(Day, Ratio_mean, Ratio_sd,
           Species, Treatment) %>%
    compose_data()
  )
# Works smoothly -> use only non-centred priors, because they are the same

dry_prior %>% 
  prior_posterior_draws(
    posterior_samples = dry_c_samples,
    group = deco_dry_summary %>%
      select(Species, Treatment),
    parameters = c("log_delta_mu", "log_delta_sigma",
                   "log_delta[Species, Treatment]", 
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

dry_prior %>% 
  prior_posterior_draws(
    posterior_samples = dry_nc_samples,
    group = deco_dry_summary %>%
      select(Species, Treatment),
    parameters = c("log_delta_mu", "log_delta_sigma",
                   "log_delta[Species, Treatment]",
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
# Posteriors are smoother for the centred model
# I will re-run the centred model

# 2.1.6 Rerun optimal model ####
dry_c_samples <- dry_c_model$sample(
          data = deco_dry_summary %>%
            select(Day, Ratio_mean, Ratio_sd,
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

dry_c_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 4% of rhat above 1.001. rhat = 1.00 ± 0.000271.

dry_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Chains are very good

dry_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[1,1]", "log_mu[1,1]", "log_tau[1,1]"))
dry_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[1,2]", "log_mu[1,2]", "log_tau[1,2]"))

dry_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[2,1]", "log_mu[2,1]", "log_tau[2,1]"))
dry_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[2,2]", "log_mu[2,2]", "log_tau[2,2]"))
# Very strong correlation due to non-identifiability between delta and
# tau for seagrasses (the effect of detrital photosynthesis on decomposition
# is not observable) but I cannot reparameterise in terms of alpha,
# otherwise exponential decay may decrease over time because -alpha > tau.

dry_prior %>% 
  prior_posterior_draws(
    posterior_samples = dry_c_samples,
    group = deco_dry_summary %>%
      select(Species, Treatment),
    parameters = c("log_delta_mu", "log_delta_sigma",
                   "log_delta[Species, Treatment]",
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
# Posteriors look smooth

# 2.1.7 Prediction ####
# Global parameters
dry_prior_posterior_global <- dry_prior %>% 
  prior_posterior_draws(
    posterior_samples = dry_c_samples,
    parameters = c("log_delta_mu", "log_delta_sigma",
                   "log_mu_mu", "log_mu_sigma",
                   "log_tau_mu", "log_tau_sigma",
                   "log_epsilon_mu", "log_epsilon_sigma",
                   "log_lambda_mu", "log_lambda_sigma",
                   "log_theta_mu", "log_theta_sigma"),
    format = "short"
  ) %>%
  mutate( # Calculate parameters for unobserved species/treatments
    delta = rnorm( n() , log_delta_mu , log_delta_sigma ) %>% exp(),
    mu = rnorm( n() , log_mu_mu , log_mu_sigma ) %>% exp(),
    tau = rnorm( n() , log_tau_mu , log_tau_sigma ) %>% exp(),
    alpha = delta - tau,
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
dry_prior_posterior <- dry_prior %>% 
  prior_posterior_draws(
    posterior_samples = dry_c_samples,
    group = deco_dry_summary %>%
      select(Species, Treatment),
    parameters = c("log_delta[Species, Treatment]", 
                   "log_mu[Species, Treatment]", 
                   "log_tau[Species, Treatment]", 
                   "log_epsilon[Species, Treatment]", 
                   "log_lambda[Species, Treatment]", 
                   "log_theta[Species, Treatment]"),
    format = "short"
  ) %>% 
  mutate(
    across( # Exponentiate all logged parameters
      starts_with("log"), ~ exp(.x), .names = "{sub('^log_', '', .col)}"
    ),
    alpha = delta - tau
  ) %>% # Remove redundant priors within species
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

# Calculate rounded values
dry_parameters <- dry_prior_posterior %>%
  select(!starts_with(".")) %>%
  filter(Treatment != "Prior") %>%
  mutate(
    # Note I am converting exponential rates to %
    alpha = alpha * 100,
    tau = tau * 100,
    # Calculate parameters on log scale
    log_mu = log(mu),
    log_tau = log(tau),
    log_epsilon = log(epsilon),
    log_lambda = log(lambda),
    log_theta = log(theta)
  ) %>%
  group_by(Species, Treatment) %>%
  summarise(
    across( everything(), list(mean = mean, sd = sd, median = median) ),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(
    alpha = glue("{signif(alpha_mean, 2)} ± {signif(alpha_sd, 2)}"),
    mu_median_rounded = if_else(mu_median < 100, signif(mu_median, 2), signif(mu_median, 3)),
    mu = glue("{mu_median_rounded} ({signif(log_mu_mean, 2)} ± {signif(log_mu_sd, 2)})"),
    tau = glue("{signif(tau_median, 2)} ({signif(log_tau_mean, 2)} ± {signif(log_tau_sd, 2)})"),
    # Only log scale for precision parameters
    epsilon = glue("{signif(log_epsilon_mean, 2)} ± {signif(log_epsilon_sd, 2)}"),
    lambda = glue("{signif(log_lambda_mean, 2)} ± {signif(log_lambda_sd, 2)}"),
    theta = glue("{signif(log_theta_mean, 2)} ± {signif(log_theta_sd, 2)}")
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# Predict across predictor range
dry_prediction <- dry_prior_posterior %>%
  spread_continuous(data = deco_dry_summary %>%
                      # Ensure predictor range starts at 0
                      bind_rows(
                        expand_grid(
                          Day = 0, 
                          Ratio_mean = 1,
                          Species = c("Ecklonia radiata", "Amphibolis griffithii"),
                          Treatment = c("Live", "Dead")
                        ) 
                      ), 
                    # all groups have the same predictor range, so no grouping
                    predictor_name = "Day") %>%
  mutate(
    r_mu = exp(
      Day * alpha - ( alpha + tau ) * mu / 5 * (
        log1p_exp( 5 / mu * ( Day - mu ) ) - log1p_exp( -5 )
      )
    ),
    k = ( alpha + tau ) / ( 1 + exp( 5 / mu * ( Day - mu ) ) ) - tau,
    nu = ( epsilon - theta ) * exp( -lambda * Day ) + theta,
    r = rbetapr( n() , r_mu * ( 1 + nu ) , 2 + nu )
  ) %>% # Summarise predictions
  group_by(Day, Species, Treatment) %>%
  median_qi(r_mu, k, nu, r, .width = c(.5, .8, .9)) %T>%
  print()

# Save progress and clean up
dry_prior_posterior %>%
  write_rds(here("Decomposition", "RDS", "dry_prior_posterior.rds"))
dry_prediction %>%
  write_rds(here("Decomposition", "RDS", "dry_prediction.rds"))

rm(list = ls(pattern = "ratio"), 
   dry_prior, dry_c_model, dry_c_samples,
   dry_nc_model, dry_nc_samples)

# 2.2 Dry mass conventional model ####
# 2.2.1 Prior simulation ####
# The classic exponential decay model for proportions is e^-k*t.
# I am again taking 0.06 d^-1, the mean k value for Ecklonia radiata 
# from Simpkins et al. 2025 (doi: 10.1002/lno.70006) as my prior,
# but with more uncertainty because only one parameter is estimated.
tibble(n = 1:1e3,
       log_k_mu = rnorm( 1e3 , log(0.06) , 0.6 ),
       log_k_sigma = rtnorm( 1e3 , 0 , 0.6 , 0 ),
       log_sigma_mu = rnorm( 1e3 , log(0.1) , 0.3 ),
       log_sigma_sigma = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       k = exp( rnorm( 1e3 , log_k_mu , log_k_sigma ) ),
       sigma = exp( rnorm( 1e3 , log_sigma_mu , log_sigma_sigma ) )) %>%
  expand_grid(Day = deco_dry_summary %$% 
                seq(min(Day), max(Day), length.out = 100)) %>%
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
dry_k_c_model <- here("Decomposition", "Stan", "dry_k_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

dry_k_nc_model <- here("Decomposition", "Stan", "dry_k_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

dry_k_c_samples <- dry_k_c_model$sample(
          data = deco_dry_summary %>%
            select(Day, Ratio_mean,
                   Species, Treatment) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

dry_k_nc_samples <- dry_k_nc_model$sample(
          data = deco_dry_summary %>%
            select(Day, Ratio_mean,
                   Species, Treatment) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

# 2.2.3 Model checks ####
# Rhat
dry_k_c_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000420.

dry_k_nc_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.000145.

# Chains
dry_k_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay()

dry_k_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# No difference in chains

# Pairs
dry_k_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k[1,1]", "log_sigma[1,1]"))
dry_k_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k[1,2]", "log_sigma[1,2]"))

dry_k_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k[2,1]", "log_sigma[2,1]"))
dry_k_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k[2,2]", "log_sigma[2,2]"))

dry_k_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k[1,1]", "log_sigma[1,1]"))
dry_k_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k[1,2]", "log_sigma[1,2]"))

dry_k_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k[2,1]", "log_sigma[2,1]"))
dry_k_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k[2,2]", "log_sigma[2,2]"))
# Correlations are absent for both models

# 2.2.4 Prior-posterior comparison ####
dry_k_prior <- prior_samples(
  model = dry_k_nc_model,
  data = deco_dry_summary %>%
    select(Day, Ratio_mean,
           Species, Treatment) %>%
    compose_data()
  )

dry_k_prior %>% 
  prior_posterior_draws(
    posterior_samples = dry_k_c_samples,
    group = deco_dry_summary %>%
      select(Species, Treatment),
    parameters = c("log_k_mu", "log_k_sigma",
                   "log_k[Species, Treatment]", 
                   "log_sigma_mu", "log_sigma_sigma",
                   "log_sigma[Species, Treatment]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Species", 
    second_group_name = "Treatment"
  )

dry_k_prior %>% 
  prior_posterior_draws(
    posterior_samples = dry_k_nc_samples,
    group = deco_dry_summary %>%
      select(Species, Treatment),
    parameters = c("log_k_mu", "log_k_sigma",
                   "log_k[Species, Treatment]", 
                   "log_sigma_mu", "log_sigma_sigma",
                   "log_sigma[Species, Treatment]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Species", 
    second_group_name = "Treatment"
  )
# Posteriors look similarly smooth
# Proceed with the centred model

# 2.2.5 Prediction ####
# Global parameters
dry_k_prior_posterior_global <- dry_k_prior %>% 
  prior_posterior_draws(
    posterior_samples = dry_k_c_samples,
    parameters = c("log_k_mu", "log_k_sigma",
                   "log_sigma_mu", "log_sigma_sigma"),
    format = "short"
  ) %>%
  mutate( # Calculate parameters for unobserved species/treatments
    k = rnorm( n() , log_k_mu , log_k_sigma ) %>% exp(),
    sigma = rnorm( n() , log_sigma_mu , log_sigma_sigma ) %>% exp()
  ) %>%
  select(starts_with("."), distribution, k, sigma) %T>%
  print()

dry_k_prior_posterior_global %>%
  pivot_longer(cols = -c(starts_with("."), distribution),
               names_to = "parameter") %>%
  group_by(distribution, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n())

# Species/treatment parameters
dry_k_prior_posterior <- dry_k_prior %>% 
  prior_posterior_draws(
    posterior_samples = dry_k_c_samples,
    group = deco_dry_summary %>%
      select(Species, Treatment),
    parameters = c("log_k[Species, Treatment]", 
                   "log_sigma[Species, Treatment]"),
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
  select(starts_with("."), Species, Treatment, k, sigma) %T>%
  print()

dry_k_prior_posterior %>%
  pivot_longer(cols = -c(starts_with("."), Species, Treatment),
               names_to = "parameter") %>%
  group_by(Species, Treatment, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n()) %>%
  print(n = 12)

# Calculate rounded values
dry_k_parameters <- dry_k_prior_posterior %>%
  select(!starts_with(".")) %>%
  filter(Treatment != "Prior") %>%
  mutate(
    t0.5 = log(2) / k, # Calculate half-life (days)
    k = k * 100, # Note I am converting k to % per day
    log_k = log(k), # Calculate on log scale
    log_t0.5 = log(t0.5)
  ) %>%
  group_by(Species, Treatment) %>%
  summarise(
    across( everything(), list(mean = mean, sd = sd, median = median) ),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(
    k = glue("{signif(k_median, 2)} ({signif(log_k_mean, 2)} ± {signif(log_k_sd, 2)})"),
    t0.5 = glue("{signif(t0.5_median, 2)} ({signif(log_t0.5_mean, 2)} ± {signif(log_t0.5_sd, 2)})"),
    sigma = glue("{signif(sigma_mean, 2)} ± {signif(sigma_sd, 2)}")
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# Calculate contrasts
dry_contrast_k <- dry_k_prior_posterior %>%
  filter(Treatment != "Prior") %>%
  select(-sigma) %>%
  pivot_wider(names_from = Treatment, values_from = k) %>%
  mutate(difference = Dead - Live,
         ratio = Dead / Live,
         log_ratio = log10(ratio)) %T>%
  print()

dry_contrast_t0.5 <- dry_k_prior_posterior %>%
  filter(Treatment != "Prior") %>%
  mutate(t0.5 = log(2) / k) %>% # Calculate half-life (days)
  select(-c(k, sigma)) %>%
  pivot_wider(names_from = Treatment, values_from = t0.5) %>%
  mutate(difference = Live - Dead,
         ratio = Live / Dead,
         log_ratio = log10(ratio)) %T>%
  print()

# Summarise contrasts
dry_contrast_k_summary <- dry_contrast_k %>%
  select(!starts_with(".")) %>%
  group_by(Species) %>%
  summarise(
    across(
      everything(), 
      list(
        median = median,
        lower = ~qi(.x, 0.9)[1], 
        upper = ~qi(.x, 0.9)[2],
        mean = mean, sd = sd
      )
    ),
    n = n(),
    P = mean( difference > 0 ) %>% signif(2)
  ) %>%
  ungroup() %>%
  mutate( # Note I am converting k to %
    Live = glue("{signif(Live_median*100, 2)} ({signif(Live_lower*100, 2)}–{signif(Live_upper*100, 2)})"),
    Dead = glue("{signif(Dead_median*100, 2)} ({signif(Dead_lower*100, 2)}–{signif(Dead_upper*100, 2)})"),
    difference = glue("{signif(difference_median*100, 2)} ({signif(difference_lower*100, 2)}–{signif(difference_upper*100, 2)})"),
    ratio = glue("{signif(ratio_median, 2)} ({signif(ratio_lower, 2)}–{signif(ratio_upper, 2)})"),
    log_ratio = glue("{signif(log_ratio_median, 2)} ({signif(log_ratio_lower, 2)}–{signif(log_ratio_upper, 2)})"),
    log_ratio_sym = glue("{signif(log_ratio_mean, 2)} ± {signif(log_ratio_sd, 2)}")
  ) %>%
  select(!(contains("median") | contains("lower") | contains("upper") | 
             contains("mean") | contains("sd"))) %T>%
  print()

dry_contrast_t0.5_summary <- dry_contrast_t0.5 %>%
  select(!starts_with(".")) %>%
  group_by(Species) %>%
  summarise(
    across(
      everything(), 
      list(
        median = median,
        lower = ~qi(.x, 0.9)[1], 
        upper = ~qi(.x, 0.9)[2],
        mean = mean, sd = sd
      )
    ),
    n = n(),
    P = mean( difference > 0 ) %>% signif(2)
  ) %>%
  ungroup() %>%
  mutate(
    Live = glue("{signif(Live_median, 2)} ({signif(Live_lower, 2)}–{signif(Live_upper, 2)})"),
    Dead = glue("{signif(Dead_median, 2)} ({signif(Dead_lower, 2)}–{signif(Dead_upper, 2)})"),
    difference = glue("{signif(difference_median, 2)} ({signif(difference_lower, 2)}–{signif(difference_upper, 2)})"),
    ratio = glue("{signif(ratio_median, 2)} ({signif(ratio_lower, 2)}–{signif(ratio_upper, 2)})"),
    log_ratio = glue("{signif(log_ratio_median, 2)} ({signif(log_ratio_lower, 2)}–{signif(log_ratio_upper, 2)})"),
    log_ratio_sym = glue("{signif(log_ratio_mean, 2)} ± {signif(log_ratio_sd, 2)}")
  ) %>%
  select(!(contains("median") | contains("lower") | contains("upper") | 
             contains("mean") | contains("sd"))) %T>%
  print()

# Predict across predictor range
dry_k_prediction <- dry_k_prior_posterior %>%
  spread_continuous(data = deco_dry_summary %>%
                      # Ensure predictor range starts at 0
                      bind_rows(
                        expand_grid(
                          Day = 0, 
                          Ratio_mean = 1,
                          Species = c("Ecklonia radiata", "Amphibolis griffithii"),
                          Treatment = c("Live", "Dead")
                        ) 
                      ), 
                    # all groups have the same predictor range
                    predictor_name = "Day") %>%
  mutate(
    r_mu = exp( -k * Day ),
    r = rnorm( n() , r_mu , sigma )
  ) %>% # Summarise predictions
  group_by(Day, Species, Treatment) %>%
  median_qi(r_mu, r, .width = c(.5, .8, .9)) %T>%
  print()

# Save progress and clean up
dry_k_prior_posterior %>%
  write_rds(here("Decomposition", "RDS", "dry_k_prior_posterior.rds"))
dry_k_prediction %>%
  write_rds(here("Decomposition", "RDS", "dry_k_prediction.rds"))

rm(dry_k_prior, dry_k_c_model, dry_k_c_samples,
   dry_k_nc_model, dry_k_nc_samples)

# 2.3 Fresh mass macroalgal model ####
# 2.3.1 Visualisation ####
ggplot() +
  geom_point(data = deco_fresh,
             aes(Day, Ratio)) +
  facet_grid(Treatment ~ Species) +
  theme_minimal()

# 2.3.2 Prior simulation ####
# Same as 2.1.2

# 2.3.3 Stan models ####
fresh_c_model <- here("Decomposition", "Stan", "fresh_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

fresh_nc_model <- here("Decomposition", "Stan", "fresh_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

fresh_c_samples <- fresh_c_model$sample(
          data = deco_fresh %>%
            select(Day, Ratio, Species, Treatment) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

fresh_nc_samples <- fresh_nc_model$sample(
          data = deco_fresh %>%
            select(Day, Ratio, Species, Treatment) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

# 2.3.4 Model checks ####
# Rhat
fresh_c_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 78% of rhat above 1.001. rhat = 1.00 ± 0.00170.

fresh_nc_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 8% of rhat above 1.001. rhat = 1.00 ± 0.000303.

# Chains
fresh_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Lost one chain, chains are not good

fresh_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# No chains lost, chains are better

# Pairs
fresh_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[1,1]", "log_mu[1,1]", "log_tau[1,1]"))
fresh_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[1,2]", "log_mu[1,2]", "log_tau[1,2]"))

fresh_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[2,1]", "log_mu[2,1]", "log_tau[2,1]"))
fresh_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[2,2]", "log_mu[2,2]", "log_tau[2,2]"))

fresh_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[3,1]", "log_mu[3,1]", "log_tau[3,1]"))
fresh_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[3,2]", "log_mu[3,2]", "log_tau[3,2]"))

fresh_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[1,1]", "log_mu[1,1]", "log_tau[1,1]"))
fresh_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[1,2]", "log_mu[1,2]", "log_tau[1,2]"))

fresh_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[2,1]", "log_mu[2,1]", "log_tau[2,1]"))
fresh_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[2,2]", "log_mu[2,2]", "log_tau[2,2]"))

fresh_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[3,1]", "log_mu[3,1]", "log_tau[3,1]"))
fresh_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_delta[3,2]", "log_mu[3,2]", "log_tau[3,2]"))
# Correlations are strong and look fairly similar across models.
# Same issue of non-identifiability between delta and tau.

# 2.3.5 Prior-posterior comparison ####
fresh_prior <- prior_samples(
  model = fresh_nc_model,
  data = deco_fresh %>%
    select(Day, Ratio, Species, Treatment) %>%
    compose_data()
  )

fresh_prior %>% 
  prior_posterior_draws(
    posterior_samples = fresh_c_samples,
    group = deco_fresh %>%
      select(Species, Treatment),
    parameters = c("log_delta_mu", "log_delta_sigma",
                   "log_delta[Species, Treatment]", 
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

fresh_prior %>% 
  prior_posterior_draws(
    posterior_samples = fresh_nc_samples,
    group = deco_fresh %>%
      select(Species, Treatment),
    parameters = c("log_delta_mu", "log_delta_sigma",
                   "log_delta[Species, Treatment]", 
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
# Posteriors are smoother for non-centred model.
# The non-centred model is chosen as optimal.

# 2.3.6 Prediction ####
# Global parameters
fresh_prior_posterior_global <- fresh_prior %>% 
  prior_posterior_draws(
    posterior_samples = fresh_nc_samples,
    parameters = c("log_delta_mu", "log_delta_sigma",
                   "log_mu_mu", "log_mu_sigma",
                   "log_tau_mu", "log_tau_sigma",
                   "log_epsilon_mu", "log_epsilon_sigma",
                   "log_lambda_mu", "log_lambda_sigma",
                   "log_theta_mu", "log_theta_sigma"),
    format = "short"
  ) %>%
  mutate( # Calculate parameters for unobserved species/treatments
    delta = rnorm( n() , log_delta_mu , log_delta_sigma ) %>% exp(),
    mu = rnorm( n() , log_mu_mu , log_mu_sigma ) %>% exp(),
    tau = rnorm( n() , log_tau_mu , log_tau_sigma ) %>% exp(),
    alpha = delta - tau,
    epsilon = rnorm( n() , log_epsilon_mu , log_epsilon_sigma ) %>% exp(),
    lambda = rnorm( n() , log_lambda_mu , log_lambda_sigma ) %>% exp(),
    theta = rnorm( n() , log_theta_mu , log_theta_sigma ) %>% exp()
  ) %>%
  select(starts_with("."), distribution, 
         alpha, mu, tau, epsilon, lambda, theta) %T>%
  print()

fresh_prior_posterior_global %>%
  pivot_longer(cols = -c(starts_with("."), distribution),
               names_to = "parameter") %>%
  group_by(distribution, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n())

# Species/treatment parameters
fresh_prior_posterior <- fresh_prior %>% 
  prior_posterior_draws(
    posterior_samples = fresh_nc_samples,
    group = deco_fresh %>%
      select(Species, Treatment),
    parameters = c("log_delta[Species, Treatment]", 
                   "log_mu[Species, Treatment]", 
                   "log_tau[Species, Treatment]", 
                   "log_epsilon[Species, Treatment]", 
                   "log_lambda[Species, Treatment]", 
                   "log_theta[Species, Treatment]"),
    format = "short"
  ) %>% 
  mutate(
    across( # Exponentiate all logged parameters
      starts_with("log"), ~ exp(.x), .names = "{sub('^log_', '', .col)}"
    ),
    alpha = delta - tau
  ) %>% # Remove redundant priors within species
  filter(!(Treatment == "Dead" & distribution == "prior")) %>%
  mutate( # Embed prior in treatment
    Treatment = if_else(
      distribution == "prior", "Prior", Treatment
    ) %>% fct()
  ) %>%
  select(starts_with("."), Species, Treatment, 
         alpha, mu, tau, epsilon, lambda, theta) %T>%
  print()

fresh_prior_posterior %>%
  pivot_longer(cols = -c(starts_with("."), Species, Treatment),
               names_to = "parameter") %>%
  group_by(Species, Treatment, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n()) %>%
  print(n = 54)

# Calculate rounded values for supplementary table
fresh_parameters <- fresh_prior_posterior %>%
  select(!starts_with(".")) %>%
  filter(Treatment != "Prior") %>%
  mutate(
    # Note I am converting exponential rates to %
    alpha = alpha * 100,
    tau = tau * 100,
    # Calculate parameters on log scale
    log_mu = log(mu),
    log_tau = log(tau),
    log_epsilon = log(epsilon),
    log_lambda = log(lambda),
    log_theta = log(theta)
  ) %>%
  group_by(Species, Treatment) %>%
  summarise(
    across( everything(), list(mean = mean, sd = sd, median = median) ),
    n = n()
  ) %>%
  ungroup() %>%
  mutate( # Note I am converting exponential rates to %
    alpha = glue("{signif(alpha_mean, 2)} ± {signif(alpha_sd, 2)}"),
    mu = glue("{signif(mu_median, 2)} ({signif(log_mu_mean, 2)} ± {signif(log_mu_sd, 2)})"),
    tau = glue("{signif(tau_median, 2)} ({signif(log_tau_mean, 2)} ± {signif(log_tau_sd, 2)})"),
    # Only log scale for precision parameters
    epsilon = glue("{signif(log_epsilon_mean, 2)} ± {signif(log_epsilon_sd, 2)}"),
    lambda = glue("{signif(log_lambda_mean, 2)} ± {signif(log_lambda_sd, 2)}"),
    theta = glue("{signif(log_theta_mean, 2)} ± {signif(log_theta_sd, 2)}")
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# Predict across predictor range
fresh_prediction <- fresh_prior_posterior %>%
  spread_continuous(data = deco_fresh %>%
                      # Ensure predictor range starts at 0
                      bind_rows(
                        expand_grid(
                          Day = 0, 
                          Ratio = 1,
                          Species = c("Laminaria digitata", "Ecklonia radiata", "Amphibolis griffithii"),
                          Treatment = c("Live", "Dead")
                        ) 
                      ),
                    group_name = "Species", # different species have different predictor ranges
                    predictor_name = "Day") %>%
  mutate(
    r_mu = exp(
      Day * alpha - ( alpha + tau ) * mu / 5 * (
        log1p_exp( 5 / mu * ( Day - mu ) ) - log1p_exp( -5 )
      )
    ),
    k = ( alpha + tau ) / ( 1 + exp( 5 / mu * ( Day - mu ) ) ) - tau,
    nu = ( epsilon - theta ) * exp( -lambda * Day ) + theta,
    r = rbetapr( n() , r_mu * ( 1 + nu ) , 2 + nu )
  ) %>% # Summarise predictions
  group_by(Day, Species, Treatment) %>%
  median_qi(r_mu, k, nu, r, .width = c(.5, .8, .9)) %T>%
  print()

# Save progress and clean up
fresh_prior_posterior %>%
  write_rds(here("Decomposition", "RDS", "fresh_prior_posterior.rds"))
fresh_prediction %>%
  write_rds(here("Decomposition", "RDS", "fresh_prediction.rds"))

rm(fresh_prior, fresh_c_model, fresh_c_samples,
   fresh_nc_model, fresh_nc_samples)

# 2.4 Fresh mass conventional model ####
# 2.4.1 Prior simulation ####
# Same as 2.2.1

# 2.4.2 Stan model ####
fresh_k_c_model <- here("Decomposition", "Stan", "fresh_k_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

fresh_k_nc_model <- here("Decomposition", "Stan", "fresh_k_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

fresh_k_c_samples <- fresh_k_c_model$sample(
          data = deco_fresh %>%
            select(Day, Ratio,
                   Species, Treatment) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

fresh_k_nc_samples <- fresh_k_nc_model$sample(
          data = deco_fresh %>%
            select(Day, Ratio,
                   Species, Treatment) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

# 2.4.3 Model checks ####
# Rhat
fresh_k_c_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000697.

fresh_k_nc_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.000103.

# Chains
fresh_k_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay()

fresh_k_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# No difference in chains

# Pairs
fresh_k_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k[1,1]", "log_sigma[1,1]"))
fresh_k_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k[1,2]", "log_sigma[1,2]"))

fresh_k_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k[2,1]", "log_sigma[2,1]"))
fresh_k_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k[2,2]", "log_sigma[2,2]"))

fresh_k_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k[3,1]", "log_sigma[3,1]"))
fresh_k_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k[3,2]", "log_sigma[3,2]"))

fresh_k_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k[1,1]", "log_sigma[1,1]"))
fresh_k_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k[1,2]", "log_sigma[1,2]"))

fresh_k_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k[2,1]", "log_sigma[2,1]"))
fresh_k_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k[2,2]", "log_sigma[2,2]"))

fresh_k_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k[3,1]", "log_sigma[3,1]"))
fresh_k_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_k[3,2]", "log_sigma[3,2]"))
# Correlations are absent for both models

# 2.4.4 Prior-posterior comparison ####
fresh_k_prior <- prior_samples(
  model = fresh_k_nc_model,
  data = deco_fresh %>%
    select(Day, Ratio,
           Species, Treatment) %>%
    compose_data()
  )

fresh_k_prior %>% 
  prior_posterior_draws(
    posterior_samples = fresh_k_c_samples,
    group = deco_fresh %>%
      select(Species, Treatment),
    parameters = c("log_k_mu", "log_k_sigma",
                   "log_k[Species, Treatment]", 
                   "log_sigma_mu", "log_sigma_sigma",
                   "log_sigma[Species, Treatment]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Species", 
    second_group_name = "Treatment"
  )

fresh_k_prior %>% 
  prior_posterior_draws(
    posterior_samples = fresh_k_nc_samples,
    group = deco_fresh %>%
      select(Species, Treatment),
    parameters = c("log_k_mu", "log_k_sigma",
                   "log_k[Species, Treatment]", 
                   "log_sigma_mu", "log_sigma_sigma",
                   "log_sigma[Species, Treatment]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Species", 
    second_group_name = "Treatment"
  )
# Posteriors look similarly smooth
# Proceed with the centred model

# 2.4.5 Prediction ####
# Global parameters
fresh_k_prior_posterior_global <- fresh_k_prior %>% 
  prior_posterior_draws(
    posterior_samples = fresh_k_c_samples,
    parameters = c("log_k_mu", "log_k_sigma",
                   "log_sigma_mu", "log_sigma_sigma"),
    format = "short"
  ) %>%
  mutate( # Calculate parameters for unobserved species/treatments
    k = rnorm( n() , log_k_mu , log_k_sigma ) %>% exp(),
    sigma = rnorm( n() , log_sigma_mu , log_sigma_sigma ) %>% exp()
  ) %>%
  select(starts_with("."), distribution, k, sigma) %T>%
  print()

fresh_k_prior_posterior_global %>%
  pivot_longer(cols = -c(starts_with("."), distribution),
               names_to = "parameter") %>%
  group_by(distribution, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n())

# Species/treatment parameters
fresh_k_prior_posterior <- fresh_k_prior %>% 
  prior_posterior_draws(
    posterior_samples = fresh_k_c_samples,
    group = deco_fresh %>%
      select(Species, Treatment),
    parameters = c("log_k[Species, Treatment]", 
                   "log_sigma[Species, Treatment]"),
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
  select(starts_with("."), Species, Treatment, k, sigma) %T>%
  print()

fresh_k_prior_posterior %>%
  pivot_longer(cols = -c(starts_with("."), Species, Treatment),
               names_to = "parameter") %>%
  group_by(Species, Treatment, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n()) %>%
  print(n = 18)

# Calculate rounded values
fresh_k_parameters <- fresh_k_prior_posterior %>%
  select(!starts_with(".")) %>%
  filter(Treatment != "Prior") %>%
  mutate(
    t0.5 = log(2) / k, # Calculate half-life (days)
    k = k * 100, # Note I am converting k to % per day
    log_k = log(k), # Calculate on log scale
    log_t0.5 = log(t0.5)
  ) %>%
  group_by(Species, Treatment) %>%
  summarise(
    across( everything(), list(mean = mean, sd = sd, median = median) ),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(
    k = glue("{signif(k_median, 2)} ({signif(log_k_mean, 2)} ± {signif(log_k_sd, 2)})"),
    t0.5_median_rounded = case_when(
      t0.5_median < 100 ~ signif(t0.5_median, 2),
      t0.5_median < 1e3 ~ signif(t0.5_median, 3),
      t0.5_median < 1e4 ~ signif(t0.5_median, 4),
      TRUE ~ signif(t0.5_median, 5)
    ),
    t0.5 = glue("{t0.5_median_rounded} ({signif(log_t0.5_mean, 2)} ± {signif(log_t0.5_sd, 2)})"),
    sigma = glue("{signif(sigma_mean, 2)} ± {signif(sigma_sd, 2)}")
  ) %>%
  select(!(contains("mean") | contains("sd") | contains("median"))) %T>%
  print()

# Calculate contrasts
fresh_contrast_k <- fresh_k_prior_posterior %>%
  filter(Treatment != "Prior") %>%
  select(-sigma) %>%
  pivot_wider(names_from = Treatment, values_from = k) %>%
  mutate(difference = Dead - Live,
         ratio = Dead / Live,
         log_ratio = log10(ratio)) %T>%
  print()

fresh_contrast_t0.5 <- fresh_k_prior_posterior %>%
  filter(Treatment != "Prior") %>%
  mutate(t0.5 = log(2) / k) %>% # Calculate half-life (days)
  select(-c(k, sigma)) %>%
  pivot_wider(names_from = Treatment, values_from = t0.5) %>%
  mutate(difference = Live - Dead,
         ratio = Live / Dead,
         log_ratio = log10(ratio)) %T>%
  print()

# Summarise contrasts
fresh_contrast_k_summary <- fresh_contrast_k %>%
  select(!starts_with(".")) %>%
  group_by(Species) %>%
  summarise(
    across(
      everything(), 
      list(
        median = median,
        lower = ~qi(.x, 0.9)[1], 
        upper = ~qi(.x, 0.9)[2],
        mean = mean, sd = sd
      )
    ),
    n = n(),
    P = mean( difference > 0 ) %>% signif(2)
  ) %>%
  ungroup() %>%
  mutate( # Note I am converting k to %
    Live = glue("{signif(Live_median*100, 2)} ({signif(Live_lower*100, 2)}–{signif(Live_upper*100, 2)})"),
    Dead = glue("{signif(Dead_median*100, 2)} ({signif(Dead_lower*100, 2)}–{signif(Dead_upper*100, 2)})"),
    difference = glue("{signif(difference_median*100, 2)} ({signif(difference_lower*100, 2)}–{signif(difference_upper*100, 2)})"),
    ratio_median_rounded = case_when(
      ratio_median < 100 ~ signif(ratio_median, 2),
      ratio_median < 1e3 ~ signif(ratio_median, 3),
      TRUE ~ signif(ratio_median, 4)
    ),
    ratio_lower_rounded = case_when(
      ratio_lower < 100 ~ signif(ratio_lower, 2),
      ratio_lower < 1e3 ~ signif(ratio_lower, 3),
      TRUE ~ signif(ratio_lower, 4)
    ),
    ratio_upper_rounded = case_when(
      ratio_upper < 100 ~ signif(ratio_upper, 2),
      ratio_upper < 1e3 ~ signif(ratio_upper, 3),
      TRUE ~ signif(ratio_upper, 4)
    ),
    ratio = glue("{ratio_median_rounded} ({ratio_lower_rounded}–{ratio_upper_rounded})"),
    log_ratio = glue("{signif(log_ratio_median, 2)} ({signif(log_ratio_lower, 2)}–{signif(log_ratio_upper, 2)})"),
    log_ratio_sym = glue("{signif(log_ratio_mean, 2)} ± {signif(log_ratio_sd, 2)}")
  ) %>%
  select(!(contains("median") | contains("lower") | contains("upper") | 
             contains("mean") | contains("sd"))) %T>%
  print()

fresh_contrast_t0.5_summary <- fresh_contrast_t0.5 %>%
  select(!starts_with(".")) %>%
  group_by(Species) %>%
  summarise(
    across(
      everything(), 
      list(
        median = median,
        lower = ~qi(.x, 0.9)[1], 
        upper = ~qi(.x, 0.9)[2],
        mean = mean, sd = sd
      )
    ),
    n = n(),
    P = mean( difference > 0 ) %>% signif(2)
  ) %>%
  ungroup() %>%
  mutate(
    Live_median_rounded = case_when(
      Live_median < 100 ~ signif(Live_median, 2),
      Live_median < 1e3 ~ signif(Live_median, 3),
      Live_median < 1e4 ~ signif(Live_median, 4),
      TRUE ~ signif(Live_median, 5)
    ),
    Live_lower_rounded = case_when(
      Live_lower < 100 ~ signif(Live_lower, 2),
      Live_lower < 1e3 ~ signif(Live_lower, 3),
      Live_lower < 1e4 ~ signif(Live_lower, 4),
      TRUE ~ signif(Live_lower, 5)
    ),
    Live_upper_rounded = case_when(
      Live_upper < 100 ~ signif(Live_upper, 2),
      Live_upper < 1e3 ~ signif(Live_upper, 3),
      Live_upper < 1e4 ~ signif(Live_upper, 4),
      TRUE ~ signif(Live_upper, 5)
    ),
    Live = glue("{Live_median_rounded} ({Live_lower_rounded}–{Live_upper_rounded})"),
    Dead_median_rounded = case_when(
      Dead_median < 100 ~ signif(Dead_median, 2),
      Dead_median < 1e3 ~ signif(Dead_median, 3),
      Dead_median < 1e4 ~ signif(Dead_median, 4),
      TRUE ~ signif(Dead_median, 5)
    ),
    Dead_lower_rounded = case_when(
      Dead_lower < 100 ~ signif(Dead_lower, 2),
      Dead_lower < 1e3 ~ signif(Dead_lower, 3),
      Dead_lower < 1e4 ~ signif(Dead_lower, 4),
      TRUE ~ signif(Dead_lower, 5)
    ),
    Dead_upper_rounded = case_when(
      Dead_upper < 100 ~ signif(Dead_upper, 2),
      Dead_upper < 1e3 ~ signif(Dead_upper, 3),
      Dead_upper < 1e4 ~ signif(Dead_upper, 4),
      TRUE ~ signif(Dead_upper, 5)
    ),
    Dead = glue("{Dead_median_rounded} ({Dead_lower_rounded}–{Dead_upper_rounded})"),
    difference_median_rounded = case_when(
      difference_median < 100 ~ signif(difference_median, 2),
      difference_median < 1e3 ~ signif(difference_median, 3),
      difference_median < 1e4 ~ signif(difference_median, 4),
      TRUE ~ signif(difference_median, 5)
    ),
    difference_lower_rounded = case_when(
      difference_lower < 100 ~ signif(difference_lower, 2),
      difference_lower < 1e3 ~ signif(difference_lower, 3),
      difference_lower < 1e4 ~ signif(difference_lower, 4),
      TRUE ~ signif(difference_lower, 5)
    ),
    difference_upper_rounded = case_when(
      difference_upper < 100 ~ signif(difference_upper, 2),
      difference_upper < 1e3 ~ signif(difference_upper, 3),
      difference_upper < 1e4 ~ signif(difference_upper, 4),
      TRUE ~ signif(difference_upper, 5)
    ),
    difference = glue("{difference_median_rounded} ({difference_lower_rounded}–{difference_upper_rounded})"),
    ratio_median_rounded = case_when(
      ratio_median < 100 ~ signif(ratio_median, 2),
      ratio_median < 1e3 ~ signif(ratio_median, 3),
      ratio_median < 1e4 ~ signif(ratio_median, 4),
      TRUE ~ signif(ratio_median, 5)
    ),
    ratio_lower_rounded = case_when(
      ratio_lower < 100 ~ signif(ratio_lower, 2),
      ratio_lower < 1e3 ~ signif(ratio_lower, 3),
      ratio_lower < 1e4 ~ signif(ratio_lower, 4),
      TRUE ~ signif(ratio_lower, 5)
    ),
    ratio_upper_rounded = case_when(
      ratio_upper < 100 ~ signif(ratio_upper, 2),
      ratio_upper < 1e3 ~ signif(ratio_upper, 3),
      ratio_upper < 1e4 ~ signif(ratio_upper, 4),
      TRUE ~ signif(ratio_upper, 5)
    ),
    ratio = glue("{ratio_median_rounded} ({ratio_lower_rounded}–{ratio_upper_rounded})"),
    log_ratio = glue("{signif(log_ratio_median, 2)} ({signif(log_ratio_lower, 2)}–{signif(log_ratio_upper, 2)})"),
    log_ratio_sym = glue("{signif(log_ratio_mean, 2)} ± {signif(log_ratio_sd, 2)}")
  ) %>%
  select(!(contains("median") | contains("lower") | contains("upper") | 
             contains("mean") | contains("sd"))) %T>%
  print()

# Predict across predictor range
fresh_k_prediction <- fresh_k_prior_posterior %>%
  spread_continuous(data = deco_fresh %>%
                      # Ensure predictor range starts at 0
                      bind_rows(
                        expand_grid(
                          Day = 0, 
                          Ratio = 1,
                          Species = c("Ecklonia radiata", "Amphibolis griffithii"),
                          Treatment = c("Live", "Dead")
                        ) 
                      ), 
                    # all groups have the same predictor range
                    predictor_name = "Day") %>%
  mutate(
    r_mu = exp( -k * Day ),
    r = rnorm( n() , r_mu , sigma )
  ) %>% # Summarise predictions
  group_by(Day, Species, Treatment) %>%
  median_qi(r_mu, r, .width = c(.5, .8, .9)) %T>%
  print()

# Save progress and clean up
fresh_k_prior_posterior %>%
  write_rds(here("Decomposition", "RDS", "fresh_k_prior_posterior.rds"))
fresh_k_prediction %>%
  write_rds(here("Decomposition", "RDS", "fresh_k_prediction.rds"))

rm(fresh_k_prior, fresh_k_c_model, fresh_k_c_samples,
   fresh_k_nc_model, fresh_k_nc_samples)

# 3. Figures ####
# 3.1 Combine predictions ####
prediction <- dry_prediction %>%
  mutate(Mass = "Dry" %>% fct()) %>%
  bind_rows(fresh_prediction %>%
              mutate(Mass = "Fresh" %>% fct())) %>%
  mutate(Species = Species %>% 
           fct_relevel("Laminaria digitata", "Ecklonia radiata")) %T>%
  print()

contrast_k <- dry_contrast_k %>%
  mutate(Mass = "Dry" %>% fct()) %>%
  bind_rows(fresh_contrast_k %>%
              mutate(Mass = "Fresh" %>% fct())) %>%
  mutate(Species = Species %>% 
           fct_relevel("Laminaria digitata", "Ecklonia radiata")) %T>%
  print()

contrast_t0.5 <- dry_contrast_t0.5 %>%
  mutate(Mass = "Dry" %>% fct()) %>%
  bind_rows(fresh_contrast_t0.5 %>%
              mutate(Mass = "Fresh" %>% fct())) %>%
  mutate(Species = Species %>% 
           fct_relevel("Laminaria digitata", "Ecklonia radiata")) %T>%
  print()

posterior <- dry_k_prior_posterior %>%
  mutate(Mass = "Dry" %>% fct()) %>%
  bind_rows(fresh_k_prior_posterior %>%
              mutate(Mass = "Fresh" %>% fct())) %>%
  mutate(
    Species = Species %>% 
      fct_relevel("Laminaria digitata", "Ecklonia radiata"),
    Treatment = Treatment %>% fct_relevel("Dead"),
    t0.5 = log(2) / k # Calculate half-life (days)
  ) %T>%
  print()

# 3.2 Combine data ####
deco <- deco_dry_summary %>%
  mutate(Mass = "Dry" %>% fct()) %>%
  bind_rows(deco_fresh %>%
              mutate(Mass = "Fresh" %>% fct()) %>%
              rename(Ratio_mean = Ratio)) %>%
  mutate(Species = Species %>% 
           fct_relevel("Laminaria digitata", "Ecklonia radiata")) %T>%
  print()

# 3.3 Figure 2 ####
# 3.3.1 Figure 2a ####
Fig_2a <- prediction %>%
  filter(
    Treatment != "Prior" & 
      !(Species %in% c("Ecklonia radiata", "Amphibolis griffithii") & Mass == "Fresh")
  ) %>%
  ggplot() +
    geom_hline(yintercept = c(0, 1)) +
    geom_pointrange(data = deco %>%
                      filter(!(Species %in% c("Ecklonia radiata", "Amphibolis griffithii") &
                                 Mass == "Fresh")),
                    aes(Day, Ratio_mean, ymin = Ratio_lwr, ymax = Ratio_upr, colour = Species),
                    shape = 16, alpha = 0.5, size = 0.5) +
    geom_ribbon(data = . %>% filter(.width == 0.9),
                aes(Day, ymin = r.lower, ymax = r.upper,
                    fill = Species), alpha = 0.3) +
    geom_ribbon(aes(Day, ymin = r_mu.lower, ymax = r_mu.upper,
                    alpha = factor(.width), fill = Species)) +
    geom_line(aes(Day, r_mu, colour = Species)) +
    scale_fill_manual(values = c("#333b08", "#c3b300", "#4a7518"), guide = "none") +
    scale_colour_manual(values = c("#333b08", "#c3b300", "#4a7518"), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_y_continuous(breaks = seq(0, 1.8, 0.6),
                       labels = scales::label_number(accuracy = c(1, rep(0.1, 3)))) +
    labs(x = "Detrital age (days)",
         y = expression("Relative detrital mass ("*italic(m)*"/"*italic(m)[0]*")")) +
    coord_cartesian(xlim = c(0, 80), ylim = c(0, 1.8),
                    expand = F, clip = "off") +
    facet_grid2(Treatment ~ Species,
                switch = "y",
                strip = strip_nested(text_y = element_text(angle = 0, hjust = 0, vjust = 1))) +
    mytheme +
    theme(strip.text.x = element_text(face = "italic", hjust = 0),
          plot.margin = margin(0, 0.5, 0, 0.2, unit = "cm"))

Fig_2a
# Safely ignore warning, which is due to Laminaria digitata not having a standard 
# deviation for geom_pointrange.

# 3.3.2 Figure 2b ####
require(ggridges)
Fig_2b <- posterior %>%
  filter(
    Treatment != "Prior" & 
      !(Species %in% c("Ecklonia radiata", "Amphibolis griffithii") & Mass == "Fresh")
  ) %>%
  mutate( # Laminaria digitata half-life is best visualised on the log scale
    t0.5 = if_else(Species == "Laminaria digitata", log10(t0.5), t0.5)
  ) %>%
  ggplot() +
    stat_density_ridges(aes(t0.5, y = Treatment, fill = Species), 
                        colour = NA, n = 2^10, from = rep(0, 3), to = c(6, 80, 120),
                        scale = 2, bandwidth = c(6*0.02, 80*0.02, 120*0.02)) +
    geom_text(
      data = tibble(
        Species = posterior %$% levels(Species) %>% rep(each = 2) %>% fct(),
        Treatment = c("Live", "Dead") %>% rep(3) %>% fct(),
        label = c("Live", "Dead", NA %>% rep(4))
      ),
      aes(x = -1.565, y = Treatment, label = label),
      family = "Futura", size = 12, size.unit = "pt",
      hjust = 0, vjust = 0
    ) +
    scale_fill_manual(values = c("#333b08", "#c3b300", "#4a7518"), guide = "none") +
    facet_grid(~ Species, scales = "free") +
    facetted_pos_scales(
      x = list(
        Species == "Laminaria digitata" ~
          scale_x_continuous(limits = c(0, 6),
                             breaks = seq(0, 6, by = 2),
                             labels = scales::label_math(10^.x),
                             oob = scales::oob_keep),
        Species == "Ecklonia radiata" ~
          scale_x_continuous(limits = c(0, 80),
                             breaks = seq(0, 80, by = 20),
                             oob = scales::oob_keep),
        Species == "Amphibolis griffithii" ~
          scale_x_continuous(limits = c(0, 120),
                             breaks = seq(0, 120, by = 30),
                             oob = scales::oob_keep)
        )
      ) +
    labs(x = "Half-life (days)") +
    coord_cartesian(expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          strip.text = element_blank())

Fig_2b
# Safely ignore warning, which is due to intentional NAs in geom_text.

# 3.3.3 Figure 2c ####
Fig_2c <- contrast_t0.5 %>%
  filter(!(Species %in% c("Ecklonia radiata", "Amphibolis griffithii") & 
             Mass == "Fresh")) %>%
  ggplot() +
    stat_density_ridges(aes(log_ratio, y = 0, fill = Species),
                        colour = NA, n = 2^10, from = c(1, 0.4, 0), to = c(4, 1.2, 0.4),
                        bandwidth = c(3*0.02, 0.8*0.02, 0.4*0.02)) +
    scale_fill_manual(values = c("#333b08", "#c3b300", "#4a7518"), guide = "none") +
    facet_wrap(~ Species, scales = "free") + # I am wrapping here so that densities are the same height
    facetted_pos_scales(
      x = list(
        Species == "Laminaria digitata" ~
          scale_x_continuous(limits = c(1, 4),
                             breaks = seq(1, 4, 1),
                             labels = scales::label_math(10^.x),
                             oob = scales::oob_keep),
        Species == "Ecklonia radiata" ~
          scale_x_continuous(limits = c(0.4, 1.2),
                             breaks = seq(0.4, 1.2, 0.2),
                             labels = scales::label_math(10^.x),
                             oob = scales::oob_keep),
        Species == "Amphibolis griffithii" ~
          scale_x_continuous(limits = c(0, 0.4),
                             breaks = seq(0, 0.4, 0.1),
                             labels = scales::label_math(10^.x),
                             oob = scales::oob_keep)
        )
      ) +
    labs(x = "Relative half-life (Live / Dead)") +
    coord_cartesian(expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          strip.text = element_blank())

Fig_2c

# 3.3.4 Combine panels ####
require(patchwork)
Fig_2 <- ( Fig_2a / Fig_2b / Fig_2c ) +
  plot_layout(heights = c(1, 0.3, 0.06))

Fig_2 %>%
  ggsave(filename = "Fig_2.pdf", path = "Figures",
         device = cairo_pdf, height = 15, width = 20, units = "cm")

# 3.4 Figure S3 ####
# 3.4.1 Figure S3a ####
Fig_S3a <- prediction %>%
  filter(Treatment != "Prior" & Mass == "Fresh") %>%
  ggplot() +
    geom_hline(yintercept = c(0, 1)) +
    geom_point(data = deco %>% filter(Mass == "Fresh"),
               aes(Day, Ratio_mean, colour = Species),
               shape = 16, alpha = 0.5, size = 2.4) +
    geom_ribbon(data = . %>% filter(.width == 0.9),
                aes(Day, ymin = r.lower, ymax = r.upper,
                    fill = Species), alpha = 0.3) +
    geom_ribbon(aes(Day, ymin = r_mu.lower, ymax = r_mu.upper,
                    alpha = factor(.width), fill = Species)) +
    geom_line(aes(Day, r_mu, colour = Species)) +
    scale_fill_manual(values = c("#333b08", "#c3b300", "#4a7518"), guide = "none") +
    scale_colour_manual(values = c("#333b08", "#c3b300", "#4a7518"), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_y_continuous(breaks = seq(0, 1.8, 0.6),
                       labels = scales::label_number(accuracy = c(1, rep(0.1, 3)))) +
    labs(x = "Detrital age (days)",
         y = expression("Relative detrital mass ("*italic(m)*"/"*italic(m)[0]*")")) +
    coord_cartesian(xlim = c(0, 80), ylim = c(0, 1.8),
                    expand = F, clip = "off") +
    facet_grid2(Treatment ~ Species,
                switch = "y",
                strip = strip_nested(text_y = element_text(angle = 0, hjust = 0, vjust = 1))) +
    mytheme +
    theme(strip.text.x = element_text(face = "italic", hjust = 0),
          plot.margin = margin(0, 0.5, 0, 0.2, unit = "cm"))

Fig_S3a

# 3.4.2 Figure S3b ####
Fig_S3b <- posterior %>%
  filter(Treatment != "Prior" & Mass == "Fresh") %>%
  mutate( # Laminaria digitata and Amphibolis griffithii half-lifes are 
          # best visualised on the log scale
    t0.5 = if_else(Species == "Ecklonia radiata", t0.5, log10(t0.5))
  ) %>%
  ggplot() +
    stat_density_ridges(aes(t0.5, y = Treatment, fill = Species, alpha = Species), 
                        colour = NA, n = 2^10, from = c(0, 0, 2), to = c(6, 100, 5),
                        scale = 2, bandwidth = c(6*0.02, 100*0.02, 3*0.02)) +
    geom_text(
      data = tibble(
        Species = posterior %$% levels(Species) %>% rep(each = 2) %>% fct(),
        Treatment = c("Live", "Dead") %>% rep(3) %>% fct(),
        label = c("Live", "Dead", NA %>% rep(4))
      ),
      aes(x = -1.565, y = Treatment, label = label),
      family = "Futura", size = 12, size.unit = "pt",
      hjust = 0, vjust = 0
    ) +
    scale_fill_manual(values = c("#333b08", "#c3b300", "#4a7518"), guide = "none") +
    scale_alpha_manual(values = c(1, 1, 0.7), guide = "none") +
    facet_grid(~ Species, scales = "free") +
    facetted_pos_scales(
      x = list(
        Species == "Laminaria digitata" ~
          scale_x_continuous(limits = c(0, 6),
                             breaks = seq(0, 6, by = 2),
                             labels = scales::label_math(10^.x),
                             oob = scales::oob_keep),
        Species == "Ecklonia radiata" ~
          scale_x_continuous(limits = c(0, 100),
                             breaks = seq(0, 100, by = 25),
                             oob = scales::oob_keep),
        Species == "Amphibolis griffithii" ~
          scale_x_continuous(limits = c(2, 5),
                             labels = scales::label_math(10^.x),
                             oob = scales::oob_keep)
        )
      ) +
    labs(x = "Half-life (days)") +
    coord_cartesian(expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          strip.text = element_blank())

Fig_S3b
# Safely ignore warning, which is due to intentional NAs in geom_text.

# 3.4.3 Figure S3c ####
Fig_S3c <- contrast_t0.5 %>%
  filter(Mass == "Fresh") %>%
  ggplot() +
    stat_density_ridges(aes(log_ratio, y = 0, fill = Species),
                        colour = NA, n = 2^10, from = c(1, 0.3, -1), to = c(4, 1.1, 2),
                        bandwidth = c(3*0.02, 0.8*0.02, 3*0.02)) +
    geom_vline(
      data = tibble(
        Species = contrast_t0.5 %$% levels(Species) %>% fct(),
        x = c(NA, NA, 0)
      ),
      aes(xintercept = x)
    ) +
    scale_fill_manual(values = c("#333b08", "#c3b300", "#4a7518"), guide = "none") +
    facet_wrap(~ Species, scales = "free") + # I am wrapping here so that densities are the same height
    facetted_pos_scales(
      x = list(
        Species == "Laminaria digitata" ~ 
          scale_x_continuous(limits = c(1, 4), 
                             breaks = seq(1, 4, by = 1),
                             labels = scales::label_math(10^.x),
                             oob = scales::oob_keep),
        Species == "Ecklonia radiata" ~ 
          scale_x_continuous(limits = c(0.3, 1.1), 
                             breaks = seq(0.3, 1.1, by = 0.2),
                             labels = scales::label_math(10^.x),
                             oob = scales::oob_keep),
        Species == "Amphibolis griffithii" ~ 
          scale_x_continuous(limits = c(-1, 2),
                             breaks = seq(-1, 2, by = 1),
                             labels = scales::label_math(10^.x), # scale_x_log10 + lable_log doesn't work
                             oob = scales::oob_keep) # so need to replace hyphen with minus afterwards
        )
      ) +
    labs(x = "Relative half-life (Live / Dead)") +
    coord_cartesian(expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          strip.text = element_blank())

Fig_S3c

# 3.4.4 Combine panels ####
Fig_S3 <- ( Fig_S3a / Fig_S3b / Fig_S3c ) +
  plot_layout(heights = c(1, 0.3, 0.06))

Fig_S3 %>%
  ggsave(filename = "Fig_S3.pdf", path = "Figures",
         device = cairo_pdf, height = 15, width = 20, units = "cm")

# 3.5 Figure S4 ####
Fig_S4 <- posterior %>%
  filter(Treatment != "Prior") %>% 
  mutate(Mass = Mass %>% fct_relevel("Fresh"),
         Transparent = Species == "Amphibolis griffithii" & Mass == "Fresh") %>%
  ggplot() +
    stat_density_ridges(aes(k, y = Treatment, fill = Species, alpha = Transparent), 
                        colour = NA, n = 2^10, from = rep(0, 3), to = c(0.12, 0.18, 0.02), 
                        scale = 2, bandwidth = c(0.12*0.02, 0.18*0.02, 0.02*0.02)) +
    geom_text(
      data = . %>% droplevels() %$%
        expand_grid(
          Species = levels(Species) %>% fct(),
          Mass = levels(Mass) %>% fct(),
          Treatment = levels(Treatment) %>% fct()
        ) %>%
        mutate(
          label = if_else(Species == "Laminaria digitata",
                          Treatment, NA)
        ),
      aes(x = -0.0421, y = Treatment, label = label),
      family = "Futura", size = 12, size.unit = "pt",
      hjust = 0, vjust = 0
    ) +
    scale_fill_manual(values = c("#333b08", "#c3b300", "#4a7518"), guide = "none") +
    scale_alpha_manual(values = c(1, 0.7), guide = "none") +
    facet_grid2(Mass ~ Species, scales = "free", switch = "y",
                strip = strip_nested(text_y = element_text(angle = 0, hjust = 0, vjust = 1))) +
    facetted_pos_scales(
      x = list(
        Species == "Laminaria digitata" ~
          scale_x_continuous(limits = c(0, 0.12),
                             breaks = seq(0, 0.12, by = 0.04),
                             labels = scales::label_number(accuracy = c(1, 0.01 %>% rep(3))),
                             oob = scales::oob_keep),
        Species == "Ecklonia radiata" ~
          scale_x_continuous(limits = c(0, 0.18),
                             breaks = seq(0, 0.18, by = 0.06),
                             labels = scales::label_number(accuracy = c(1, 0.01 %>% rep(3))),
                             oob = scales::oob_keep),
        Species == "Amphibolis griffithii" ~
          scale_x_continuous(limits = c(0, 0.02),
                             breaks = seq(0, 0.02, by = 0.01),
                             labels = scales::label_number(accuracy = c(1, 0.01 %>% rep(2))),
                             oob = scales::oob_keep)
        )
      ) +
    labs(x = expression("Exponential decay ("*italic(k)*", day"^-1*")")) +
    coord_cartesian(expand = F, clip = "off") +
    mytheme +
    theme(plot.margin = margin(0.2, 0.5, 0.2, 0, unit = "cm"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          strip.text.x = element_text(face = "italic"),
          strip.text.y = element_text(face = "bold"))

Fig_S4
# First warning is benign, just informs that a panel in the grid is missing 
# (data for Laminaria digitata dry mass don't exist).
# Safely ignore second warning, which is due to intentional NAs in geom_text.

Fig_S4 %>%
  ggsave(filename = "Fig_S4.pdf", path = "Figures",
         device = cairo_pdf, height = 8, width = 16, units = "cm")

# 4. Tables ####
# 4.1 Combine estimates ####
k_parameters <- dry_k_parameters %>%
  mutate(Mass = "Dry") %>%
  bind_rows(fresh_k_parameters %>%
              mutate(Mass = "Fresh")) %>%
  mutate(Species = Species %>% 
           fct_relevel("Laminaria digitata", "Ecklonia radiata"),
         Mass = Mass %>% fct_relevel("Fresh")) %T>%
  print()
  
parameters <- dry_parameters %>%
  mutate(Mass = "Dry") %>%
  bind_rows(fresh_parameters %>%
              mutate(Mass = "Fresh")) %>%
  mutate(Species = Species %>% 
           fct_relevel("Laminaria digitata", "Ecklonia radiata"),
         Mass = Mass %>% fct_relevel("Fresh")) %T>%
  print()

contrast_k_summary <- dry_contrast_k_summary %>%
  mutate(Mass = "Dry" %>% fct()) %>%
  bind_rows(fresh_contrast_k_summary %>%
              mutate(Mass = "Fresh" %>% fct())) %>%
  mutate(Species = Species %>% 
           fct_relevel("Laminaria digitata", "Ecklonia radiata")) %T>%
  print()
  
contrast_t0.5_summary <- dry_contrast_t0.5_summary %>%
  mutate(Mass = "Dry" %>% fct()) %>%
  bind_rows(fresh_contrast_t0.5_summary %>%
              mutate(Mass = "Fresh" %>% fct())) %>%
  mutate(Species = Species %>% 
           fct_relevel("Laminaria digitata", "Ecklonia radiata")) %T>%
  print()

# 4.2 Table 1 ####
# I favour reporting half-life contrasts because they are more intuitive
# (besides, the log_ratio of the t0.5 contrast is identical to the k contrast).
Table_1.1 <- contrast_t0.5_summary %>%
  filter(!(Species %in% c("Ecklonia radiata", "Amphibolis griffithii") & 
             Mass == "Fresh")) %>%
  select(Species, Live, Dead, difference, log_ratio, log_ratio_sym, P) %>%
  arrange(Species) %T>%
  print()

Table_1.1 %>%
  write_csv(here("Tables", "Table_1.1.csv"))

require(officer)
read_docx() %>%
  body_add_table(value = Table_1.1) %>%
  print(target = here("Tables", "Table_1.1.docx"))

# Save alternative anyway
Table_1.1_k <- contrast_k_summary %>%
  filter(!(Species %in% c("Ecklonia radiata", "Amphibolis griffithii") & 
             Mass == "Fresh")) %>%
  select(Species, Live, Dead, difference, log_ratio, log_ratio_sym, P) %>%
  arrange(Species) %T>%
  print()

Table_1.1_k %>%
  write_csv(here("Tables", "Table_1.1_k.csv"))

read_docx() %>%
  body_add_table(value = Table_1.1_k) %>%
  print(target = here("Tables", "Table_1.1_k.docx"))

# 4.3 Table S3 ####
Table_S3 <- k_parameters %>% 
  select(Species, Mass, Treatment, k, t0.5) %>%
  left_join(
    parameters %>% 
      select(Species, Mass, Treatment, alpha, mu, tau),
    by = c("Species", "Mass", "Treatment")
  ) %>%
  arrange(Species, Mass) %T>%
  print()

Table_S3 %>%
  write_csv(here("Tables", "Table_S3.csv"))

read_docx() %>%
  body_add_table(value = Table_S3) %>%
  print(target = here("Tables", "Table_S3.docx"))

# 4.4 Table S4 ####
Table_S4 <- contrast_t0.5_summary %>%
  filter(Mass == "Fresh") %>%
  select(Species, Live, Dead, difference, log_ratio, log_ratio_sym, P) %>%
  arrange(Species) %T>%
  print()

Table_S4 %>%
  write_csv(here("Tables", "Table_S4.csv"))

read_docx() %>%
  body_add_table(value = Table_S4) %>%
  print(target = here("Tables", "Table_S4.docx"))

Table_S4_k <- contrast_k_summary %>%
  filter(Mass == "Fresh") %>%
  select(Species, Live, Dead, difference, log_ratio, log_ratio_sym, P) %>%
  arrange(Species) %T>%
  print()

Table_S4_k %>%
  write_csv(here("Tables", "Table_S4_k.csv"))

read_docx() %>%
  body_add_table(value = Table_S4_k) %>%
  print(target = here("Tables", "Table_S4_k.docx"))