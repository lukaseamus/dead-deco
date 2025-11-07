#### Dead deco: kelp decomposition without physiology ####
#### Luka Seamus Wright                               ####

# 1. Prepare data ####
# 1.1 Load data ####
require(tidyverse)
require(magrittr)

decouk <- read_csv("Decomposition_UK.csv", col_types = list("f", "f", "f")) %>%
  mutate(Days = Hours / 24,
         # Replace 0 with small constant within balance measurement error
         # because the model is undefined for y = 0
         Final = if_else(Final == 0, 1e-5, Final)) %T>%
  print()

decoau <- read_csv("Decomposition_AU.csv", col_types = list("f", "f", "f", "f")) %>%
  mutate(Deployment = Deployment %>% dmy(),
         Retrieval = Retrieval %>% dmy(),
         Days = Deployment %--% Retrieval / ddays(),
         Final = if_else(Final == 0, 1e-5, Final),
         Dry = if_else(Dry == 0, 1e-6, Dry)) %T>%
  print()

ratio <- read_csv("Ratio.csv", col_types = list("f", "f")) %>%
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
require(here)
ratio_model <- here("Stan", "ratio.stan") %>% 
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
          adapt_delta = 0.99 # force sampler to slow down
        )

# 1.2.4 Model checks ####
# Rhat
ratio_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000876.

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
# Posteriors are not constrained by priors.

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
# Species               beta_mean beta_sd cv_mean cv_sd     n
# Unobserved                 0.28  0.11     0.36  0.71  80000
# Prior                      0.29  0.14     0.97  1.4   80000
# Ecklonia radiata           0.27  0.0063   0.13  0.018 80000
# Amphibolis griffithii      0.25  0.004    0.087 0.012 80000

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

# 1.2.7 Visualisation ####
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
                 panel.spacing = unit(1, "cm"),
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
               size = 2, alpha = 0.5, shape = 16) +
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
    mytheme

Fig_S1 %>%
  ggsave(filename = "Fig_S1_unedited.pdf", path = "Figures",
         device = cairo_pdf, height = 10, width = 21, units = "cm")
# For some reason annotation font turns bold, so needs to be edited.

# 1.3 Remaining proportion ####
# 1.3.1 Dry proportion ####
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
    Proportion = Dry / Initial_dry
  ) %T>%
  print()

# Summarise
deco_dry_summary <- deco_dry %>%
  group_by(ID, Species, Treatment, Tank, Deployment, Initial, 
           Retrieval, Final, Dry, Initials, Notes, Days) %>%
  summarise(
    beta_mean = mean(beta),
    beta_sd = sd(beta),
    Initial_dry_mean = mean(Initial_dry),
    Initial_dry_sd = sd(Initial_dry),
    Proportion_mean = mean(Proportion),
    Proportion_sd = sd(Proportion),
    Proportion_lwr = qi(Proportion, .width = 0.99)[1],
    Proportion_upr = qi(Proportion, .width = 0.99)[2],
    n = n()
  ) %>%
  ungroup() %T>%
  print()

# 1.3.2 Fresh proportion ####
# Combine and calculate for both datasets 
deco_wet <- decoau %>%
  select(ID, Species, Treatment, Initial, 
         Final, Days, Initials) %>%
  bind_rows(decouk %>% select(-Hours)) %>%
  mutate(Proportion = Final / Initial) %T>%
  print()

# 2. Model data ####
# 2.1 Dry proportion ####
# 2.1.1 Visualisation ####
ggplot() +
  geom_violin(data = deco_dry,
              aes(Days, Proportion, group = ID)) +
  facet_grid(Treatment ~ Species) +
  theme_minimal()
# This won't do for visualisation.

ggplot() +
  geom_pointrange(data = deco_dry_summary,
                  aes(Days, Proportion_mean,
                      ymin = Proportion_mean - Proportion_sd,
                      ymax = Proportion_mean + Proportion_sd),
                  position = position_dodge2(width = 1)) +
  facet_grid(Treatment ~ Species) +
  theme_minimal()
# But ± s.d. clearly doesn't represent the skew.

ggplot() +
  geom_pointrange(data = deco_dry_summary,
                  aes(Days, Proportion_mean,
                      ymin = Proportion_lwr,
                      ymax = Proportion_upr),
                  position = position_dodge2(width = 1)) +
  facet_grid(Treatment ~ Species) +
  theme_minimal()
# Better. It is clear that modelling the measurment error
# will dramatically increase the weight of values near zero.
# This seems a somewhat artificial inflation of uncertainty
# since a scenario where only wet mass is measured as in the
# case of Laminaria digitata yields similar data with more
# apparent certainty for large values. Conversely, if we had 
# been able to measure initial dry mass, the trajectory would
# likely look similar, again without the skew in weighting.

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
# include the smallest positive values.
tibble(n = 1:1e3,
       alpha_mu = rnorm( 1e3 , -0.02 , 0.01 ), 
       log_mu_mu = rnorm( 1e3 , log(40) , 0.3 ),
       log_tau_mu = rnorm( 1e3 , log(0.06) , 0.4 ),
       alpha_sigma = rtnorm( 1e3 , 0 , 0.01 , 0 ), 
       log_mu_sigma = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       log_tau_sigma = rtnorm( 1e3 , 0 , 0.4 , 0 ),
       log_epsilon_mu = rnorm( 1e3 , log(4e4) , 0.3 ),
       log_lambda_mu = rnorm( 1e3 , log(0.3) , 0.3 ),
       log_theta_mu = rnorm( 1e3 , log(500) , 0.3 ),
       log_epsilon_sigma = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       log_lambda_sigma = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       log_theta_sigma = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       alpha = rnorm( 1e3 , alpha_mu , alpha_sigma ),
       mu = exp( rnorm( 1e3 , log_mu_mu , log_mu_sigma ) ),
       tau = exp( rnorm( 1e3 , log_tau_mu , log_tau_sigma ) ),
       epsilon = exp( rnorm( 1e3 , log_epsilon_mu , log_epsilon_sigma ) ),
       lambda = exp( rnorm( 1e3 , log_lambda_mu , log_lambda_sigma ) ),
       theta = exp( rnorm( 1e3 , log_theta_mu , log_theta_sigma ) )) %>%
  expand_grid(Days = deco_dry_summary %$% 
                seq(min(Days), max(Days), length.out = 100)) %>%
  mutate(
    p_mu = exp(
      Days * alpha - ( alpha + tau ) * mu / 5 * (
        log1p_exp( 5 / mu * ( Days - mu ) ) - log1p_exp( -5 )
      )
    ),
    k = ( alpha + tau ) / ( 1 + exp( 5 / mu * ( Days - mu ) ) ) - tau,
    nu = theta + exp( log(epsilon - theta) - lambda * Days ),
    p = rbetapr( n() , p_mu * ( 1 + nu ) , 2 + nu )
  ) %>%
  pivot_longer(cols = c(p_mu, k, nu, p),
               names_to = "parameter") %>%
  ggplot(aes(Days, value, group = n)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    facet_wrap(~parameter, scale = "free", nrow = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 2.1.3 Stan model ####
dry_c_model <- here("Stan", "dry_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

dry_nc_model <- here("Stan", "dry_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

dry_c_samples <- dry_c_model$sample(
          data = deco_dry_summary %>%
            select(Days, Proportion_mean, Proportion_sd,
                   Species, Treatment) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4,
          adapt_delta = 0.99,
          treedepth = 15
        ) %T>%
  print(max_rows = 200)

dry_nc_samples <- dry_nc_model$sample(
          data = deco_dry_summary %>%
            select(Days, Proportion_mean, Proportion_sd,
                   Species, Treatment) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4,
          adapt_delta = 0.99,
          treedepth = 15
        ) %T>%
  print(max_rows = 200)

# 2.1.4 Model checks ####
# Rhat
dry_c_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 2% of rhat above 1.001. rhat = 1.00 ± 0.000217.

dry_nc_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.000121.

# Chains
dry_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay()

dry_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Chains are far better for non-centred model

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


# 2.1.5 Prior-posterior comparison ####
dry_c_prior <- prior_samples(
  model = dry_c_model,
  data = deco_dry_summary %>%
    select(Days, Proportion_mean, Proportion_sd,
           Species, Treatment) %>%
    compose_data(),
  adapt_delta = 0.99,
  max_treedepth = 15
  )
# cannot effectively sample centred hierarchical priors

dry_nc_prior <- prior_samples(
  model = dry_nc_model,
  data = deco_dry_summary %>%
    select(Days, Proportion_mean, Proportion_sd,
           Species, Treatment) %>%
    compose_data()
)

# use only non-centred priors
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

# 2.1.6 Prediction ####
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

# Predict across predictor range
dry_prediction <- dry_prior_posterior %>%
  spread_continuous(data = deco_dry_summary %>%
                      # Ensure predictor range starts at 0
                      bind_rows(
                        expand_grid(
                          Days = 0, 
                          Proportion_mean = 1,
                          Species = c("Ecklonia radiata", "Amphibolis griffithii"),
                          Treatment = c("Live", "Dead")
                        ) 
                      ), 
                    # all groups have the same predictor range
                    predictor_name = "Days", length = 200) %>%
  mutate(
    p_mu = exp(
      Days * alpha - ( alpha + tau ) * mu / 5 * (
        log1p_exp( 5 / mu * ( Days - mu ) ) - log1p_exp( -5 )
      )
    ),
    k = ( alpha + tau ) / ( 1 + exp( 5 / mu * ( Days - mu ) ) ) - tau,
    nu = ( epsilon - theta ) * exp( -lambda * Days ) + theta,
    p = rbetapr( n() , p_mu * ( 1 + nu ) , 2 + nu )
  ) %>% # Summarise predictions
  group_by(Days, Species, Treatment) %>%
  median_qi(p_mu, k, nu, p, .width = c(.5, .8, .9)) %T>%
  print()

Fig_1a <- dry_prediction %>%
  filter(Treatment != "Prior") %>%
  ggplot() +
    geom_hline(yintercept = c(0, 1)) +
    geom_line(aes(Days, p, colour = Species)) +
    geom_ribbon(aes(Days, ymin = p.lower, ymax = p.upper,
                    alpha = factor(.width), fill = Species)) +
    geom_pointrange(data = deco_dry_summary,
                    aes(Days, Proportion_mean,
                        ymin = Proportion_lwr,
                        ymax = Proportion_upr,
                        colour = Species),
                    shape = 16, alpha = 0.5) +
    # geom_point(data = deco_dry_summary,
    #            aes(Days, Proportion_mean, colour = Species),
    #            shape = 16, alpha = 0.5) +
    scale_fill_manual(values = c("#c3b300", "#4a7518"), guide = "none") +
    scale_colour_manual(values = c("#c3b300", "#4a7518"), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_y_continuous(breaks = seq(0, 1.8, 0.6),
                       labels = scales::label_number(accuracy = c(1, rep(0.1, 3)))) +
    labs(x = "Detrital age (days)",
         y = "Remaining proportion of mass") +
    coord_cartesian(xlim = c(0, 80), ylim = c(0, 1.8),
                    expand = F, clip = "off") +
    facet_grid2(Treatment ~ Species,
                switch = "y",
                strip = strip_nested(text_y = element_text(angle = 0, hjust = 0, vjust = 1))) +
    mytheme +
    theme(strip.text.x = element_text(face = "italic", hjust = 0),
          plot.margin = margin(0, 0.5, 0, 0.2, unit = "cm"))



Fig_1a %>%
  ggsave(filename = "Fig_1a.pdf", path = "Figures",
         device = cairo_pdf, height = 10, width = 20, units = "cm")



# 2.2 Wet proportion ####
# 2.2.1 Visualisation ####

ggplot() +
  geom_point(data = deco_wet,
             aes(Days, Proportion)) +
  facet_grid(Treatment ~ Species) +
  theme_minimal()

# 3. Visualisation ####


# 4. Naive models ####


