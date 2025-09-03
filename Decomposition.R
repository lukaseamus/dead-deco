# 1. Prepare data ####
# 1.1 Load data ####
require(tidyverse)
require(magrittr)

deco <- read_csv("Decomposition.csv", col_types = list("f", "f", "f", "f")) %>%
  mutate(Deployment = Deployment %>% dmy(),
         Retrieval = Retrieval %>% dmy(),
         Days = Deployment %--% Retrieval / ddays()) %T>%
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
# No rhat above 1.001. rhat = 1.00 ± 0.0000531.

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

# Predict across predictor range
require(truncnorm) # R doesn't have a native truncated normal
ratio_prediction <- ratio_prior_posterior %>%
  spread_continuous(data = ratio, 
                    predictor_name = "Fresh",
                    group_name = "Species") %>%
  mutate(mu = beta * Fresh,
         sigma = cv * mu,
         Dry = rtruncnorm( n() , mean = mu , sd = sigma , a = 0 )) %T>%
  print()

# Summarise predictions
ratio_prediction_summary <- ratio_prediction %>%
  group_by(Species, Fresh) %>%
  mean_qi(mu, Dry, .width = c(.5, .8, .9)) %T>%
  print()

# 1.2.7 Visualisation ####
# Summarise parameters for annotation
require(glue)
ratio_annotation <- ratio_prior_posterior %>%
  group_by(Species) %>%
  summarise(beta_mean = mean(beta),
            beta_sd = sd(beta),
            cv_mean = mean(cv),
            cv_sd = sd(cv),
            n = n()) %>%
  mutate(label = glue(
        "μ = {signif(beta_mean, 2)} ± {signif(beta_sd, 2)} × x\n
         σ = {signif(cv_mean, 2)} ± {signif(cv_sd, 2)} × μ"
      )
    ) %T>%
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
                 strip.text = element_text(size = 12, hjust = 0, face = "italic"),
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
    # geom_ribbon(data = ratio_prediction_summary %>%
    #               filter(Species == "Prior" &
    #                        .width == 0.9) %>%
    #               select(-Species),
    #             aes(Fresh, ymin = Dry.lower, ymax = Dry.upper),
    #             alpha = 0.1) +
    geom_line(data = ratio_prediction_summary %>%
                filter(!Species %in% c("Unobserved", "Prior")),
              aes(Fresh, mu, colour = Species)) +
    # geom_ribbon(data = ratio_prediction_summary %>%
    #               filter(!Species %in% c("Unobserved", "Prior")),
    #             aes(Fresh, ymin = mu.lower, ymax = mu.upper,
    #                 fill = Species, alpha = factor(.width))) +
    geom_ribbon(data = ratio_prediction_summary %>%
                  filter(!Species %in% c("Unobserved", "Prior")),
                aes(Fresh, ymin = Dry.lower, ymax = Dry.upper,
                    fill = Species, alpha = factor(.width))) +
    geom_text(data = ratio_annotation %>%
                filter(!Species %in% c("Unobserved", "Prior")),
              aes(x = c(33, 4.95), y = c(2.5, 0.335), label = label),
              family = "Futura", size = 3.5, hjust = 0,
              lineheight = 0.8) +
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
         device = cairo_pdf, height = 10, width = 21, 
         units = "cm")
# For some reason annotation font turns bold, so needs to be edited.

# 1.3 Remaining proportion ####
# 1.3.1 Dry proportion ####
deco_dry <- deco %>%
  left_join(
    ratio_prior_posterior %>%
      select(Species, beta, cv),
    by = "Species",
    relationship = "many-to-many"
  ) %>%
  mutate(
    mu = beta * Initial,
    sigma = cv * mu,
    Initial_dry = rtruncnorm( n() , mean = mu , sd = sigma , a = 0 ),
    Proportion = Dry / Initial_dry
  ) %T>%
  print()

deco_dry_summary <- deco_dry %>%
  group_by(ID, Species, Treatment) %>%
  summarise(
    Initial = unique(Initial), 
    Dry = unique(Dry), 
    Days = unique(Days), 
    beta_mean = mean(beta),
    beta_sd = sd(beta),
    mu_mean = mean(mu),
    mu_sd = sd(mu),
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

ggplot() +
  geom_pointrange(data = deco_dry_summary,
                  aes(Days, Proportion_mean,
                      ymin = Proportion_lwr,
                      ymax = Proportion_upr),
                  position = position_dodge2(width = 1)) +
  facet_grid(Treatment ~ Species) +
  theme_minimal()

# 1.3.2 Fresh proportion ####
deco %<>%
  mutate(Proportion = Final / Initial) %T>%
  print()

ggplot() +
  geom_point(data = deco,
             aes(Days, Proportion)) +
  facet_grid(Treatment ~ Species) +
  theme_minimal()

# 2. Model data ####
# 2.1 Dry proportion ####
# 2.1.1 Prior simulation ####
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




