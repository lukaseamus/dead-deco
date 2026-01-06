#### Dead deco: kelp decomposition without physiology ####
#### Luka Seamus Wright                               ####

# 1. Prepare data ####
# 1.1 Turnover relationship ####
# 1.1.1 Load data ####
require(tidyverse)
require(magrittr)
require(here)

carbon <- here("Turnover", "Carbon.csv") %>% 
  read_csv(col_types = list("f", "f", "f", "f", "f")) %T>%
  print()

# 1.1.2 Select relevant variables ####
PB <- carbon %>%
  mutate(
    PB_plant = `Net primary production (g C m^-2 y^-1)` / `Plant biomass (g C m^-2)`,
    # In steady state all these measures of detrital P/B are the same:
    PB_detritus = `Detrital production (g C m^-2 y^-1)` / `Detrital biomass (g C m^-2)`,
    PB_detritus_k = `k (d^-1)` * 365,
    PB_detritus_deco = `Decomposition (g C m^-2 y^-1)` / `Detrital biomass (g C m^-2)`
  ) %>%
  select(Reference, Community, PB_plant, PB_detritus, PB_detritus_k, PB_detritus_deco) %T>%
  print(n = 941)
# In many cases detrital P/B, k and absolute decomposition were clearly derived
# from the same observations and are not independent. See Cebrián & Lartigue 2004
# (doi: 10.1890/03-4019) for details, especially methods on pp. 240ff.

# 1.1.3 Clean up ####
# I favour detrital P/B derived from detrital production since this was used in
# the original study (Cebrián & Duarte 1995, doi: 10.1126/science.268.5217.1606),
# but also want to include other values if they can replace NAs. I favour P/B
# derived from k over P/B derived from absolute decomposition because k is
# already in the relative form required for P/B and is what I am ultimately predicting.
PB %<>%
  mutate(PB_detritus = coalesce(PB_detritus, PB_detritus_k, PB_detritus_deco)) %>%
  select(Reference, Community, PB_plant, PB_detritus) %T>%
  print(n = 941)

# Filter for pairs
PB %<>% 
  drop_na(PB_plant, PB_detritus) %>%
  droplevels() %T>%
  print()

# Check that pairs are unique per Reference and Community
PB %>%
  group_by(Reference, Community, PB_plant, PB_detritus) %>%
  filter(n() > 1)
# No duplication

# 1.1.4 Check grouping ####
PB %>% nrow()
PB %$% nlevels(Community)
PB %$% nlevels(Reference)
PB %>% count(Community)

# 1.1.5 Visualise ####
PB %>%
  ggplot(aes(PB_plant, PB_detritus)) +
    geom_point() +
    theme_minimal()
# Hence the log-log plot in the original study (Figure 3b
# in doi: 10.1126/science.268.5217.1606)

PB %>%
  ggplot(aes(log10(PB_plant), log10(PB_detritus))) +
  geom_point() +
  theme_minimal()
# Clearer relationship

PB %>%
  ggplot(aes(log10(PB_plant), log10(PB_detritus),
             colour = Community)) +
  geom_point() +
  theme_minimal()
# Relationship breaks down for some groups

# 1.2 Kelp turnover ####
# 1.2.1 Load data ####
PB_kelp <- here("Turnover", "Turnover.csv") %>% 
  read_csv(col_types = list("f", "c", "c", "f", "f")) %T>%
  print()

# 1.2.2 Check grouping ####
PB_kelp %>% nrow() # 223 observations
PB_kelp %$% nlevels(Species) # 33 species
PB_kelp %$% nlevels(Reference) # 77 references
PB_kelp %$% nlevels(Level) # 3 measurement levels

PB_kelp %>%
  distinct(Reference, Species) %>%
  count(Species) %>%
  print(n = 33)
# Most of the 33 species are just in one reference

PB_kelp %>% 
  distinct(Reference, Species) %>%
  count(Reference) %>%
  print(n = 77)
# Most of the 77 references just report one species

PB_kelp %>% 
  group_by(Level) %>%
  summarise(n_Species = n_distinct(Species),
            n_Reference = n_distinct(Reference)) %>%
  ungroup()
# Plenty of species and references in each measurement level

# 1.2.3 Visualise ####
PB_kelp %>%
  ggplot(aes(Turnover)) +
    geom_density() +
    geom_jitter(aes(y = 0), 
                height = 0.01, 
                shape = 16, 
                alpha = 0.2) +
    theme_minimal()

# 1.3 Forest turnover ####
# 1.3.1 Load ForC data ####
# These data can be downloaded from github.com/forc-db/ForC/tree/master/data.
# The data publciation is Anderson-Teixeira 2018 (doi: 10.1002/ecy.2229).
PB_forest_ForC <- here("Turnover", "ForC_measurements.csv") %>% 
  read_csv() %T>%
  print()
# Some variables are problematic
problems(PB_forest_ForC) %>% print(n = 30)
PB_forest_ForC[,c(37, 40)] # Not important

# 1.3.2 Select relevant variables #### 
# Consult ForC_variables.csv from github.com/forc-db/ForC/tree/master/data
PB_forest_ForC %<>%
  select(citation.ID, measurement.ID, sites.sitename, plot.name, 
         stand.age, dominant.veg, dominant.life.form, scientific.name, 
         variable.name, date, mean, required.citations) %>%
  filter(
    variable.name %in% 
      c("NPP_1_C", "NPP_0_C", "biomass_C")
  ) %>%
  droplevels() %>%
  mutate(
    Reference = citation.ID %>% str_remove("_[^_]*$") %>% fct_relevel(sort),
    Site = sites.sitename %>% fct(),
    Plot = plot.name %>% fct(),
    Review = required.citations %>% str_remove("_[^_]*$") %>% fct()
  ) %>%
  select(-c(citation.ID, sites.sitename, plot.name, required.citations)) %T>%
  print()

# 1.3.3 Tidy data #### 
PB_forest_ForC %<>%
  select(-measurement.ID) %>% # Measurment ID is unique per variable
  pivot_wider(names_from = variable.name,
              values_from = mean,
              # Summarise duplicates as mean per site
              values_fn = mean) %T>%
  print()

# 1.3.4 Coalesce NPP and calculate P/B #### 
PB_forest_ForC %<>%
  mutate(
    NPP = coalesce(NPP_1_C, NPP_0_C),
    Turnover = NPP / biomass_C
  ) %>%
  select(-c(NPP_1_C, NPP_0_C)) %>%
  drop_na(Turnover) %>%
  filter(Turnover > 0) %>%
  droplevels() %T>%
  print()

PB_forest_ForC %>% nrow() # 205 observations
PB_forest_ForC %$% nlevels(Reference) # 36 references
PB_forest_ForC %$% levels(Review) 
# Anderson-Teixeira 2018 (doi: 10.1002/ecy.2229) is the reference to cite

# 1.3.5 Load Just's data #### 
# These data are contained in the object "carbon".
PB_forest_Just <- carbon %>%
  filter(Biome == "Woodland") %>%
  mutate(
    Turnover = `Net primary production (g C m^-2 y^-1)` / `Plant biomass (g C m^-2)`
  ) %>%
  select(Reference, Turnover) %>%
  distinct(Reference, Turnover) %>% # Filter out duplicates within reference
  drop_na(Turnover) %>% # Filter out NAs
  droplevels() %T>%
  print()

PB_forest_Just %>% nrow() # 67 observations
PB_forest_Just %$% nlevels(Reference) # 31 references

# 1.3.6 Check redundancy and merge ####
PB_forest_ForC %$% levels(Reference)
PB_forest_Just %$% levels(Reference)
# Don't see any duplicate references

PB_forest_ForC %>%
  distinct(Reference) %>%
  mutate(Reference = Reference %>% str_extract("^[^_ ]+")) %>%
  inner_join(
    PB_forest_Just %>%
      distinct(Reference) %>%
      mutate(Reference = Reference %>% str_extract("^[^_ ]+"))
  )
# Indeed no redundancy

PB_forest <- PB_forest_ForC %>%
  select(Reference, Turnover) %>%
  bind_rows(PB_forest_Just) %T>%
  print()

PB_forest %>% nrow() # 272 observations
PB_forest %$% nlevels(Reference) # 67 references

# 1.3.7 Visualise ####
PB_forest %>%
  ggplot(aes(Turnover)) +
    geom_density() +
    geom_jitter(aes(y = 0), 
                height = 0.01, 
                shape = 16, 
                alpha = 0.2) +
    theme_minimal()

# 1.4 Kelp k ####
# 1.4.1 Load data ####
k <- here("Turnover", "k.csv") %>% 
  read_csv(col_types = list("f", "c", "c", "f", "c", "f", "f")) %T>%
  print()

# 1.4.2 Select relevant data ####
# I am only interested in data on natural subtidal kelp decomposition.
# This excludes wrack and dissolved organic matter, and experiments that
# prekilled, buried, or otherwise killed (anoxia, darkness) detritus.
# Lab experiments are included but only when they are not smallscale
# but mesocosm scale.
k %<>%
  filter(Kelp, !Wrack, !Dissolved, !Prekilled, !Buried, 
         !Anoxic, !Darkness, !Smallscale) %T>%
  print()

k %>% nrow() # 358 observations after filtering
k %$% range(k)

# 1.4.3 Check grouping ####
# This is meta-meta-data with repeated references within reviews. I am
# counting Filbee-Dexter et al. 2024 and Moten Foldager Pedersen unpublished
# as one merged review since they were not conducted independently.
k %<>%
  mutate(
    Review = if_else(
      Review %>% str_detect("Filbee-Dexter et al. 2024"),
      "Filbee-Dexter et al. 2024", Review
    ) %>% fct()
  ) %T>%
  print()

k %>% 
  distinct(Reference, Review) %>%
  count(Review)
# The number of references per review is highly skewed. 
# Luka Seamus Wright unpublished has most references and 
# two reviews only have one reference. Let's look closer:
k %>% 
  filter(Review %in% c("Enríquez et al. 1993", "Cebrián & Lartigue 2004"))
# The reference is the same for both reviews and both report
# the same values. These reviews are not worthwhile including.
k %<>% 
  filter(!Review %in% c("Enríquez et al. 1993", "Cebrián & Lartigue 2004")) %>%
  droplevels() %T>%
  print()

k %>% 
  distinct(Reference, Review) %>%
  count(Reference) %>%
  print(n = 30)
# There are 22 references across the 4 remaining reviews and most
# are contained in more than one review.

k %>% 
  count(Review, Reference, Species) %>%
  print(n = 100)
# There are usually several observations per Review-Reference-Species
# combination. All species names are complete.
k %$% nlevels(Species) # There are 14 species.

# 1.4.4 Visualise ####
k %>%
  ggplot(aes(k)) +
  geom_density() +
  geom_jitter(aes(y = 0), 
              height = 0.01, 
              shape = 16, 
              alpha = 0.2) +
  theme_minimal()
# Some values are negative which is not allowed if we assume decay.
# I could replace these values with a positive constant below minimum 
# positive k value
k %>% 
  filter(k > 0) %$%
  min(k)
# 1.40143e-05 is the minimum positive value so 1e-5 seems reasonable.
# However, this would collapse variation in the negative values and
# since the underlying frequentist models that k was derived from
# are unconstrained, a normal distribution would likely be the most
# representative likelihood. The alternatives would be shifting the
# data before analysis and reversing the shift in the posteriors or 
# simply removing negative k values. The first is arbitrary and the
# second removes data that can guide inference.

# 2. Turnover relationship ####
# 2.1 Prior simulation ####
# I am using the same log10-log10 model as in the original study,
# but with partial pooling across plant ecosystems.
require(extraDistr)
tibble(n = 1:1e3,
       alpha_mu = rnorm( 1e3 , 0 , 1 ), # one order of magnitude standard deviations
       alpha_sigma = rtnorm( 1e3 , 0 , 1 , 0 ), # half-normal distribution
       beta_mu = rnorm( 1e3 , 0 , 1 ),
       beta_sigma = rtnorm( 1e3 , 0 , 1 , 0 ),
       alpha = rnorm( 1e3 , alpha_mu , alpha_sigma ),
       beta = rnorm( 1e3 , beta_mu , beta_sigma ),
       sigma = rexp( 1e3 , 1 )) %>%
  expand_grid(log_PB_plant = PB %$% 
                seq(min(log10(PB_plant)), max(log10(PB_plant)), length.out = 100)) %>%
  mutate(
    mu = alpha + beta * log_PB_plant,
    log_PB_detritus = rnorm( n() , mu , sigma )
  ) %>%
  pivot_longer(cols = c(mu, log_PB_detritus),
               names_to = "parameter") %>%
  ggplot(aes(log_PB_plant, value, group = n)) +
    geom_hline(yintercept = PB %$% range(log10(PB_detritus))) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    facet_wrap(~parameter, scale = "free", nrow = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Covers all possible values.

# 2.2 Stan models ####
require(cmdstanr)
PB_c_model <- here("Turnover", "Stan", "PB_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

PB_nc_model <- here("Turnover", "Stan", "PB_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
PB_c_samples <- PB_c_model$sample(
          data = PB %>%
            select(PB_plant, PB_detritus, Community) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

PB_nc_samples <- PB_nc_model$sample(
          data = PB %>%
            select(PB_plant, PB_detritus, Community) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

# 2.3 Model checks ####
# Rhat
PB_c_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 21% of rhat above 1.001. rhat = 1.00 ± 0.000822.

PB_nc_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.000111.

# Chains
require(bayesplot)
PB_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay()

PB_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Non-centred chains are better

# Pairs
PB_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_mu", "beta_mu"))
PB_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[1]", "beta[1]"))
PB_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[2]", "beta[2]"))
# Some correlation between alpha and beta but
# not bad and not evident for global parameters.

PB_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_mu", "beta_mu"))
PB_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[1]", "beta[1]"))
PB_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha[2]", "beta[2]"))
# Correlations are similar

# 2.4 Prior-posterior comparison ####
source("functions.R")
PB_prior <- prior_samples(
  model = PB_nc_model,
  data = PB %>%
    select(PB_plant, PB_detritus, Community) %>%
    compose_data()
  )

PB_prior %>% 
  prior_posterior_draws(
    posterior_samples = PB_c_samples,
    group = PB %>% select(Community),
    parameters = c("alpha_mu", "alpha_sigma",
                   "alpha[Community]", 
                   "beta_mu", "beta_sigma",
                   "beta[Community]",
                   "sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Community"
  )

PB_prior %>% 
  prior_posterior_draws(
    posterior_samples = PB_nc_samples,
    group = PB %>% select(Community),
    parameters = c("alpha_mu", "alpha_sigma",
                   "alpha[Community]", 
                   "beta_mu", "beta_sigma",
                   "beta[Community]",
                   "sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Community"
  )
# Posteriors look somewhat more constrained
# for the non-centred model. Proceed with
# non-centred model.

# 2.5 Prediction ####
# 2.5.1 Global parameters ####
PB_prior_posterior_global <- PB_prior %>% 
  prior_posterior_draws(
    posterior_samples = PB_nc_samples,
    parameters = c("alpha_mu", "alpha_sigma",
                   "beta_mu", "beta_sigma",
                   "sigma"),
    format = "short"
  ) %>%
  mutate( # Calculate parameters for unobserved plant ecosystems
    alpha = rnorm( n() , alpha_mu , alpha_sigma ),
    beta = rnorm( n() , beta_mu , beta_sigma )
  ) %>%
  select(starts_with("."), distribution, 
         alpha_mu, beta_mu, alpha, beta, sigma) %T>%
  print()

PB_prior_posterior_global %>%
  pivot_longer(cols = -c(starts_with("."), distribution),
               names_to = "parameter") %>%
  group_by(distribution, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n())

# 2.5.2 Community parameters ####
PB_prior_posterior_community <- PB_prior %>% 
  prior_posterior_draws(
    posterior_samples = PB_nc_samples,
    group = PB %>% select(Community),
    parameters = c("alpha[Community]", 
                   "beta[Community]",
                   "sigma"),
    format = "short"
  ) %>% # Remove redundant priors within community
  filter(distribution == "prior" & Community == "marine benthic microalgal beds" |
           distribution == "posterior") %>%
  mutate( # Embed prior in community
    Community = if_else(
      distribution == "prior", "Prior", Community
    ) %>% fct()
  ) %>%
  select(starts_with("."), Community, 
         alpha, beta, sigma) %T>%
  print()

PB_prior_posterior_community %>%
  pivot_longer(cols = -c(starts_with("."), Community),
               names_to = "parameter") %>%
  group_by(Community, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n()) %>%
  print(n = 30)

# 2.5.3 Combine parameters ####
PB_prior_posterior <- PB_prior_posterior_global %>%
  filter(distribution == "posterior") %>% # priors are redundant
  select(!c(distribution, alpha_mu, beta_mu)) %>%
  mutate(Community = "Unobserved" %>% fct()) %>%
  bind_rows(PB_prior_posterior_community) %T>%
  print()
  
# 2.5.4 Summarise ####
require(glue)
PB_parameters <- PB_prior_posterior %>%
  filter(Community != "Prior") %>%
  select(!starts_with(".")) %>%
  group_by(Community) %>%
  summarise(
    across( everything(), list(mean = mean, sd = sd) ),
    P = mean( beta > 0 ),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(
    alpha = glue("{signif(alpha_mean, 2)} ± {signif(alpha_sd, 2)}"),
    beta = glue("{signif(beta_mean, 2)} ± {signif(beta_sd, 2)}"),
    sigma = glue("{signif(sigma_mean, 2)} ± {signif(sigma_sd, 2)}"),
    P = signif(P, 2)
  ) %>%
  select(!(contains("mean") | contains("sd"))) %T>%
  print()

# Export parameter table
PB_parameters %>%
  write_csv(here("Tables", "PB_parameters.csv"))

require(officer)
read_docx() %>%
  body_add_table(value = PB_parameters) %>%
  print(target = here("Tables", "PB_parameters.docx"))

# 2.5.5 Predict across predictor range ####
PB_prediction <- PB_prior_posterior %>%
  spread_continuous(data = PB %>% # Add transformed predictor
                      mutate(log_PB_plant = log10(PB_plant)), 
                    group_name = "Community",
                    predictor_name = "log_PB_plant") %>%
  mutate(
    mu = alpha + beta * log_PB_plant,
    log_PB_detritus = rnorm( n() , mu , sigma )
  ) %>% # Summarise predictions
  group_by(log_PB_plant, Community) %>%
  median_qi(mu, log_PB_detritus, .width = c(.5, .8, .9)) %T>%
  print()

# Save progress and clean up
PB_prior_posterior %>%
  write_rds(here("Turnover", "RDS", "PB_prior_posterior.rds"))
PB_prediction %>%
  write_rds(here("Turnover", "RDS", "PB_prediction.rds"))

rm(PB_prior, PB_c_model, PB_c_samples,
   PB_nc_model, PB_nc_samples)

# 3. Kelp turnover ####
# 3.1 Prior simulation ####
# I could use a gamma generalised linear model because turnover
# is strictly positive. Its exact distribution is unknown. Beta
# prime would best describe a ratio of gamma distributions which
# is how turnover is often calculated (P/B) but the turnover rate
# itself isn't necessarily beta prime. Lognormal is another option.
# I generally prefer gamma, which would make the turnover time
# inverse gamma distributed. But since I want to visualise these
# data on the log scale, the lognormal makes sense. I will partially 
# pool across references, species and measurement levels. 

# I will take the mean marine macroalgal turnover from Cebrián and
# Duarte 1994 (doi: 10.2307/2390077), Cebrián 1999 (doi: 10.1086/303244)
# etc. as my prior. These data are all contained in object "carbon".
carbon %>%
  filter(Community == "marine macroalgal beds") %>%
  mutate(PB = `Net primary production (g C m^-2 y^-1)` / `Plant biomass (g C m^-2)`) %>%
  drop_na(PB) %>%
  summarise(PB_min = min(PB),
            PB_median = median(PB),
            PB_max = max(PB),
            PB_mean = mean(PB),
            PB_sd = sd(PB),
            n = n())
# Some very extreme and unlikely values in there. Clearly the distribution
# is right-skewed. I will go with the median and large uncertainty.

tibble(n = 1:1e4,
       alpha_mu = rnorm( 1e4 , log(7.54) , 1 ),
       alpha_sigma_s = rtnorm( 1e4 , 0 , 0.3 , 0 ), # half-normal distribution
       alpha_sigma_l = rtnorm( 1e4 , 0 , 0.3 , 0 ),
       alpha_sigma_r = rtnorm( 1e4 , 0 , 0.3 , 0 ),
       alpha_s = rnorm( 1e4 , alpha_mu , alpha_sigma_s ),
       alpha_l = rnorm( 1e4 , 0 , alpha_sigma_l ),
       alpha_r = rnorm( 1e4 , 0 , alpha_sigma_r ),
       mu = alpha_s + alpha_l + alpha_r,
       exp_mu = exp(mu),
       sigma = rexp( 1e4 , 1 ),
       Turnover = rlnorm( 1e4 , mu , sigma )) %>%
  pivot_longer(cols = c(exp_mu, Turnover),
               names_to = "parameter") %>%
  ggplot(aes(value)) +
    geom_vline(xintercept = PB_kelp %$% range(Turnover)) +
    geom_density(alpha = 0.05) +
    scale_x_continuous(limits = c(0, 30),
                       oob = scales::oob_keep) +
    coord_cartesian(expand = F, clip = "off") +
    facet_wrap(~parameter, scale = "free", nrow = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Covers all possible values.

# 3.2 Stan models ####
PB_kelp_c_model <- here("Turnover", "Stan", "PB_kelp_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

PB_kelp_nc_model <- here("Turnover", "Stan", "PB_kelp_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

PB_kelp_c_samples <- PB_kelp_c_model$sample(
          data = PB_kelp %>%
            select(Turnover, Species, Level, Reference) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

PB_kelp_nc_samples <- PB_kelp_nc_model$sample(
          data = PB_kelp %>%
            select(Turnover, Species, Level, Reference) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

# 3.3 Model checks ####
# Rhat
PB_kelp_c_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 3% of rhat above 1.001. rhat = 1.00 ± 0.000270.

PB_kelp_nc_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000730.

# Chains
PB_kelp_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay()

PB_kelp_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Chains are better for non-centred model

# Pairs
PB_kelp_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_mu", "alpha_s[1]", "alpha_l[1]", "alpha_r[1]"))
PB_kelp_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_mu", "alpha_s[2]", "alpha_l[2]", "alpha_r[2]"))
# Some non-identifiability in alpha.

PB_kelp_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_mu", "alpha_s[1]", "alpha_l[1]", "alpha_r[1]"))
PB_kelp_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_mu", "alpha_s[2]", "alpha_l[2]", "alpha_r[2]"))
# Correlations are similar

# 3.4 Prior-posterior comparison ####
PB_kelp_prior <- prior_samples(
  model = PB_kelp_nc_model,
  data = PB_kelp %>%
    select(Turnover, Species, Level, Reference) %>%
    compose_data()
  )

PB_kelp_prior %>% 
  prior_posterior_draws(
    posterior_samples = PB_kelp_c_samples,
    group = PB_kelp %>% select(Species, Level),
    parameters = c("alpha_mu", "alpha_sigma_s",
                   "alpha_s[Species]", 
                   "alpha_sigma_l", "alpha_l[Level]",
                   "sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Species",
    second_group_name = "Level"
  )

PB_kelp_prior %>% 
  prior_posterior_draws(
    posterior_samples = PB_kelp_c_samples,
    group = PB_kelp %>% select(Reference),
    parameters = c("alpha_mu", "alpha_sigma_r",
                   "alpha_r[Reference]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Reference"
  )

PB_kelp_prior %>% 
  prior_posterior_draws(
    posterior_samples = PB_kelp_nc_samples,
    group = PB_kelp %>% select(Species, Level),
    parameters = c("alpha_mu", "alpha_sigma_s",
                   "alpha_s[Species]", 
                   "alpha_sigma_l", "alpha_l[Level]",
                   "sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Species",
    second_group_name = "Level"
  )

PB_kelp_prior %>% 
  prior_posterior_draws(
    posterior_samples = PB_kelp_nc_samples,
    group = PB_kelp %>% select(Reference),
    parameters = c("alpha_mu", "alpha_sigma_r",
                   "alpha_r[Reference]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Reference"
  )
# Posteriors look similar. Proceed with
# non-centred model because of marginally
# better rhat.

# 3.5 Prediction ####
# 3.5.1 Global parameters ####
PB_kelp_prior_posterior_global <- PB_kelp_prior %>% 
  prior_posterior_draws(
    posterior_samples = PB_kelp_nc_samples,
    parameters = c("alpha_mu", "alpha_sigma_s",
                   "alpha_sigma_l", "alpha_sigma_r",
                   "sigma"),
    format = "short"
  ) %>%
  mutate(
    # Calculate mu for unobserved species, levels and references
    mu = rnorm( n() , alpha_mu , alpha_sigma_s ) + 
          rnorm( n() , 0 , alpha_sigma_l ) +
          rnorm( n() , 0 , alpha_sigma_r ),
    median = exp(mu),
    Turnover = rlnorm( n() , mu , sigma )
  ) %>%
  select(starts_with("."), distribution, 
         mu, median, Turnover) %T>%
  print()

PB_kelp_prior_posterior_global %>%
  pivot_longer(cols = -c(starts_with("."), distribution),
               names_to = "parameter") %>%
  group_by(distribution, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n())

# 3.5.2 Species parameters ####
PB_kelp_prior_posterior_species <- PB_kelp_prior %>% 
  prior_posterior_draws(
    posterior_samples = PB_kelp_nc_samples,
    group = PB_kelp %>% select(Species),
    parameters = c("alpha_s[Species]", "alpha_sigma_l",
                   "alpha_sigma_r", "sigma"),
    format = "short"
  ) %>%
  mutate(
    # Calculate mu for observed species, but unobserved levels and references
    mu = rnorm( n() , alpha_s , alpha_sigma_l ) +
          rnorm( n() , 0 , alpha_sigma_r ),
    median = exp(mu),
    Turnover = rlnorm( n() , mu , sigma )
  ) %>% # Remove redundant priors within species
  filter(distribution == "prior" & Species == "Laminaria hyperborea" |
           distribution == "posterior") %>%
  mutate( # Embed prior in species
    Species = if_else(
      distribution == "prior", "Prior", Species
    ) %>% fct()
  ) %>%
  select(starts_with("."), Species, 
         mu, median, Turnover) %T>%
  print()

PB_kelp_prior_posterior_species %>%
  pivot_longer(cols = -c(starts_with("."), Species),
               names_to = "parameter") %>%
  group_by(Species, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n()) %>%
  print(n = 200)

# 3.5.3 Combine parameters ####
PB_kelp_prior_posterior <- PB_kelp_prior_posterior_global %>%
  filter(distribution == "posterior") %>% # priors are redundant
  select(!distribution) %>%
  mutate(Species = "Unobserved" %>% fct()) %>%
  bind_rows(PB_kelp_prior_posterior_species) %T>%
  print()

# Save progress and clean up
PB_kelp_prior_posterior %>%
  write_rds(here("Turnover", "RDS", "PB_kelp_prior_posterior.rds"))

rm(PB_kelp_prior, PB_kelp_c_model, PB_kelp_c_samples,
   PB_kelp_nc_model, PB_kelp_nc_samples)

# 4. Forest turnover ####
# 4.1 Prior simulation ####
# Same as 3.1. I am not aware of an independent source for a prior,
# so am using the same prior as before. Here the only partial pooling
# variable is reference because species and measurement level are not 
# identified, so I can increase prior variance for reference.

# 4.2 Stan models ####
PB_forest_c_model <- here("Turnover", "Stan", "PB_forest_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

PB_forest_nc_model <- here("Turnover", "Stan", "PB_forest_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

PB_forest_c_samples <- PB_forest_c_model$sample(
          data = PB_forest %>% compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

PB_forest_nc_samples <- PB_forest_nc_model$sample(
          data = PB_forest %>% compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

# 4.3 Model checks ####
# Rhat
PB_forest_c_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000718.

PB_forest_nc_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000801.

# Chains
PB_forest_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay()

PB_forest_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Chains are similar

# Pairs
PB_forest_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_mu", "alpha_sigma", "alpha[1]"))
PB_forest_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_mu", "alpha_sigma", "alpha[2]"))

PB_forest_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_mu", "alpha_sigma", "alpha[1]"))
PB_forest_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_mu", "alpha_sigma", "alpha[2]"))
# Lack of correlation is similar

# 4.4 Prior-posterior comparison ####
PB_forest_prior <- prior_samples(
  model = PB_forest_nc_model,
  data = PB_forest %>% compose_data()
  )

PB_forest_prior %>% 
  prior_posterior_draws(
    posterior_samples = PB_forest_c_samples,
    group = PB_forest %>% select(Reference),
    parameters = c("alpha_mu", "alpha_sigma",
                   "alpha[Reference]", "sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Reference"
  )

PB_forest_prior %>% 
  prior_posterior_draws(
    posterior_samples = PB_forest_nc_samples,
    group = PB_forest %>% select(Reference),
    parameters = c("alpha_mu", "alpha_sigma",
                   "alpha[Reference]", "sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Reference"
  )
# Posteriors look similar. Proceed with
# non-centred model.

# 4.5 Prediction ####
# 4.5.1 Global parameters ####
PB_forest_prior_posterior <- PB_forest_prior %>% 
  prior_posterior_draws(
    posterior_samples = PB_forest_nc_samples,
    parameters = c("alpha_mu", "alpha_sigma", "sigma"),
    format = "short"
  ) %>%
  mutate( 
    # Calculate parameters for unobserved references
    mu = rnorm( n() , alpha_mu , alpha_sigma ),
    median = exp(mu),
    Turnover = rlnorm( n() , mu , sigma )
  ) %>%
  select(starts_with("."), distribution, 
         mu, median, Turnover) %T>%
  print()

PB_forest_prior_posterior %>%
  pivot_longer(cols = -c(starts_with("."), distribution),
               names_to = "parameter") %>%
  group_by(distribution, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n())

# Save progress and clean up
PB_forest_prior_posterior %>%
  write_rds(here("Turnover", "RDS", "PB_forest_prior_posterior.rds"))

rm(PB_forest_prior, PB_forest_c_model, PB_forest_c_samples,
   PB_forest_nc_model, PB_forest_nc_samples)

# 5. Kelp k ####
# 5.1 Prior simulation ####
# I am using a normal distribution because the k values reported
# in the literature originate from frequentist nonlinear models
# or were calculated as ( ln(m0) - ln(mt) ) / t and can therefore
# take any value on the real line. In fact several observations 
# report negative k values, which should be impossible according
# to classic exponential decay. To represent this unconstrained
# variation, I am using a model with a normal likelihood. I will 
# partially pool across reviews, references and species. 

# I will take the mean marine macroalgal k from Cebrián & Lartigue 2004
# (doi: 10.1890/03-4019), Cebrián 1999 (doi: 10.1086/303244) etc. as my 
# prior. These data are all contained in object "carbon".
carbon %>%
  filter(Community == "marine macroalgal beds") %>%
  drop_na(`k (d^-1)`) %>%
  summarise(k_min = min(`k (d^-1)`),
            k_median = median(`k (d^-1)`),
            k_max = max(`k (d^-1)`),
            k_mean = mean(`k (d^-1)`),
            k_sd = sd(`k (d^-1)`),
            n = n())
# I will go with 0.01 which is between the median and mean.

tibble(n = 1:1e4,
       alpha_mu = rnorm( 1e4 , 0.01 , 0.02 ),
       alpha_sigma_s = rtnorm( 1e4 , 0 , 0.02 , 0 ), # half-normal distribution
       alpha_sigma_m = rtnorm( 1e4 , 0 , 0.02 , 0 ),
       alpha_sigma_r = rtnorm( 1e4 , 0 , 0.02 , 0 ),
       alpha_s = rnorm( 1e4 , alpha_mu , alpha_sigma_s ),
       alpha_m = rnorm( 1e4 , 0 , alpha_sigma_m ),
       alpha_r = rnorm( 1e4 , 0 , alpha_sigma_r ),
       mu = alpha_s + alpha_m + alpha_r,
       sigma = rexp( 1e4 , 1 ),
       k = rnorm( 1e4 , mu , sigma )) %>%
  pivot_longer(cols = c(mu, k),
               names_to = "parameter") %>%
  ggplot(aes(value)) +
    geom_vline(xintercept = k %$% range(k)) +
    geom_density(alpha = 0.05) +
    scale_x_continuous(limits = c(-1, 1),
                       oob = scales::oob_keep) +
    coord_cartesian(expand = F, clip = "off") +
    facet_wrap(~parameter, scale = "free", nrow = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())
# Covers all possible values.

# 5.2 Stan models ####
k_c_model <- here("Turnover", "Stan", "k_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

k_nc_model <- here("Turnover", "Stan", "k_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

k_c_samples <- k_c_model$sample(
          data = k %>%
            select(k, Species, Review, Reference) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

k_nc_samples <- k_nc_model$sample(
          data = k %>%
            select(k, Species, Review, Reference) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

# 5.3 Model checks ####
# Rhat
k_c_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 94% of rhat above 1.001. rhat = 1.00 ± 0.00108.

k_nc_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000796.

# Chains
k_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay()

k_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Non-centred chains are better

# Pairs
k_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_mu", "alpha_s[1]", "alpha_m[1]", "alpha_r[1]"))
k_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_mu", "alpha_s[2]", "alpha_m[2]", "alpha_r[2]"))
# Some non-identifiability in alpha.

k_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_mu", "alpha_s[1]", "alpha_m[1]", "alpha_r[1]"))
k_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("alpha_mu", "alpha_s[2]", "alpha_m[2]", "alpha_r[2]"))
# Correlations are similar

# 5.4 Prior-posterior comparison ####
k_prior <- prior_samples(
  model = k_nc_model,
  data = k %>%
    select(k, Species, Review, Reference) %>%
    compose_data()
  )

k_prior %>% 
  prior_posterior_draws(
    posterior_samples = k_c_samples,
    group = k %>% select(Species, Review),
    parameters = c("alpha_mu", "alpha_sigma_s",
                   "alpha_s[Species]", 
                   "alpha_sigma_m", "alpha_m[Review]",
                   "sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Species",
    second_group_name = "Review"
  )

k_prior %>% 
  prior_posterior_draws(
    posterior_samples = k_c_samples,
    group = k %>% select(Reference),
    parameters = c("alpha_mu", "alpha_sigma_r",
                   "alpha_r[Reference]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Reference"
  )

k_prior %>% 
  prior_posterior_draws(
    posterior_samples = k_nc_samples,
    group = k %>% select(Species, Review),
    parameters = c("alpha_mu", "alpha_sigma_s",
                   "alpha_s[Species]", 
                   "alpha_sigma_m", "alpha_m[Review]",
                   "sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Species",
    second_group_name = "Review"
  )

k_prior %>% 
  prior_posterior_draws(
    posterior_samples = k_nc_samples,
    group = k %>% select(Reference),
    parameters = c("alpha_mu", "alpha_sigma_r",
                   "alpha_r[Reference]"),
    format = "long"
    ) %>%
  prior_posterior_plot(
    group_name = "Reference"
  )
# Posteriors are smoother for non-centred model.
# Proceed with non-centred model.

# 5.5 Prediction ####
# 5.5.1 Global parameters ####
k_prior_posterior_global <- k_prior %>% 
  prior_posterior_draws(
    posterior_samples = k_nc_samples,
    parameters = c("alpha_mu", "alpha_sigma_s",
                   "alpha_sigma_m", "alpha_sigma_r",
                   "sigma"),
    format = "short"
  ) %>%
  mutate( 
    # Calculate mu for unobserved species, reviews and references
    mu = rnorm( n() , alpha_mu , alpha_sigma_s ) + 
          rnorm( n() , 0 , alpha_sigma_m ) +
          rnorm( n() , 0 , alpha_sigma_r ),
    k = rnorm( n() , mu , sigma )
  ) %>%
  select(starts_with("."), distribution, 
         mu, k) %T>%
  print()

k_prior_posterior_global %>%
  pivot_longer(cols = -c(starts_with("."), distribution),
               names_to = "parameter") %>%
  group_by(distribution, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n())

# 5.5.2 Species parameters ####
k_prior_posterior_species <- k_prior %>% 
  prior_posterior_draws(
    posterior_samples = k_nc_samples,
    group = k %>% select(Species),
    parameters = c("alpha_s[Species]", "alpha_sigma_m",
                   "alpha_sigma_r", "sigma"),
    format = "short"
  ) %>%
  mutate(
    # Calculate mu for observed species, but unobserved reviews and references
    mu = rnorm( n() , alpha_s , alpha_sigma_m ) +
          rnorm( n() , 0 , alpha_sigma_r ),
    k = rnorm( n() , mu , sigma )
  ) %>% # Remove redundant priors within species
  filter(distribution == "prior" & Species == "Macrocystis pyrifera" |
           distribution == "posterior") %>%
  mutate( # Embed prior in species
    Species = if_else(
      distribution == "prior", "Prior", Species
    ) %>% fct()
  ) %>%
  select(starts_with("."), Species, mu, k) %T>%
  print()

k_prior_posterior_species %>%
  pivot_longer(cols = -c(starts_with("."), Species),
               names_to = "parameter") %>%
  group_by(Species, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n()) %>%
  print(n = 30)

# 5.5.3 Combine parameters ####
k_prior_posterior <- k_prior_posterior_global %>%
  filter(distribution == "posterior") %>% # priors are redundant
  select(!distribution) %>%
  mutate(Species = "Unobserved" %>% fct()) %>%
  bind_rows(k_prior_posterior_species) %T>%
  print()
  
# 5.5.4 Summarise ####
k_parameters <- k_prior_posterior %>%
  filter(Species != "Prior") %>%
  select(!starts_with(".")) %>%
  group_by(Species) %>%
  summarise(
    across( everything(), list(mean = mean, sd = sd) ),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(
    mu = glue("{signif(mu_mean, 2)} ± {signif(mu_sd, 2)}"),
    k = glue("{signif(k_mean, 2)} ± {signif(k_sd, 2)}")
  ) %>%
  select(!(contains("mean") | contains("sd"))) %T>%
  print()
# All estimates are way too uncertain to make comparisons.

# Save progress and clean up
k_prior_posterior %>%
  write_rds(here("Turnover", "RDS", "k_prior_posterior.rds"))

rm(k_prior, k_c_model, k_c_samples,
   k_nc_model, k_nc_samples)

# 6. Comparison ####
# 6.1 Merge kelp and forest turnover ####
# Merge kelp and forest turnover posteriors
PB_kelp_forest_posterior <- PB_kelp_prior_posterior %>%
  filter(Species != "Prior") %>%
  droplevels() %>%
  mutate(
    Species = if_else(
      Species == "Unobserved",
      "Kelps", Species 
    ) %>% fct()
  ) %>%
  bind_rows(
    PB_forest_prior_posterior %>%
      filter(distribution == "posterior") %>%
      select(-distribution) %>%
      mutate(Species = "Trees" %>% fct())
  ) %T>%
  print()

# 6.2 Predict decomposition from turnover ####
PB_D_posterior <- PB_prior_posterior_global %>%
  filter(distribution == "posterior") %>%
  select(!distribution) %>%
  full_join(
    PB_kelp_forest_posterior,
    by = c(".chain", ".iteration", ".draw")
  ) %>%
  mutate(
    # Convert turnover rate (y^-1) to turnover time (d)
    Time = 1 / Turnover * 365,
    log_Time = log10(Time),
    log_Turnover = log10(Turnover),
    # Predict log10 P/B of detritus, i.e. decomposition, using relationship
    log_Decomposition = rnorm( n() , alpha + beta * log_Turnover , sigma ),
    # Convert decomposition back to original scale and convert y^-1 to d^-1
    k = 10^log_Decomposition / 365,
    t0.5 = log(2)/k,
    log_k = log10(k),
    log_t0.5 = log10(t0.5)
  ) %T>%
  print()

# 6.3 Create summary table ####
PB_D_parameters <- PB_D_posterior %>%
  select(Species, Turnover, log_Turnover, Time, log_Time, 
         k, log_k, t0.5, log_t0.5) %>%
  group_by(Species) %>%
  summarise(
    across( starts_with("log"), list(mean = mean, sd = sd) ),
    across( everything(), list(median = median) ),
    n = n()
  ) %>%
  ungroup() %>%
  mutate( 
    Turnover = glue("{signif(Turnover_median, 2)} ({signif(log_Turnover_mean, 2)} ± {signif(log_Turnover_sd, 2)})"),
    Time_median_rounded = case_when(
      Time_median < 100 ~ signif(Time_median, 2),
      Time_median < 1e3 ~ signif(Time_median, 3),
      Time_median < 1e4 ~ signif(Time_median, 4)
    ),
    Time = glue("{Time_median_rounded} ({signif(log_Time_mean, 2)} ± {signif(log_Time_sd, 2)})"),
    k = glue("{signif(k_median, 2)} ({signif(log_k_mean, 2)} ± {signif(log_k_sd, 2)})"),
    t0.5_median_rounded = case_when(
      t0.5_median < 100 ~ signif(t0.5_median, 2),
      t0.5_median < 1e3 ~ signif(t0.5_median, 3),
      t0.5_median < 1e4 ~ signif(t0.5_median, 4)
    ),
    t0.5 = glue("{t0.5_median_rounded} ({signif(log_t0.5_mean, 2)} ± {signif(log_t0.5_sd, 2)})")
  ) %>%
  select(!(contains("median") | contains("mean") | contains("sd"))) %T>%
  print(n = 35)

# 6.4 Contrasts ####
PB_D_contrast <- PB_D_posterior %>%
  filter(Species %in% c("Kelps", "Trees")) %>%
  droplevels() %>%
  select(starts_with("."), Species, Turnover, k) %>%
  pivot_longer(cols = c(Turnover, k),
               names_to = "Variable") %>%
  pivot_wider(names_from = Species,
              values_from = value) %>%
  mutate(difference = Kelps - Trees,
         ratio = Kelps / Trees,
         log_ratio = log10(ratio)) %T>%
  print()

# Summarise
PB_D_contrast_summary <- PB_D_contrast %>%
  group_by(Variable) %>%
  summarise(
    across( everything(), list(mean = mean, sd = sd) ),
    P = mean( difference > 0 ) %>% signif(2),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(
    Kelps = glue("{signif(Kelps_mean, 2)} ± {signif(Kelps_sd, 2)}"),
    Trees = glue("{signif(Trees_mean, 2)} ± {signif(Trees_sd, 2)}"),
    difference = glue("{signif(difference_mean, 2)} ± {signif(difference_sd, 2)}"),
    ratio = glue("{signif(ratio_mean, 2)} ± {signif(ratio_sd, 2)}"),
    log_ratio = glue("{signif(log_ratio_mean, 2)} ± {signif(log_ratio_sd, 2)}")
  ) %>%
  select(!(contains("mean") | contains("sd"))) %T>%
  print()
# Probabilities for turnover time and half-life are the exact complements
# to probabilities for turnover and decomposition rate because these
# posteriors are related via constants.

# 6.5 Table S1 ####
Table_S1 <- PB_D_parameters %>%
  mutate(Species = Species %>% 
           fct_relevel(sort) %>%
           fct_relevel("Trees", "Kelps")) %>%
  arrange(Species)

Table_S1 %>%
  write_csv(here("Tables", "Table_S1.csv"))

require(officer)
read_docx() %>%
  body_add_table(value = Table_S1) %>%
  print(target = here("Tables", "Table_S1.docx"))

# 6.6 Figure 1 ####
# 6.6.1 Merge data ####
PB_kelp_forest <- PB_kelp %>%
  mutate(Plant = "Kelps" %>% fct()) %>%
  bind_rows(
    PB_forest %>% 
      mutate(Plant = "Trees" %>% fct())
  ) %T>%
  print()

# 6.6.2 Define custom theme ####
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

# 6.6.3 Figure 1a ####
require(ggridges)
Fig_1a <- PB_D_posterior %>%
  filter(Species %in% c("Kelps", "Trees")) %>%
  droplevels() %>%
  mutate(Species = Species %>%
           fct_recode("Kelp forests" = "Kelps",
                      "Terrestrial forests" = "Trees")) %>%
  ggplot() +
    geom_point(data = PB_kelp_forest,
               aes(Turnover, Plant %>% as.numeric() - 0.5, colour = Plant),
               shape = 16, alpha = 0.2, size = 2.4,
               position = position_jitter(height = 0.3)) +
    stat_density_ridges(aes(Turnover, Species %>% as.numeric(), fill = Species),
                        colour = NA, n = 2^10, from = -3, to = 2,
                        bandwidth = 5*0.02, scale = 1.5, alpha = 0.7) +
    scale_colour_manual(values = c("#a29400", "#004237"), guide = "none") +
    scale_fill_manual(values = c("#a29400", "#004237"), 
                      guide = guide_legend(reverse = TRUE)) +
    scale_x_log10(limits = c(10^-3, 10^2),
                  breaks = 10^(-3:2),
                  labels = scales::label_log(), # scales::label_math(10^.x) is the alternative
                  oob = scales::oob_keep) + # when using the normal scale but does not use minus
    labs(x = expression("Turnover (year"^-1*")")) +
    coord_cartesian(ylim = c(0, NA), expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.margin = margin(0.2, 0.5, 0.2, 0, unit = "cm"))

Fig_1a

# 6.6.3 Figure 1b ####
PB_annotation <- PB_parameters %>%
  mutate(
    label_mu = glue(
      "italic(μ)*' = {alpha} + {beta} × '*italic(x)"
    ),
    label_sigma = glue(
      "italic(σ)*' = {sigma}'"
    )
  ) %>%
  pivot_longer(cols = contains("label"),
               names_to = "parameter",
               values_to = "label",
               names_prefix = "label_") %T>%
  print()

require(geomtextpath)
Fig_1b <- PB_prediction %>%
  filter(Community != "Prior") %>%
  ggplot() +
    geom_textline(data = tibble(x = c(10^-3, 10^3.4), y = c(10^-3, 10^3.4)), aes(x, y),
                  label = "1:1", family = "Futura", size = 3.5, hjust = 1) +
    geom_point(data = PB, aes(PB_plant, PB_detritus),
               shape = 16, alpha = 0.5, size = 2.4) +
    geom_ribbon(data = . %>% filter(Community == "Unobserved"),
                aes(10^log_PB_plant, ymin = 10^log_PB_detritus.lower, 
                    ymax = 10^log_PB_detritus.upper,
                    alpha = factor(.width))) +
    geom_line(data = . %>% filter(Community != "Unobserved"),
              aes(10^log_PB_plant, 10^mu, group = Community),
              linewidth = 0.8) +
    geom_text(data = PB_annotation %>%
                filter(Community == "Unobserved"),
              aes(x = 10^-2.85, 
                  y = c(10^5.5, 10^4.25), 
                  label = label),
              family = "Futura", size = 10, size.unit = "pt", 
              hjust = 0, parse = TRUE) +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_x_log10(limits = c(10^-3, 10^3),
                  breaks = 10^(-3:3),
                  labels = scales::label_log(),
                  oob = scales::oob_keep) +
    scale_y_log10(limits = c(10^-4, 10^6),
                  breaks = 10^(seq(-4, 6, 2)),
                  labels = scales::label_log(),
                  oob = scales::oob_keep) +
    labs(x = expression("Turnover (year"^-1*")"),
         y = expression("Decomposition (year"^-1*")")) +
    coord_cartesian(expand = F, clip = "off") +
    mytheme +
    theme(plot.margin = margin(0.2, 0.5, 0.2, 0, unit = "cm"))

Fig_1b

# 6.6.4 Figure 1c ####
Fig_1c <- PB_D_posterior %>%
  filter(Species %in% c("Kelps", "Trees")) %>%
  droplevels() %>%
  ggplot() +
    stat_density_ridges(aes(k, Species %>% as.numeric(), fill = Species),
                        colour = NA, n = 2^10, from = -8, to = 2,
                        bandwidth = 10*0.02, scale = 1.5, alpha = 0.7) +
    scale_fill_manual(values = c("#a29400", "#004237"), guide = "none") +
    scale_x_log10(limits = c(10^-8, 10^2),
                  breaks = 10^(seq(-8, 2, 2)),
                  labels = scales::label_log(),
                  oob = scales::oob_keep) +
    labs(x = expression("Predicted decomposition (day"^-1*")")) +
    coord_cartesian(expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.margin = margin(0.2, 0.5, 0.2, 0, unit = "cm"))

Fig_1c # the 70% probability of kelp decomposition being
# greater is not obvious. Visualise difference instead.

# Add labels to PB_D_contrast
PB_D_contrast %<>%
  left_join(PB_D_contrast_summary %>%
              select(Variable, P), 
            by = "Variable") %>%
  mutate(label_Kelps = ( P * 100 ) %>% str_c("%"),
         label_Trees = ( (1 - P) * 100 ) %>% str_c("%")) %T>%
  print()

Fig_1c <- PB_D_contrast %>%
  filter(Variable == "k") %>%
  ggplot() +
    stat_density_ridges(aes(10^log_ratio, 0, 
                            fill = if_else(after_stat(x) < 10^0,
                                           "Trees", "Kelps")),
                        geom = "density_ridges_gradient",
                        colour = NA, n = 2^10, from = -4, to = 6,
                        bandwidth = 10*0.02, scale = 1) +
    geom_textdensity(aes(x = 10^log_ratio, y = after_stat(density),
                         label = label_Kelps),
                     colour = "#a29400", family = "Futura",
                     size = 3.5, hjust = 0.6, vjust = 0,
                     n = 2^10, bw = 10*0.02, text_only = TRUE) +
    geom_textdensity(aes(x = 10^log_ratio, y = after_stat(density), 
                         label = label_Trees),
                     colour = "#004237", family = "Futura",
                     size = 3.5, hjust = 0.3, vjust = 0,
                     n = 2^10, bw = 10*0.02, text_only = TRUE) +
    geom_vline(xintercept = 10^0) +
    scale_fill_manual(values = c("#a29400", "#004237"), 
                      guide = "none") +
    scale_x_log10(limits = c(10^-4, 10^6),
                  breaks = 10^(seq(-4, 6, 2)),
                  labels = scales::label_log(),
                  oob = scales::oob_keep) +
    labs(x = "Predicted relative decomposition") +
    coord_cartesian(expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.margin = margin(0.2, 0.5, 0.2, 0, unit = "cm"))

Fig_1c

# 6.6.5 Assemble panels ####
require(patchwork)
Fig_1 <- ( Fig_1a / Fig_1b / Fig_1c ) +
  plot_layout(heights = c(0.5, 1, 0.25))

Fig_1 %>%
  ggsave(filename = "Fig_1.pdf", path = "Figures",
         device = cairo_pdf, height = 15, width = 10, units = "cm")