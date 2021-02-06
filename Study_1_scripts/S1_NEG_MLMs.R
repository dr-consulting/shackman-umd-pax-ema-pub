###################################################################################################
# Study 1 Modeling Script: Negative Mood Models

# Description: 
#   The analyses below involve a series of increasingly complex Bayesian multilevel regression 
#   models. The analyses addressed three broad research questions designed to provide a better 
#   understanding of the association between dispositional negativity and momentary 
#   negative affect:
#     1. What is a reasonable estimate of the tonic or "unique" association between dispositional 
#     negativity and momentary negative affect? 
#     2. What is a reasonable estimate of the association between dispositional negativity and 
#     momentary negative affect that can attributed to differences in overall emotional context? 
#     3. What is a reasonable estimate of the association between dispositional negativity and 
#     momentary negative affect that can be attributed to reactivity to recent emotionally salient 
#     events?

# Modeling Notes: 
#   * Exploratory analyses revealed that negative mood ratings were positively skewed, thus a 
#   weakly informative lognormal prior was chosen for these negative mood scores
#   * Missingness was addressed at runtime by taking draws from the posterior predictive 
#   distribution - for both continuous predictors and momentary negative mood scores
#   * To generate a more informative posterior predictive distribution in the missingness models, 
#   we included summary scores of the EMA - (see: S1_PosEvnt_miss and S1_NegEvnt_miss)
#   * NEG is mainly an aggregate of anxious mood - a limitation addressed in Study 2 
#   * Event ratings were individually mean-centered to maintain separation of between- and within-
#   subjects sources of variation in negative mood
###################################################################################################

#--------------------------------------------------------------------------------------------------
# Loading relevant packages
library(brms)
library(rstan)
library(bayesplot)
library(tidyverse)
library(glue)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
# Location of repo stored locally
wd <- '~/github/ATNL/shackman-umd-pax-ema-pub'
data.folder <- '{wd}/Data' %>% glue()
study1.model <- '{wd}/Study_1_model_summaries' %>% glue()

# Will save very large posterior files from analyses (not recommended for git repo)
# For anyone attempting to reproduce these analyses be sure to identify a storage location with sufficient memory
posterior_save_dir <- "/media/dr-owner/HDD1"
study1.out <- '{posterior_save_dir}/EMA_S1_Bayesian_Posteriors/gamma' %>% glue()

# Also generally not recommended to store image files on GH... 
study1.graphics <- '{study1.out}/diagnostic_plots' %>% glue()
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
# Loading Study 1 Data (from Emotion MS - Shackman et al. 2017)
load('{data.folder}/study1_data.RData' %>% glue())

# source utility functions
source('{wd}/utils.R' %>%  glue())

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# Initial null intercept model - will be necessary to generate final variance estimates
# Note the lognormal() prior for the intercept - due to positive skew of negative mood ratings NEG

prior_config <- c(
  set_prior(
    'lognormal(2, .5)', 
    class = 'shape'
  ), 
  set_prior(
    'lognormal(-1.5, .5)', 
    class = 'sd'
  ), 
  set_prior(
    'lognormal(-.75, .55)', 
    class = "Intercept"
  )
)

S1_NEG_ucm_form <- bf(
  NEG ~ 1 + (1|ID)
) + Gamma(link = "log")

# Running model with priors (see above)
S1_NEG_ucm <- brm(
  S1_NEG_ucm_form,
  data = dat.study1_list[[1]], 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S1_NEG_ucm", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S1_NEG_ucm", "S1_NEG_ucm_form", "dat.study1_list", "dat.study1_model"), 
     file = '{study1.out}/S1_NEG_ucm.RData' %>% glue()) 

sink(paste0(study1.model, '/S1_NEG_ucm.txt'))
print(summary(S1_NEG_ucm), digits = 5)
sink()

create_diagnostic_plots(S1_NEG_ucm, dat.study1_list, n_samples = 100, dir_path = study1.graphics)

remove(list=c("S1_NEG_ucm", "S1_NEG_ucm_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# Model with a level-1 random effect for negative event rating (individually mean centered)
prior_config <- c(
  prior_config, 
  set_prior(
    'normal(0, 2)', 
    class = 'b'
  )
)

S1_NEG_NegEvnt_form <- bf(
  NEG ~ 1 + c.NegEvnt + (1 + c.NegEvnt|ID)
) + Gamma(link = "log")

# Running model with priors (see above)
S1_NEG_NegEvnt <- brm(
  S1_NEG_NegEvnt_form,
  data = dat.study1_list[[1]], 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S1_NEG_NegEvnt", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S1_NEG_NegEvnt", "S1_NEG_NegEvnt_form", "dat.study1_list", "dat.study1_model"), 
     file = '{study1.out}/S1_NEG_NegEvnt.RData' %>% glue()) 

sink(paste0(study1.model, '/S1_NEG_NegEvnt.txt'))
print(summary(S1_NEG_NegEvnt), digits = 5)
sink()

create_diagnostic_plots(S1_NEG_NegEvnt, dat.study1_list, n_samples = 100, dir_path = study1.graphics)

remove(list=c("S1_NEG_NegEvnt", "S1_NEG_NegEvnt_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# Model with a level-1 random effect for positive event rating (individually mean-centered)
S1_NEG_PosEvnt_form <- bf(
  NEG ~ 1 + c.PosEvnt + (1 + c.PosEvnt|ID)
) + Gamma(link = "log")

# Running model with priors (see above)
S1_NEG_PosEvnt <- brm(
  S1_NEG_PosEvnt_form,
  data = dat.study1_list[[1]], 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S1_NEG_PosEvnt", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S1_NEG_PosEvnt", "S1_NEG_PosEvnt_form", "dat.study1_list", "dat.study1_model"), 
     file = '{study1.out}/S1_NEG_PosEvnt.RData' %>% glue()) 

sink(paste0(study1.model, '/S1_NEG_PosEvnt.txt'))
print(summary(S1_NEG_PosEvnt), digits = 5)
sink()

create_diagnostic_plots(S1_NEG_PosEvnt, dat.study1_list, n_samples = 100, dir_path = study1.graphics)

remove(list=c("S1_NEG_PosEvnt", "S1_NEG_PosEvnt_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Model with level 1 random effect for recent negative event ratings and DN 
#   Note that due to the individually mean-centered level-1 effect, this model should capture an 
#   estimate of DN's total between-subjects "effect" on momentary mood. 
S1_NEG_NegEvnt_DN_form <- bf(
  NEG ~ 1 + c.NegEvnt + c.DN + (1 + c.NegEvnt|ID)
) + Gamma(link = "log")

# Running model with priors (see above)
S1_NEG_NegEvnt_DN <- brm(
  S1_NEG_NegEvnt_DN_form,
  data = dat.study1_list[[1]], 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S1_NEG_NegEvnt_DN", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S1_NEG_NegEvnt_DN", "S1_NEG_NegEvnt_DN_form", "dat.study1_list", "dat.study1_model"), 
     file = '{study1.out}/S1_NEG_NegEvnt_DN.RData' %>% glue()) 

sink(paste0(study1.model, '/S1_NEG_NegEvnt_DN.txt'))
print(summary(S1_NEG_NegEvnt_DN), digits = 5)
sink()

create_diagnostic_plots(S1_NEG_NegEvnt_DN, dat.study1_list, n_samples = 100, dir_path = study1.graphics)

remove(list=c("S1_NEG_NegEvnt_DN", "S1_NEG_NegEvnt_DN_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Model with level 1 random effect for recent positive event ratings and DN 
#   Note that due to the individually mean-centered level-1 effect, this model should capture an 
#   estimate of DN's total between-subjects "effect" on momentary mood. Relatedly, accounting for 
#   slight differences related to the Bayesian approach, this effect should be effectively the same
#   as the one observed in the previous model. 
S1_NEG_PosEvnt_DN_form <- bf(
  NEG ~ 1 + c.PosEvnt + c.DN + (1 + c.PosEvnt|ID)
) + Gamma(link = "log")

# Running model with priors (see above)
S1_NEG_PosEvnt_DN <- brm(
  S1_NEG_PosEvnt_DN_form,
  data = dat.study1_list[[1]], 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S1_NEG_PosEvnt_DN", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S1_NEG_PosEvnt_DN", "S1_NEG_PosEvnt_DN_form", "dat.study1_list", "dat.study1_model"), 
     file = '{study1.out}/S1_NEG_PosEvnt_DN.RData' %>% glue()) 

sink(paste0(study1.model, '/S1_NEG_PosEvnt_DN.txt'))
print(summary(S1_NEG_PosEvnt_DN), digits = 5)
sink()

create_diagnostic_plots(S1_NEG_PosEvnt_DN, dat.study1_list, n_samples = 100, dir_path = study1.graphics)

remove(list=c("S1_NEG_PosEvnt_DN", "S1_NEG_PosEvnt_DN_form"))
gc()


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# This model includes *average* negative event ratings as a measure of an individual's overall 
# negative event context. This is our proxy for overall exposure to more intense negative events. 
S1_NEG_NegEvnt_Exp_form <- bf(
  NEG ~ 1 + c.NegEvnt + prop_NegEvnt + (1 + c.NegEvnt|ID)
) + Gamma(link = "log")

# Running model with priors (see above)
S1_NEG_NegEvnt_Exp <- brm(
  S1_NEG_NegEvnt_Exp_form,
  data = dat.study1_list[[1]], 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S1_NEG_NegEvnt_Exp", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S1_NEG_NegEvnt_Exp", "S1_NEG_NegEvnt_Exp_form", "dat.study1_list", "dat.study1_model"), 
     file = '{study1.out}/S1_NEG_NegEvnt_Exp.RData' %>% glue()) 

sink(paste0(study1.model, '/S1_NEG_NegEvnt_Exp.txt'))
print(summary(S1_NEG_NegEvnt_Exp), digits = 5)
sink()

create_diagnostic_plots(S1_NEG_NegEvnt_Exp, dat.study1_list, n_samples = 100, dir_path = study1.graphics)

remove(list=c("S1_NEG_NegEvnt_Exp", "S1_NEG_NegEvnt_Exp_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Model is the same as the one above, with the exception that the focus is on inclusion of between-
# subjects differences in average positive mood rating. The idea is similar though in that this is 
# our proxy for the overall contextual effect of reporting more intense positive events on average.
S1_NEG_PosEvnt_Exp_form <- bf(
  NEG ~ 1 + c.PosEvnt + prop_PosEvnt + (1 + c.PosEvnt|ID)
) + Gamma(link = "log")

# Running model with priors (see above)
S1_NEG_PosEvnt_Exp <- brm(
  S1_NEG_PosEvnt_Exp_form,
  data = dat.study1_list[[1]], 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S1_NEG_PosEvnt_Exp", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S1_NEG_PosEvnt_Exp", "S1_NEG_PosEvnt_Exp_form", "dat.study1_list", "dat.study1_model"), 
     file = '{study1.out}/S1_NEG_PosEvnt_Exp.RData' %>% glue()) 

sink(paste0(study1.model, '/S1_NEG_PosEvnt_Exp.txt'))
print(summary(S1_NEG_PosEvnt_Exp), digits = 5)
sink()

create_diagnostic_plots(S1_NEG_PosEvnt_Exp, dat.study1_list, n_samples = 100, dir_path = study1.graphics)

remove(list=c("S1_NEG_PosEvnt_Exp", "S1_NEG_PosEvnt_Exp_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# This model in conjunction with the models including just DN or just average negative event ratings
# will be used to isolate the amount of "unique" variance attributable to DN, and to overall 
# negative contexts
S1_NEG_NegEvnt_DN_Exp_form <- bf(
  NEG ~ 1 + c.NegEvnt + prop_NegEvnt + c.DN + (1 + c.NegEvnt|ID)
) + Gamma(link = "log")

# Running model with priors (see above)
S1_NEG_NegEvnt_DN_Exp <- brm(
  S1_NEG_NegEvnt_DN_Exp_form,
  data = dat.study1_list[[1]], 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S1_NEG_NegEvnt_DN_Exp", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S1_NEG_NegEvnt_DN_Exp", "S1_NEG_NegEvnt_DN_Exp_form", "dat.study1_list", "dat.study1_model"), 
     file = '{study1.out}/S1_NEG_NegEvnt_DN_Exp.RData' %>% glue()) 

sink(paste0(study1.model, '/S1_NEG_NegEvnt_DN_Exp.txt'))
print(summary(S1_NEG_NegEvnt_DN_Exp), digits = 5)
sink()

create_diagnostic_plots(S1_NEG_NegEvnt_DN_Exp, dat.study1_list, n_samples = 100, dir_path = study1.graphics)

remove(list=c("S1_NEG_NegEvnt_DN_Exp", "S1_NEG_NegEvnt_DN_Exp_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# This model in conjunction with the models including just DN or just average positive event ratings
# will be used to isolate the amount of "unique" variance attributable to DN, and to overall 
# positive contexts
S1_NEG_PosEvnt_DN_Exp_form <- bf(
  NEG ~ 1 + c.PosEvnt + prop_PosEvnt + c.DN + (1 + c.PosEvnt|ID)
) + Gamma(link = "log")

# Running model with priors (see above)
S1_NEG_PosEvnt_DN_Exp <- brm(
  S1_NEG_PosEvnt_DN_Exp_form,
  data = dat.study1_list[[1]], 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S1_NEG_PosEvnt_DN_Exp", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S1_NEG_PosEvnt_DN_Exp", "S1_NEG_PosEvnt_DN_Exp_form", "dat.study1_list", "dat.study1_model"), 
     file = '{study1.out}/S1_NEG_PosEvnt_DN_Exp.RData' %>% glue()) 

sink(paste0(study1.model, '/S1_NEG_PosEvnt_DN_Exp.txt'))
print(summary(S1_NEG_PosEvnt_DN_Exp), digits = 5)
sink()

create_diagnostic_plots(S1_NEG_PosEvnt_DN_Exp, dat.study1_list, n_samples = 100, dir_path = study1.graphics)

remove(list=c("S1_NEG_PosEvnt_DN_Exp", "S1_NEG_PosEvnt_DN_Exp_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# This model actually represents the "final" model for this modeling tree (DN, negative events, 
# negative mood). All variance components obtained in the present set of analyses stem from 
# isolated R2 differences moving from the unconditional model to this model. 
S1_NEG_NegEvnt_Rct_form <- bf(
  NEG ~ 1 + c.NegEvnt * c.DN + prop_NegEvnt + c.DN + (1 + c.NegEvnt|ID)
) + Gamma(link = "log")

# Running model with priors (see above)
S1_NEG_NegEvnt_Rct <- brm(
  S1_NEG_NegEvnt_Rct_form,
  data = dat.study1_list[[1]], 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S1_NEG_NegEvnt_Rct", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S1_NEG_NegEvnt_Rct", "S1_NEG_NegEvnt_Rct_form", "dat.study1_list", "dat.study1_model"), 
     file = '{study1.out}/S1_NEG_NegEvnt_Rct.RData' %>% glue()) 

sink(paste0(study1.model, '/S1_NEG_NegEvnt_Rct.txt'))
print(summary(S1_NEG_NegEvnt_Rct), digits = 5)
sink()

create_diagnostic_plots(S1_NEG_NegEvnt_Rct, dat.study1_list, n_samples = 100, dir_path = study1.graphics)

remove(list=c("S1_NEG_NegEvnt_Rct", "S1_NEG_NegEvnt_Rct_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# This model actually represents the "final" model for this modeling tree (DN, positive events, 
# negative mood). All variance components obtained in the present set of analyses stem from 
# isolated R2 differences moving from the unconditional model to this model. 
S1_NEG_PosEvnt_Rct_form <- bf(
  NEG ~ 1 + c.PosEvnt * c.DN + prop_PosEvnt + c.DN + (1 + c.PosEvnt|ID)
) + Gamma(link = "log")

# Running model with priors (see above)
S1_NEG_PosEvnt_Rct <- brm(
  S1_NEG_PosEvnt_Rct_form,
  data = dat.study1_list[[1]], 
  chains = 3,
  cores = 3,
  iter = 15000,
  warmup = 10000, 
  control = list(adapt_delta = .99, 
                 max_treedepth = 15), 
  seed = 7201969, 
  save_all_pars = TRUE, 
  save_model = "S1_NEG_PosEvnt_Rct", 
  open_progress = FALSE,
  refresh = 0, 
  prior = prior_config
)

save(list = c("S1_NEG_PosEvnt_Rct", "S1_NEG_PosEvnt_Rct_form", "dat.study1_list", "dat.study1_model"), 
     file = '{study1.out}/S1_NEG_PosEvnt_Rct.RData' %>% glue()) 

sink(paste0(study1.model, '/S1_NEG_PosEvnt_Rct.txt'))
print(summary(S1_NEG_PosEvnt_Rct), digits = 5)
sink()

create_diagnostic_plots(S1_NEG_PosEvnt_Rct, dat.study1_list, n_samples = 100, dir_path = study1.graphics)

remove(list=c("S1_NEG_PosEvnt_Rct", "S1_NEG_PosEvnt_Rct_form"))
gc()
