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

prior_config <- c(
    prior_config, 
    set_prior(
        'normal(0, 2)', 
        class = 'b'
    )
)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# This model includes *average* negative event ratings as a measure of an individual's overall 
# negative event context. This is our proxy for overall exposure to more intense negative events. 
S1_NEG_NegEvnt_Exp_form <- bf(
    NEG ~ 1 + c.NegEvnt + prop_NegEvnt + (1 + c.NegEvnt|ID)
) + Gamma(link = "log")

# Running model with priors (see above)
S1_NEG_NegEvnt_Exp <- brm_multiple(
    S1_NEG_NegEvnt_Exp_form,
    data = dat.study1_list, 
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
S1_NEG_PosEvnt_Exp <- brm_multiple(
    S1_NEG_PosEvnt_Exp_form,
    data = dat.study1_list, 
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
S1_NEG_NegEvnt_DN_Exp <- brm_multiple(
    S1_NEG_NegEvnt_DN_Exp_form,
    data = dat.study1_list, 
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
S1_NEG_PosEvnt_DN_Exp <- brm_multiple(
    S1_NEG_PosEvnt_DN_Exp_form,
    data = dat.study1_list, 
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
S1_NEG_NegEvnt_Rct <- brm_multiple(
    S1_NEG_NegEvnt_Rct_form,
    data = dat.study1_list, 
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
S1_NEG_PosEvnt_Rct <- brm_multiple(
    S1_NEG_PosEvnt_Rct_form,
    data = dat.study1_list, 
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
