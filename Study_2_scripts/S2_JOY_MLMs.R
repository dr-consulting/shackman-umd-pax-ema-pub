####################################################################################################
# Study 2 Modeling Script: Cheerful Mood Models

# Description: 
#   The analyses below involve a series of increasingly complex Bayesian multilevel regression 
#   models. The analyses addressed three broad research questions designed to provide a better 
#   understanding of the association between dispositional negativity and momentary 
#   cheerful affect:
#     1. What is a reasonable estimate of the tonic or "unique" association between dispositional 
#     negativity and momentary cheerful affect? 
#     2. What is a reasonable estimate of the association between dispositional negativity and 
#     momentary cheerful affect that can be attributed to differences in overall emotional context? 
#     3. What is a reasonable estimate of the association between dispositional negativity and 
#     momentary cheerful affect that can be attributed to reactivity to recent emotionally salient 
#     events?

# Modeling Notes: 
#   * Exploratory analyses revealed that cheerful mood ratings were positively skewed in their
#   distribution, thus a weakly informative lognormal prior was chosen for these cheerful mood 
#   scores
#   * Missingness was addressed via multiple imputation (see imputation script elsewhere)
#   * To generate a more informative posterior predictive distribution in the missingness models, 
#   we included summary scores of the EMA
#   * JOY is mainly an aggregate of cheerful mood - and resembles the main positive mood outcome
#   from Study 1
#   * Event ratings were individually mean-centered to maintain separation of between- and within-
#   subjects sources of variation in cheerful mood
####################################################################################################

#---------------------------------------------------------------------------------------------------
library(brms)
library(rstan)
library(rstanarm)
library(bayesplot)
library(tidyverse)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
user<-ifelse(Sys.getenv("USERPROFILE")=="", "~", Sys.getenv("USERPROFILE"))
wd<-paste0(user, '/Dropbox/UMD/Shackman Lab/EMA_MS')
data.folder<-paste0(wd, '/Data')
study2.out<-paste0(wd, '/Study 2 output')
study2.graphics<-paste0(study2.out, '/Graphics')
study2.model<-paste0(study2.out, '/Model summaries')
stan.code<-paste0(wd, '/Stan_code')
EDA.folder<-paste0(study2.out, '/EDA')
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
# Loading Study 2 data with imputed values
load(paste0(data.folder, '/Study2_Clean_w_Impute.RData'))

dat.study2_list_stacked <- data.frame()
for(m in 1:M)  # Assumes M is still defined in the imputed data environment loaded from .RData
  dat.study2_list_stacked <- rbind(dat.study2_list_stacked, dat.study2_list[[i]])

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Creating utility functions to gather priors based on imputed data sets
get_imputed_mean <- function(data, variable, logged=TRUE){
  M <- length(data)  # data needs to be a list object
  tmp_vec <- c()
  for(m in 1:M){
    var <- unlist(data[[m]][variable])
    if(logged){
      var <- log(var)
    }
    tmp_vec <- c(tmp_vec, mean(var))
  }
  return(mean(tmp_vec))
}

get_imputed_sd <- function(data, variable, logged=TRUE){
  M <- length(data)  # data needs to be a list object
  tmp_vec <- c()
  for(m in 1:M){
    var <- unlist(data[[m]][variable])
    if(logged){
      var <- log(var)
    }
    tmp_vec <- c(tmp_vec, sd(var))
  }
  return(mean(tmp_vec))
}

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Null Model - Intercept Only
S2_JOY_ucm_form <- bf(
  JOY ~ 1 + (1|ID)
)+lognormal()

# Creating priors for intercept to help constrain final model in response space
mu <- get_imputed_mean(dat.study2_list, "JOY")
sigma <- get_imputed_sd(dat.study2_list, "JOY")

Int_prior <- set_prior(paste0("normal(", mu, ",", sigma, ")"),
                       class = "Intercept")

# Running model with priors (see above)
S2_JOY_ucm <- brm_multiple(S2_JOY_ucm_form,
                           data = dat.study2_list,
                           chains = 3,
                           prior = c(Int_prior),
                           iter = 25000,
                           warmup = 20000,
                           control = list(adapt_delta = .99,
                                          max_treedepth = 15))

sink(paste0(study2.model, '/S2_JOY_ucm.txt'))
print(summary(S2_JOY_ucm), digits = 5)
sink()

# Simple model check plotting:
ppc_density <-
  pp_check(S2_JOY_ucm,
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_JOY_ucm Posterior Predictive Distribution")

ppc_hist <-
  pp_check(S2_JOY_ucm,
           newdata = dat.study2_list_stacked,
           nsamples = 12,
           type = "error_hist")+
  ggtitle("S2_JOY_ucm Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_JOY_ucm_ppc.png'),
    units = "in",
    height = 5.5,
    width = 11,
    res = 900)
cowplot::plot_grid(ppc_hist,
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_JOY_ucm', 'S2_JOY_ucm_form'),
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_JOY_ucm.RData')
remove(S2_JOY_ucm)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Simple model with weakly informative prior - DN as sole predictor 
beta_prior <- set_prior('normal(0, 2)', class='b')

S2_JOY_DN_form <- bf(
  JOY ~ 1 + c.DN + (1|ID)
)+lognormal()

S2_JOY_DN <- brm_multiple(S2_JOY_DN_form,
                          data = dat.study2_list,
                          chains = 3,
                          prior = c(Int_prior,
                                    beta_prior),
                          iter = 25000,
                          warmup = 20000,
                          control = list(adapt_delta = .99,
                                         max_treedepth = 15))

sink(paste0(study2.model, '/S2_JOY_DN.txt'))
print(summary(S2_JOY_DN), digits = 5)
sink()

# Simple model check plotting:
ppc_density <-
  pp_check(S2_JOY_DN,
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_JOY_DN Posterior Predictive Distribution")

ppc_hist <-
  pp_check(S2_JOY_DN,
           newdata = dat.study2_list_stacked,
           nsamples = 12,
           type = "error_hist")+
  ggtitle("S2_JOY_DN Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_JOY_DN_ppc.png'),
    units = "in",
    height = 5.5,
    width = 11,
    res = 900)
cowplot::plot_grid(ppc_hist,
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_JOY_DN', 'S2_JOY_DN_form'),
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_JOY_DN.RData')
remove(S2_JOY_DN)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Simple model with a weakly informative prior for individual differences in the proportion of 
# negative events reported during the EMA period
S2_JOY_prop.NegEvnt_form <- bf(
  JOY ~ 1 + prop.NegEvnt + (1|ID)
)+lognormal()

S2_JOY_prop.NegEvnt <- brm_multiple(S2_JOY_prop.NegEvnt_form,
                                    data = dat.study2_list,
                                    chains = 3,
                                    prior = c(Int_prior,
                                              beta_prior),
                                    iter = 25000,
                                    warmup = 20000,
                                    control = list(adapt_delta = .99,
                                                   max_treedepth = 15))

sink(paste0(study2.model, '/S2_JOY_prop.NegEvnt.txt'))
print(summary(S2_JOY_prop.NegEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <-
  pp_check(S2_JOY_prop.NegEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_JOY_prop.NegEvnt Posterior Predictive Distribution")

ppc_hist <-
  pp_check(S2_JOY_prop.NegEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 12,
           type = "error_hist")+
  ggtitle("S2_JOY_prop.NegEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_JOY_prop.NegEvnt_ppc.png'),
    units = "in",
    height = 5.5,
    width = 11,
    res = 900)
cowplot::plot_grid(ppc_hist,
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_JOY_prop.NegEvnt', 'S2_JOY_prop.NegEvnt_form'),
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_JOY_prop.NegEvnt.RData')
remove(S2_JOY_prop.NegEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Simple model with a weakly informative prior for individual differences in the proportion of 
# positive events reported during the EMA period
S2_JOY_prop.PosEvnt_form <- bf(
  JOY ~ 1 + prop.PosEvnt + (1|ID)
)+lognormal()

S2_JOY_prop.PosEvnt <- brm_multiple(S2_JOY_prop.PosEvnt_form,
                                    data = dat.study2_list,
                                    chains = 3,
                                    prior = c(Int_prior,
                                              beta_prior),
                                    iter = 25000,
                                    warmup = 20000,
                                    control = list(adapt_delta = .99,
                                                   max_treedepth = 15))

sink(paste0(study2.model, '/S2_JOY_prop.PosEvnt.txt'))
print(summary(S2_JOY_prop.PosEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <-
  pp_check(S2_JOY_prop.PosEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_JOY_prop.PosEvnt Posterior Predictive Distribution")

ppc_hist <-
  pp_check(S2_JOY_prop.PosEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 12,
           type = "error_hist")+
  ggtitle("S2_JOY_prop.PosEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_JOY_prop.PosEvnt_ppc.png'),
    units = "in",
    height = 5.5,
    width = 11,
    res = 900)
cowplot::plot_grid(ppc_hist,
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_JOY_prop.PosEvnt', 'S2_JOY_prop.PosEvnt_form'),
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_JOY_prop.PosEvnt.RData')
remove(S2_JOY_prop.PosEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Model that includes two "level 2" predictors - DN and proportion of EMA reports that included the 
# occurence of a negative event. Used in disentangling between-subject sources of variability 
# attributable to DN vs. overall negative emotional context (i.e., variation in the number of 
# negative events reported)
S2_JOY_DN_prop.NegEvnt_form <- bf(
  JOY ~ 1 + c.DN + prop.NegEvnt + (1|ID)
)+lognormal()

S2_JOY_DN_prop.NegEvnt <- brm_multiple(S2_JOY_DN_prop.NegEvnt_form,
                                       data = dat.study2_list,
                                       chains = 3,
                                       prior = c(Int_prior,
                                                 beta_prior),
                                       iter = 25000,
                                       warmup = 20000,
                                       control = list(adapt_delta = .99,
                                                      max_treedepth = 15))

sink(paste0(study2.model, '/S2_JOY_DN_prop.NegEvnt.txt'))
print(summary(S2_JOY_DN_prop.NegEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <-
  pp_check(S2_JOY_DN_prop.NegEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_JOY_DN_prop.NegEvnt Posterior Predictive Distribution")

ppc_hist <-
  pp_check(S2_JOY_DN_prop.NegEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 12,
           type = "error_hist")+
  ggtitle("S2_JOY_DN_prop.NegEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_JOY_DN_prop.NegEvnt_ppc.png'),
    units = "in",
    height = 5.5,
    width = 11,
    res = 900)
cowplot::plot_grid(ppc_hist,
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_JOY_DN_prop.NegEvnt', 'S2_JOY_DN_prop.NegEvnt_form'),
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_JOY_DN_prop.NegEvnt.RData')
remove(S2_JOY_DN_prop.NegEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Model that includes two "level 2" predictors - DN and proportion of EMA reports that included the 
# occurence of positive event. Used in disentangling between-subject sources of variability 
# attributable to DN vs. overall positive emotional context (i.e., variation in the number of 
# positive events reported)
S2_JOY_DN_prop.PosEvnt_form <- bf(
  JOY ~ 1 + c.DN + prop.PosEvnt + (1|ID)
)+lognormal()

S2_JOY_DN_prop.PosEvnt <- brm_multiple(S2_JOY_DN_prop.PosEvnt_form,
                                       data = dat.study2_list,
                                       chains = 3,
                                       prior = c(Int_prior,
                                                 beta_prior),
                                       iter = 25000,
                                       warmup = 20000,
                                       control = list(adapt_delta = .99,
                                                      max_treedepth = 15))

sink(paste0(study2.model, '/S2_JOY_DN_prop.PosEvnt.txt'))
print(summary(S2_JOY_DN_prop.PosEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <-
  pp_check(S2_JOY_DN_prop.PosEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_JOY_DN_prop.PosEvnt Posterior Predictive Distribution")

ppc_hist <-
  pp_check(S2_JOY_DN_prop.PosEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 12,
           type = "error_hist")+
  ggtitle("S2_JOY_DN_prop.PosEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_JOY_DN_prop.PosEvnt_ppc.png'),
    units = "in",
    height = 5.5,
    width = 11,
    res = 900)
cowplot::plot_grid(ppc_hist,
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_JOY_DN_prop.PosEvnt', 'S2_JOY_DN_prop.PosEvnt_form'),
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_JOY_DN_prop.PosEvnt.RData')
remove(S2_JOY_DN_prop.PosEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Simple level 1 model with individually centered negative events at "level 1"
# Includes a weakly informative distributional prior for random effects correlations
S2_JOY_NegEvnt_form <- bf(
  JOY ~ 1 + c.NegEvnt + (1 + c.NegEvnt|ID)
)+lognormal()

cor_prior <- set_prior('lkj(2)', class='cor')

# Running model with priors (see above)
S2_JOY_NegEvnt <- brm_multiple(S2_JOY_NegEvnt_form,
                               data = dat.study2_list,
                               chains = 3,
                               prior = c(Int_prior,
                                         beta_prior,
                                         cor_prior),
                               iter = 25000,
                               warmup = 20000,
                               control = list(adapt_delta = .99,
                                              max_treedepth = 15))

sink(paste0(study2.model, '/S2_JOY_NegEvnt.txt'))
print(summary(S2_JOY_NegEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <-
  pp_check(S2_JOY_NegEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_JOY_NegEvnt Posterior Predictive Distribution")

ppc_hist <-
  pp_check(S2_JOY_NegEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 12,
           type = "error_hist")+
  ggtitle("S2_JOY_NegEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_JOY_NegEvnt_ppc.png'),
    units = "in",
    height = 5.5,
    width = 11,
    res = 900)
cowplot::plot_grid(ppc_hist,
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_JOY_NegEvnt', 'S2_JOY_NegEvnt_form'),
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_JOY_NegEvnt.RData')
remove(S2_JOY_NegEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Simple level 1 model with individually centered positive events at "level 1"
# Includes a weakly informative distributional prior for random effects correlations
S2_JOY_PosEvnt_form <- bf(
  JOY ~ 1 + c.PosEvnt + (1 + c.PosEvnt|ID)
)+lognormal()

# Running model with priors (see above)
S2_JOY_PosEvnt <- brm_multiple(S2_JOY_PosEvnt_form,
                               data = dat.study2_list,
                               chains = 3,
                               prior = c(Int_prior,
                                         beta_prior,
                                         cor_prior),
                               iter = 25000,
                               warmup = 20000,
                               control = list(adapt_delta = .99,
                                              max_treedepth = 15))

sink(paste0(study2.model, '/S2_JOY_PosEvnt.txt'))
print(summary(S2_JOY_PosEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <-
  pp_check(S2_JOY_PosEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_JOY_PosEvnt Posterior Predictive Distribution")

ppc_hist <-
  pp_check(S2_JOY_PosEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 12,
           type = "error_hist")+
  ggtitle("S2_JOY_PosEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_JOY_PosEvnt_ppc.png'),
    units = "in",
    height = 5.5,
    width = 11,
    res = 900)
cowplot::plot_grid(ppc_hist,
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_JOY_PosEvnt', 'S2_JOY_PosEvnt_form'),
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_JOY_PosEvnt.RData')
remove(S2_JOY_PosEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# A model that includes effects on the intercept at both "levels". Used in the calculations to 
# parse source of variation 
S2_JOY_NegEvnt_DN_form <- bf(
  JOY ~ 1 + c.NegEvnt + c.DN + (1 + c.NegEvnt|ID)
)+lognormal()

# Running model with priors (see above)
S2_JOY_NegEvnt_DN <- brm_multiple(S2_JOY_NegEvnt_DN_form,
                                 data = dat.study2_list,
                                 chains = 3,
                                 prior = c(Int_prior,
                                           beta_prior,
                                           cor_prior),
                                 iter = 25000,
                                 warmup = 20000,
                                 control = list(adapt_delta = .99,
                                                max_treedepth = 15))

sink(paste0(study2.model, '/S2_JOY_NegEvnt_DN.txt'))
print(summary(S2_JOY_NegEvnt_DN), digits = 5)
sink()

# Simple model check plotting:
ppc_density <-
  pp_check(S2_JOY_NegEvnt_DN,
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_JOY_NegEvnt_DN Posterior Predictive Distribution")

ppc_hist <-
  pp_check(S2_JOY_NegEvnt_DN,
           newdata = dat.study2_list_stacked,
           nsamples = 12,
           type = "error_hist")+
  ggtitle("S2_JOY_NegEvnt_DN Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_JOY_NegEvnt_DN_ppc.png'),
    units = "in",
    height = 5.5,
    width = 11,
    res = 900)
cowplot::plot_grid(ppc_hist,
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_JOY_NegEvnt_DN', 'S2_JOY_NegEvnt_DN_form'),
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_JOY_NegEvnt_DN.RData')
remove(S2_JOY_NegEvnt_DN)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# A repeat of the model above, but with positive events 
S2_JOY_PosEvnt_DN_form <- bf(
  JOY ~ 1 + c.PosEvnt + c.DN + (1 + c.PosEvnt|ID)
)+lognormal()

# Running model with priors (see above)
S2_JOY_PosEvnt_DN <- brm_multiple(S2_JOY_PosEvnt_DN_form,
                                 data = dat.study2_list,
                                 chains = 3,
                                 prior = c(Int_prior,
                                           beta_prior,
                                           cor_prior),
                                 iter = 25000,
                                 warmup = 20000,
                                 control = list(adapt_delta = .99,
                                                max_treedepth = 15))

sink(paste0(study2.model, '/S2_JOY_PosEvnt_DN.txt'))
print(summary(S2_JOY_PosEvnt_DN), digits = 5)
sink()

# Simple model check plotting:
ppc_density <-
  pp_check(S2_JOY_PosEvnt_DN,
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_JOY_PosEvnt_DN Posterior Predictive Distribution")

ppc_hist <-
  pp_check(S2_JOY_PosEvnt_DN,
           newdata = dat.study2_list_stacked,
           nsamples = 12,
           type = "error_hist")+
  ggtitle("S2_JOY_PosEvnt_DN Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_JOY_PosEvnt_DN_ppc.png'),
    units = "in",
    height = 5.5,
    width = 11,
    res = 900)
cowplot::plot_grid(ppc_hist,
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_JOY_PosEvnt_DN', 'S2_JOY_PosEvnt_DN_form'),
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_JOY_PosEvnt_DN.RData')
remove(S2_JOY_PosEvnt_DN)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Swapping proportion of negative events for DN in this "two-level" model. Used to estimate unique
# and shared variance accounted for by proportion of negative events and DN in momentary mood.
S2_JOY_NegEvnt_prop.NegEvnt_form <- bf(
  JOY ~ 1 + c.NegEvnt + prop.NegEvnt + (1 + c.NegEvnt|ID)
)+lognormal()

# Running model with priors (see above)
S2_JOY_NegEvnt_prop.NegEvnt <- brm_multiple(S2_JOY_NegEvnt_prop.NegEvnt_form,
                                  data = dat.study2_list,
                                  chains = 3,
                                  prior = c(Int_prior,
                                            beta_prior,
                                            cor_prior),
                                  iter = 25000,
                                  warmup = 20000,
                                  control = list(adapt_delta = .99,
                                                 max_treedepth = 15))

sink(paste0(study2.model, '/S2_JOY_NegEvnt_prop.NegEvnt.txt'))
print(summary(S2_JOY_NegEvnt_prop.NegEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <-
  pp_check(S2_JOY_NegEvnt_prop.NegEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_JOY_NegEvnt_prop.NegEvnt Posterior Predictive Distribution")

ppc_hist <-
  pp_check(S2_JOY_NegEvnt_prop.NegEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 12,
           type = "error_hist")+
  ggtitle("S2_JOY_NegEvnt_prop.NegEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_JOY_NegEvnt_prop.NegEvnt_ppc.png'),
    units = "in",
    height = 5.5,
    width = 11,
    res = 900)
cowplot::plot_grid(ppc_hist,
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_JOY_NegEvnt_prop.NegEvnt', 'S2_JOY_NegEvnt_prop.NegEvnt_form'),
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_JOY_NegEvnt_prop.NegEvnt.RData')
remove(S2_JOY_NegEvnt_prop.NegEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Swapping proportion of positive events for DN in this "two-level" model. Used to estimate unique
# and shared variance accounted for by proportion of positive events and DN in momentary mood.
S2_JOY_PosEvnt_prop.PosEvnt_form <- bf(
  JOY ~ 1 + c.PosEvnt + prop.PosEvnt + (1 + c.PosEvnt|ID)
)+lognormal()

# Running model with priors (see above)
S2_JOY_PosEvnt_prop.PosEvnt <- brm_multiple(S2_JOY_PosEvnt_prop.PosEvnt_form,
                                  data = dat.study2_list,
                                  chains = 3,
                                  prior = c(Int_prior,
                                            beta_prior,
                                            cor_prior),
                                  iter = 25000,
                                  warmup = 20000,
                                  control = list(adapt_delta = .99,
                                                 max_treedepth = 15))

sink(paste0(study2.model, '/S2_JOY_PosEvnt_prop.PosEvnt.txt'))
print(summary(S2_JOY_PosEvnt_prop.PosEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <-
  pp_check(S2_JOY_PosEvnt_prop.PosEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_JOY_PosEvnt_prop.PosEvnt Posterior Predictive Distribution")

ppc_hist <-
  pp_check(S2_JOY_PosEvnt_prop.PosEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 12,
           type = "error_hist")+
  ggtitle("S2_JOY_PosEvnt_prop.PosEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_JOY_PosEvnt_prop.PosEvnt_ppc.png'),
    units = "in",
    height = 5.5,
    width = 11,
    res = 900)
cowplot::plot_grid(ppc_hist,
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_JOY_PosEvnt_prop.PosEvnt', 'S2_JOY_PosEvnt_prop.PosEvnt_form'),
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_JOY_PosEvnt_prop.PosEvnt.RData')
remove(S2_JOY_PosEvnt_prop.PosEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Model with both DN and prorportion of negative events - used to parse out % of variance accounted 
# for by each of the level 2 variables. 
S2_JOY_NegEvnt_DN_prop.NegEvnt_form <- bf(
  JOY ~ 1 + c.NegEvnt + c.DN + prop.NegEvnt + (1 + c.NegEvnt|ID)
)+lognormal()

# Running model with priors (see above)
S2_JOY_NegEvnt_DN_prop.NegEvnt <- brm_multiple(S2_JOY_NegEvnt_DN_prop.NegEvnt_form,
                                            data = dat.study2_list,
                                            chains = 3,
                                            prior = c(Int_prior,
                                                      beta_prior,
                                                      cor_prior),
                                            iter = 25000,
                                            warmup = 20000,
                                            control = list(adapt_delta = .99,
                                                           max_treedepth = 15))

sink(paste0(study2.model, '/S2_JOY_NegEvnt_DN_prop.NegEvnt.txt'))
print(summary(S2_JOY_NegEvnt_DN_prop.NegEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <-
  pp_check(S2_JOY_NegEvnt_DN_prop.NegEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_JOY_NegEvnt_DN_prop.NegEvnt Posterior Predictive Distribution")

ppc_hist <-
  pp_check(S2_JOY_NegEvnt_DN_prop.NegEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 12,
           type = "error_hist")+
  ggtitle("S2_JOY_NegEvnt_DN_prop.NegEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_JOY_NegEvnt_DN_prop.NegEvnt_ppc.png'),
    units = "in",
    height = 5.5,
    width = 11,
    res = 900)
cowplot::plot_grid(ppc_hist,
                   ppc_density)
dev.off()

save(list = c('S2_JOY_NegEvnt_DN_prop.NegEvnt', 'S2_JOY_NegEvnt_DN_prop.NegEvnt_form'),
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_JOY_NegEvnt_DN_prop.NegEvnt.RData')
remove(S2_JOY_NegEvnt_DN_prop.NegEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Model with both DN and prorportion of positive events - used to parse out % of variance accounted 
# for by each of the level 2 variables. 
S2_JOY_PosEvnt_DN_prop.PosEvnt_form <- bf(
  JOY ~ 1 + c.PosEvnt + c.DN + prop.PosEvnt + (1 + c.PosEvnt|ID)
)+lognormal()

#Running model with priors (see above)
S2_JOY_PosEvnt_DN_prop.PosEvnt <- brm_multiple(S2_JOY_PosEvnt_DN_prop.PosEvnt_form,
                                            data = dat.study2_list,
                                            chains = 3,
                                            prior = c(Int_prior,
                                                      beta_prior,
                                                      cor_prior),
                                            iter = 25000,
                                            warmup = 20000,
                                            control = list(adapt_delta = .99,
                                                           max_treedepth = 15))

sink(paste0(study2.model, '/S2_JOY_PosEvnt_DN_prop.PosEvnt.txt'))
print(summary(S2_JOY_PosEvnt_DN_prop.PosEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <-
  pp_check(S2_JOY_PosEvnt_DN_prop.PosEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_JOY_PosEvnt_DN_prop.PosEvnt Posterior Predictive Distribution")

ppc_hist <-
  pp_check(S2_JOY_PosEvnt_DN_prop.PosEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 12,
           type = "error_hist")+
  ggtitle("S2_JOY_PosEvnt_DN_prop.PosEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_JOY_PosEvnt_DN_prop.PosEvnt_ppc.png'),
    units = "in",
    height = 5.5,
    width = 11,
    res = 900)
cowplot::plot_grid(ppc_hist,
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_JOY_PosEvnt_DN_prop.PosEvnt', 'S2_JOY_PosEvnt_DN_prop.PosEvnt_form'),
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_JOY_PosEvnt_DN_prop.PosEvnt.RData')
remove(S2_JOY_PosEvnt_DN_prop.PosEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Final model - with reactivity effect included (in theory accounting for level-1 variance)
S2_JOY_NegEvnt_x_DN_prop.NegEvnt_form <- bf(
  JOY ~ 1 + c.NegEvnt * c.DN + prop.NegEvnt + (1 + c.NegEvnt|ID)
)+lognormal()

# Running model with priors (see above)
S2_JOY_NegEvnt_x_DN_prop.NegEvnt <- brm_multiple(S2_JOY_NegEvnt_x_DN_prop.NegEvnt_form,
                                            data = dat.study2_list,
                                            chains = 3,
                                            prior = c(Int_prior,
                                                      beta_prior,
                                                      cor_prior),
                                            iter = 25000,
                                            warmup = 20000,
                                            control = list(adapt_delta = .99,
                                                           max_treedepth = 15))

sink(paste0(study2.model, '/S2_JOY_NegEvnt_x_DN_prop.NegEvnt.txt'))
print(summary(S2_JOY_NegEvnt_x_DN_prop.NegEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <-
  pp_check(S2_JOY_NegEvnt_x_DN_prop.NegEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_JOY_NegEvnt_x_DN_prop.NegEvnt Posterior Predictive Distribution")

ppc_hist <-
  pp_check(S2_JOY_NegEvnt_x_DN_prop.NegEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 12,
           type = "error_hist")+
  ggtitle("S2_JOY_NegEvnt_x_DN_prop.NegEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_JOY_NegEvnt_x_DN_prop.NegEvnt_ppc.png'),
    units = "in",
    height = 5.5,
    width = 11,
    res = 900)
cowplot::plot_grid(ppc_hist,
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_JOY_NegEvnt_x_DN_prop.NegEvnt', 'S2_JOY_NegEvnt_x_DN_prop.NegEvnt_form'),
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_JOY_NegEvnt_x_DN_prop.NegEvnt.RData')
remove(S2_JOY_NegEvnt_DN_prop.NegEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Final model - with reactivity effect included (in theory accounting for level-1 variance)
S2_JOY_PosEvnt_x_DN_prop.PosEvnt_form <- bf(
  JOY ~ 1 + c.PosEvnt * c.DN + prop.PosEvnt + (1 + c.PosEvnt|ID)
)+lognormal()

#Running model with priors (see above)
S2_JOY_PosEvnt_x_DN_prop.PosEvnt <- brm_multiple(S2_JOY_PosEvnt_x_DN_prop.PosEvnt_form,
                                               data = dat.study2_list,
                                               chains = 3,
                                               prior = c(Int_prior,
                                                         beta_prior,
                                                         cor_prior),
                                               iter = 25000,
                                               warmup = 20000,
                                               control = list(adapt_delta = .99,
                                                              max_treedepth = 15))

sink(paste0(study2.model, '/S2_JOY_PosEvnt_x_DN_prop.PosEvnt.txt'))
print(summary(S2_JOY_PosEvnt_x_DN_prop.PosEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <-
  pp_check(S2_JOY_PosEvnt_x_DN_prop.PosEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_JOY_PosEvnt_x_DN_prop.PosEvnt Posterior Predictive Distribution")

ppc_hist <-
  pp_check(S2_JOY_PosEvnt_x_DN_prop.PosEvnt,
           newdata = dat.study2_list_stacked,
           nsamples = 12,
           type = "error_hist")+
  ggtitle("S2_JOY_PosEvnt_x_DN_prop.PosEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_JOY_PosEvnt_x_DN_prop.PosEvnt_ppc.png'),
    units = "in",
    height = 5.5,
    width = 11,
    res = 900)
cowplot::plot_grid(ppc_hist,
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_JOY_PosEvnt_x_DN_prop.PosEvnt', 'S2_JOY_PosEvnt_x_DN_prop.PosEvnt_form'),
     file = '/media/matthew/My Book/EMA_S2_Bayesian_Posteriors/S2_JOY_PosEvnt_x_DN_prop.PosEvnt.RData')
remove(S2_JOY_PosEvnt_x_DN_prop.PosEvnt)
gc()
