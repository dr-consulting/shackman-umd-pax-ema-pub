####################################################################################################
# Study 2 Modeling Script: Positive Mood Models

# Description: 
#   The analyses below involve a series of increasingly complex Bayesian multilevel regression 
#   models. The analyses addressed three broad research questions designed to provide a better 
#   understanding of the association between dispositional negativity and momentary 
#   positive affect:
#     1. What is a reasonable estimate of the tonic or "unique" association between dispositional 
#     negativity and momentary positive affect? 
#     2. What is a reasonable estimate of the association between dispositional negativity and 
#     momentary positive affect that can be attributed to differences in overall emotional context? 
#     3. What is a reasonable estimate of the association between dispositional negativity and 
#     momentary positive affect that can be attributed to reactivity to recent emotionally salient 
#     events?

# Modeling Notes: 
#   * Exploratory analyses revealed that positive mood ratings were approximately symmetrical in 
#   their distribution, thus a weakly informative normal prior was chosen for positive mood scores
#   * Missingness was addressed via multiple imputation (see imputation script elsewhere)
#   * To generate a more informative posterior predictive distribution in the missingness models, 
#   we included summary scores of the EMA
#   * POS is mainly an aggregate of positive mood and combines high energy (cheerful) and low 
#   energy (calm) affective states. A multilevel factor analysis supported the combination of these
#   two dimensions of positive momentary mood in a single composite. 
#   * Event ratings were individually mean-centered to maintain separation of between- and within-
#   subjects sources of variation in positive mood
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

#--------------------------------------------------------------------------------------------------
# Location of repo stored locally
wd<-paste0('~/dr-consulting_GH/shackman-umd-pax-ema-pub')
data.folder<-paste0(wd, '/Data')
study2.model<-paste0(wd, '/Study2_model_summaries')

# Will save very large posterior files from analyses (not recommended for git repo)
# For anyone attempting to reproduce these analyses be sure to identify a storage location with sufficient memory
posterior_save_dir <- "/media/matthew/My Book"
study2.out<-paste0(posterior_save_dir, '/EMA_S2_Bayesian_Posteriors')

# Also generally not recommended to store image files on GH... 
study2.graphics<-paste0(posterior_save_dir, '/EMA_S2_Graphics')
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
# Loading data from study 2 - PAX
load(paste0(data.folder, '/study2_data.RData'))

dat.study2_list_stacked <- data.frame() 
for(m in 1:M) # Assumes M is still defined in the imputed data set loaded from .RData above
  dat.study2_list_stacked <- rbind(dat.study2_list_stacked, dat.study2_list[[i]])

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Creating utility functions to gather priors based on imputed data sets
get_imputed_mean <- function(data, variable, logged=TRUE){
  M <- length(data) # data needs to be a list object
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
  M <- length(data) # data needs to be a list object
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
S2_POS_ucm_form <- bf(
  POS ~ 1 + (1|ID)
)+gaussian()

# Creating priors for intercept to help constrain final model in response space
mu <- get_imputed_mean(dat.study2_list, "POS", logged=FALSE)
sigma <- get_imputed_sd(dat.study2_list, "POS", logged=FALSE)

Int_prior <- set_prior(paste0("normal(", mu, ",", sigma, ")"), 
                       class = "Intercept")

# Running model with priors (see above)
S2_POS_ucm <- brm_multiple(S2_POS_ucm_form,
                           data = dat.study2_list,
                           chains = 3,
                           prior = c(Int_prior),
                           iter = 25000,
                           warmup = 20000,
                           control = list(adapt_delta = .99,
                                          max_treedepth = 15))

sink(paste0(study2.model, '/S2_POS_ucm.txt'))
print(summary(S2_POS_ucm), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S2_POS_ucm, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_POS_ucm Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_POS_ucm, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_POS_ucm Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_POS_ucm_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_POS_ucm', 'S2_POS_ucm_form'), 
     file = paste0(study2.out, "/S2_POS_ucm.RData"))
remove(S2_POS_ucm)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Simple model with weakly informative prior - DN as sole predictor 
beta_prior <- set_prior('normal(0, 2)', class='b')

S2_POS_DN_form <- bf(
  POS ~ 1 + c.DN + (1|ID)
)+gaussian()

S2_POS_DN <- brm_multiple(S2_POS_DN_form,
                          data = dat.study2_list,
                          chains = 3,
                          prior = c(Int_prior, 
                                    beta_prior),
                          iter = 25000,
                          warmup = 20000,
                          control = list(adapt_delta = .99,
                                         max_treedepth = 15))

sink(paste0(study2.model, '/S2_POS_DN.txt'))
print(summary(S2_POS_DN), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S2_POS_DN, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_POS_DN Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_POS_DN, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_POS_DN Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_POS_DN_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_POS_DN', 'S2_POS_DN_form'), 
     file = paste0(study2.out, "/S2_POS_DN.RData"))
remove(S2_POS_DN)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Simple model with a weakly informative prior for individual differences in the proportion of 
# negative events reported during the EMA period
S2_POS_prop.NegEvnt_form <- bf(
  POS ~ 1 + prop.NegEvnt + (1|ID)
)+gaussian()

S2_POS_prop.NegEvnt <- brm_multiple(S2_POS_prop.NegEvnt_form,
                                    data = dat.study2_list,
                                    chains = 3,
                                    prior = c(Int_prior, 
                                              beta_prior),
                                    iter = 25000,
                                    warmup = 20000,
                                    control = list(adapt_delta = .99,
                                                   max_treedepth = 15))

sink(paste0(study2.model, '/S2_POS_prop.NegEvnt.txt'))
print(summary(S2_POS_prop.NegEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S2_POS_prop.NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_POS_prop.NegEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_POS_prop.NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_POS_prop.NegEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_POS_prop.NegEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

save(list = c('S2_POS_prop.NegEvnt', 'S2_POS_prop.NegEvnt_form'), 
     file = paste0(study2.out, "/S2_POS_prop.NegEvnt.RData"))
remove(S2_POS_prop.NegEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Simple model with a weakly informative prior for individual differences in the proportion of 
# positive events reported during the EMA period
S2_POS_prop.PosEvnt_form <- bf(
  POS ~ 1 + prop.PosEvnt + (1|ID)
)+gaussian()

S2_POS_prop.PosEvnt <- brm_multiple(S2_POS_prop.PosEvnt_form,
                                    data = dat.study2_list,
                                    chains = 3,
                                    prior = c(Int_prior, 
                                              beta_prior),
                                    iter = 25000,
                                    warmup = 20000,
                                    control = list(adapt_delta = .99,
                                                   max_treedepth = 15))

sink(paste0(study2.model, '/S2_POS_prop.PosEvnt.txt'))
print(summary(S2_POS_prop.PosEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S2_POS_prop.PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_POS_prop.PosEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_POS_prop.PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_POS_prop.PosEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_POS_prop.PosEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_POS_prop.PosEvnt', 'S2_POS_prop.PosEvnt_form'), 
     file = paste0(study2.out, "/S2_POS_prop.PosEvnt.RData"))
remove(S2_POS_prop.PosEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Model that includes two "level 2" predictors - DN and proportion of EMA reports that included the 
# occurence of a negative event. Used in disentangling between-subject sources of variability 
# attributable to DN vs. overall negative emotional context (i.e., variation in the number of 
# negative events reported)
S2_POS_DN_prop.NegEvnt_form <- bf(
  POS ~ 1 + c.DN + prop.NegEvnt + (1|ID)
)+gaussian()

S2_POS_DN_prop.NegEvnt <- brm_multiple(S2_POS_DN_prop.NegEvnt_form,
                                       data = dat.study2_list,
                                       chains = 3,
                                       prior = c(Int_prior, 
                                                 beta_prior),
                                       iter = 25000,
                                       warmup = 20000,
                                       control = list(adapt_delta = .99,
                                                      max_treedepth = 15))

sink(paste0(study2.model, '/S2_POS_DN_prop.NegEvnt.txt'))
print(summary(S2_POS_DN_prop.NegEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S2_POS_DN_prop.NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_POS_DN_prop.NegEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_POS_DN_prop.NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_POS_DN_prop.NegEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_POS_DN_prop.NegEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_POS_DN_prop.NegEvnt', 'S2_POS_DN_prop.NegEvnt_form'), 
     file = paste0(study2.out, "/S2_POS_DN_prop.NegEvnt.RData"))
remove(S2_POS_DN_prop.NegEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Model that includes two "level 2" predictors - DN and proportion of EMA reports that included the 
# occurence of positive event. Used in disentangling between-subject sources of variability 
# attributable to DN vs. overall positive emotional context (i.e., variation in the number of 
# positive events reported)
S2_POS_DN_prop.PosEvnt_form <- bf(
  POS ~ 1 + c.DN + prop.PosEvnt + (1|ID)
)+gaussian()

S2_POS_DN_prop.PosEvnt <- brm_multiple(S2_POS_DN_prop.PosEvnt_form,
                                       data = dat.study2_list,
                                       chains = 3,
                                       prior = c(Int_prior, 
                                                 beta_prior),
                                       iter = 25000,
                                       warmup = 20000,
                                       control = list(adapt_delta = .99,
                                                      max_treedepth = 15))

sink(paste0(study2.model, '/S2_POS_DN_prop.PosEvnt.txt'))
print(summary(S2_POS_DN_prop.PosEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S2_POS_DN_prop.PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_POS_DN_prop.PosEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_POS_DN_prop.PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_POS_DN_prop.PosEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_POS_DN_prop.PosEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_POS_DN_prop.PosEvnt', 'S2_POS_DN_prop.PosEvnt_form'), 
     file = paste0(study2.out, "/S2_POS_DN_prop.PosEvnt.RData"))
remove(S2_POS_DN_prop.PosEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Simple level 1 model with individually centered negative events at "level 1"
# Includes a weakly informative distributional prior for random effects correlations
S2_POS_NegEvnt_form <- bf(
  POS ~ 1 + c.NegEvnt + (1 + c.NegEvnt|ID)
)+gaussian()

cor_prior <- set_prior('lkj(2)', class='cor')

# Running model with priors (see above)
S2_POS_NegEvnt <- brm_multiple(S2_POS_NegEvnt_form,
                               data = dat.study2_list, 
                               chains = 3,
                               prior = c(Int_prior, 
                                         beta_prior, 
                                         cor_prior),
                               iter = 25000,
                               warmup = 20000,
                               control = list(adapt_delta = .99, 
                                              max_treedepth = 15))

sink(paste0(study2.model, '/S2_POS_NegEvnt.txt'))
print(summary(S2_POS_NegEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S2_POS_NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_POS_NegEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_POS_NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_POS_NegEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_POS_NegEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_POS_NegEvnt', 'S2_POS_NegEvnt_form'), 
     file = paste0(study2.out, "/S2_POS_NegEvnt.RData"))
remove(S2_POS_NegEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Simple level 1 model with individually centered positive events at "level 1"
# Includes a weakly informative distributional prior for random effects correlations
S2_POS_PosEvnt_form <- bf(
  POS ~ 1 + c.PosEvnt + (1 + c.PosEvnt|ID)
)+gaussian()

# Running model with priors (see above)
S2_POS_PosEvnt <- brm_multiple(S2_POS_PosEvnt_form,
                               data = dat.study2_list, 
                               chains = 3,
                               prior = c(Int_prior, 
                                         beta_prior, 
                                         cor_prior),
                               iter = 25000,
                               warmup = 20000,
                               control = list(adapt_delta = .99, 
                                              max_treedepth = 15))

sink(paste0(study2.model, '/S2_POS_PosEvnt.txt'))
print(summary(S2_POS_PosEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S2_POS_PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_POS_PosEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_POS_PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_POS_PosEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_POS_PosEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_POS_PosEvnt', 'S2_POS_PosEvnt_form'), 
     file = paste0(study2.out, "/S2_POS_PosEvnt.RData"))
remove(S2_POS_PosEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# A model that includes effects on the intercept at both "levels". Used in the calculations to 
# parse source of variation 
S2_POS_NegEvnt_DN_form <- bf(
  POS ~ 1 + c.NegEvnt + c.DN + (1 + c.NegEvnt|ID)
)+gaussian()

# Running model with priors (see above)
S2_POS_NegEvnt_DN <- brm_multiple(S2_POS_NegEvnt_DN_form,
                                  data = dat.study2_list, 
                                  chains = 3,
                                  prior = c(Int_prior, 
                                            beta_prior, 
                                            cor_prior),
                                  iter = 25000,
                                  warmup = 20000,
                                  control = list(adapt_delta = .99, 
                                                 max_treedepth = 15))

sink(paste0(study2.model, '/S2_POS_NegEvnt_DN.txt'))
print(summary(S2_POS_NegEvnt_DN), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S2_POS_NegEvnt_DN, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_POS_NegEvnt_DN Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_POS_NegEvnt_DN, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_POS_NegEvnt_DN Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_POS_NegEvnt_DN_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_POS_NegEvnt_DN', 'S2_POS_NegEvnt_DN_form'), 
     file = paste0(study2.out, "/S2_POS_NegEvnt_DN.RData"))
remove(S2_POS_NegEvnt_DN)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# A repeat of the model above, but with positive events 
S2_POS_PosEvnt_DN_form <- bf(
  POS ~ 1 + c.PosEvnt + c.DN + (1 + c.PosEvnt|ID)
)+gaussian()

# Running model with priors (see above)
S2_POS_PosEvnt_DN <- brm_multiple(S2_POS_PosEvnt_DN_form,
                                  data = dat.study2_list, 
                                  chains = 3,
                                  prior = c(Int_prior, 
                                            beta_prior, 
                                            cor_prior),
                                  iter = 25000,
                                  warmup = 20000,
                                  control = list(adapt_delta = .99, 
                                                 max_treedepth = 15))

sink(paste0(study2.model, '/S2_POS_PosEvnt_DN.txt'))
print(summary(S2_POS_PosEvnt_DN), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S2_POS_PosEvnt_DN, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_POS_PosEvnt_DN Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_POS_PosEvnt_DN, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_POS_PosEvnt_DN Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_POS_PosEvnt_DN_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_POS_PosEvnt_DN', 'S2_POS_PosEvnt_DN_form'), 
     file = paste0(study2.out, "/S2_POS_PosEvnt_DN.RData"))
remove(S2_POS_PosEvnt_DN)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Swapping proportion of negative events for DN in this "two-level" model. Used to estimate unique
# and shared variance accounted for by proportion of negative events and DN in momentary mood.
S2_POS_NegEvnt_prop.NegEvnt_form <- bf(
  POS ~ 1 + c.NegEvnt + prop.NegEvnt + (1 + c.NegEvnt|ID)
)+gaussian()

# Running model with priors (see above)
S2_POS_NegEvnt_prop.NegEvnt <- brm_multiple(S2_POS_NegEvnt_prop.NegEvnt_form,
                                            data = dat.study2_list, 
                                            chains = 3,
                                            prior = c(Int_prior, 
                                                      beta_prior, 
                                                      cor_prior),
                                            iter = 25000,
                                            warmup = 20000,
                                            control = list(adapt_delta = .99, 
                                                           max_treedepth = 15))

sink(paste0(study2.model, '/S2_POS_NegEvnt_prop.NegEvnt.txt'))
print(summary(S2_POS_NegEvnt_prop.NegEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S2_POS_NegEvnt_prop.NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_POS_NegEvnt_prop.NegEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_POS_NegEvnt_prop.NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_POS_NegEvnt_prop.NegEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_POS_NegEvnt_prop.NegEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_POS_NegEvnt_prop.NegEvnt', 'S2_POS_NegEvnt_prop.NegEvnt_form'), 
     file = paste0(study2.out, "/S2_POS_NegEvnt_prop.NegEvnt.RData"))
remove(S2_POS_NegEvnt_prop.NegEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Swapping proportion of positive events for DN in this "two-level" model. Used to estimate unique
# and shared variance accounted for by proportion of positive events and DN in momentary mood.
S2_POS_PosEvnt_prop.PosEvnt_form <- bf(
  POS ~ 1 + c.PosEvnt + prop.PosEvnt + (1 + c.PosEvnt|ID)
)+gaussian()

# Running model with priors (see above)
S2_POS_PosEvnt_prop.PosEvnt <- brm_multiple(S2_POS_PosEvnt_prop.PosEvnt_form,
                                            data = dat.study2_list, 
                                            chains = 3,
                                            prior = c(Int_prior, 
                                                      beta_prior, 
                                                      cor_prior),
                                            iter = 25000,
                                            warmup = 20000,
                                            control = list(adapt_delta = .99, 
                                                           max_treedepth = 15))

sink(paste0(study2.model, '/S2_POS_PosEvnt_prop.PosEvnt.txt'))
print(summary(S2_POS_PosEvnt_prop.PosEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S2_POS_PosEvnt_prop.PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_POS_PosEvnt_prop.PosEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_POS_PosEvnt_prop.PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_POS_PosEvnt_prop.PosEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_POS_PosEvnt_prop.PosEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_POS_PosEvnt_prop.PosEvnt', 'S2_POS_PosEvnt_prop.PosEvnt_form'), 
     file = paste0(study2.out, "/S2_POS_PosEvnt_prop.PosEvnt.RData"))
remove(S2_POS_PosEvnt_prop.PosEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Model with both DN and prorportion of negative events - used to parse out % of variance accounted 
# for by each of the level 2 variables. 
S2_POS_NegEvnt_DN_prop.NegEvnt_form <- bf(
  POS ~ 1 + c.NegEvnt + c.DN + prop.NegEvnt + (1 + c.NegEvnt|ID)
)+gaussian()

#Running model with priors (see above)
S2_POS_NegEvnt_DN_prop.NegEvnt <- brm_multiple(S2_POS_NegEvnt_DN_prop.NegEvnt_form,
                                               data = dat.study2_list, 
                                               chains = 3,
                                               prior = c(Int_prior, 
                                                         beta_prior, 
                                                         cor_prior),
                                               iter = 25000,
                                               warmup = 20000,
                                               control = list(adapt_delta = .99, 
                                                              max_treedepth = 15))

sink(paste0(study2.model, '/S2_POS_NegEvnt_DN_prop.NegEvnt.txt'))
print(summary(S2_POS_NegEvnt_DN_prop.NegEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S2_POS_NegEvnt_DN_prop.NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_POS_NegEvnt_DN_prop.NegEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_POS_NegEvnt_DN_prop.NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_POS_NegEvnt_DN_prop.NegEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_POS_NegEvnt_DN_prop.NegEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_POS_NegEvnt_DN_prop.NegEvnt', 'S2_POS_NegEvnt_DN_prop.NegEvnt_form'), 
     file = paste0(study2.out, "/S2_POS_NegEvnt_DN_prop.NegEvnt.RData"))
remove(S2_POS_NegEvnt_DN_prop.NegEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Model with both DN and prorportion of positive events - used to parse out % of variance accounted 
# for by each of the level 2 variables. 
S2_POS_PosEvnt_DN_prop.PosEvnt_form <- bf(
  POS ~ 1 + c.PosEvnt + c.DN + prop.PosEvnt + (1 + c.PosEvnt|ID)
)+gaussian()

# Running model with priors (see above)
S2_POS_PosEvnt_DN_prop.PosEvnt <- brm_multiple(S2_POS_PosEvnt_DN_prop.PosEvnt_form,
                                               data = dat.study2_list, 
                                               chains = 3,
                                               prior = c(Int_prior, 
                                                         beta_prior, 
                                                         cor_prior),
                                               iter = 25000,
                                               warmup = 20000,
                                               control = list(adapt_delta = .99, 
                                                              max_treedepth = 15))

sink(paste0(study2.model, '/S2_POS_PosEvnt_DN_prop.PosEvnt.txt'))
print(summary(S2_POS_PosEvnt_DN_prop.PosEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S2_POS_PosEvnt_DN_prop.PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_POS_PosEvnt_DN_prop.PosEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_POS_PosEvnt_DN_prop.PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_POS_PosEvnt_DN_prop.PosEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_POS_PosEvnt_DN_prop.PosEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_POS_PosEvnt_DN_prop.PosEvnt', 'S2_POS_PosEvnt_DN_prop.PosEvnt_form'), 
     file = paste0(study2.out, "/S2_POS_PosEvnt_DN_prop.PosEvnt.RData"))
remove(S2_POS_PosEvnt_DN_prop.PosEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Final model - with reactivity effect included (in theory accounting for level-1 variance)
S2_POS_NegEvnt_x_DN_prop.NegEvnt_form <- bf(
  POS ~ 1 + c.NegEvnt * c.DN + prop.NegEvnt + (1 + c.NegEvnt|ID)
)+gaussian()

# Running model with priors (see above)
S2_POS_NegEvnt_x_DN_prop.NegEvnt <- brm_multiple(S2_POS_NegEvnt_x_DN_prop.NegEvnt_form,
                                                 data = dat.study2_list, 
                                                 chains = 3,
                                                 prior = c(Int_prior, 
                                                           beta_prior, 
                                                           cor_prior),
                                                 iter = 25000,
                                                 warmup = 20000,
                                                 control = list(adapt_delta = .99, 
                                                                max_treedepth = 15))

sink(paste0(study2.model, '/S2_POS_NegEvnt_x_DN_prop.NegEvnt.txt'))
print(summary(S2_POS_NegEvnt_x_DN_prop.NegEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S2_POS_NegEvnt_x_DN_prop.NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_POS_NegEvnt_x_DN_prop.NegEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_POS_NegEvnt_x_DN_prop.NegEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_POS_NegEvnt_x_DN_prop.NegEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_POS_NegEvnt_x_DN_prop.NegEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_POS_NegEvnt_x_DN_prop.NegEvnt', 'S2_POS_NegEvnt_x_DN_prop.NegEvnt_form'), 
     file = paste0(study2.out, "/S2_POS_NegEvnt_x_DN_prop.NegEvnt.RData"))
remove(S2_POS_NegEvnt_DN_prop.NegEvnt)
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Final model - with reactivity effect included (in theory accounting for level-1 variance)
S2_POS_PosEvnt_x_DN_prop.PosEvnt_form <- bf(
  POS ~ 1 + c.PosEvnt * c.DN + prop.PosEvnt + (1 + c.PosEvnt|ID)
)+gaussian()

# Running model with priors (see above)
S2_POS_PosEvnt_x_DN_prop.PosEvnt <- brm_multiple(S2_POS_PosEvnt_x_DN_prop.PosEvnt_form,
                                                 data = dat.study2_list, 
                                                 chains = 3,
                                                 prior = c(Int_prior, 
                                                           beta_prior, 
                                                           cor_prior),
                                                 iter = 25000,
                                                 warmup = 20000, 
                                                 control = list(adapt_delta = .99, 
                                                                max_treedepth = 15))

sink(paste0(study2.model, '/S2_POS_PosEvnt_x_DN_prop.PosEvnt.txt'))
print(summary(S2_POS_PosEvnt_x_DN_prop.PosEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S2_POS_PosEvnt_x_DN_prop.PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 100)+
  ggtitle("S2_POS_PosEvnt_x_DN_prop.PosEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S2_POS_PosEvnt_x_DN_prop.PosEvnt, 
           newdata = dat.study2_list_stacked,
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S2_POS_PosEvnt_x_DN_prop.PosEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study2.graphics, '/S2_POS_PosEvnt_x_DN_prop.PosEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)
dev.off()

# Given size of files saved off posteriors to external hard drive
save(list = c('S2_POS_PosEvnt_x_DN_prop.PosEvnt', 'S2_POS_PosEvnt_x_DN_prop.PosEvnt_form'), 
     file = paste0(study2.out, "/S2_POS_PosEvnt_x_DN_prop.PosEvnt.RData"))
remove(S2_POS_PosEvnt_x_DN_prop.PosEvnt)
gc()
