####################################################################################################
# Study 1 Modeling Script: Positive Mood Models

# Description: 
#   The analyses below involve a series of increasingly complex Bayesian multilevel regression 
#   models. The analyses addressed three broad research questions designed to provide a better 
#   understanding of the association between dispositional negativity and momentary 
#   positive affect:
#     1. What is a reasonable estimate of the tonic or "unique" association between dispositional 
#     negativity and momentary positive affect? 
#     2. What is a reasonable estimate of the association between dispositional negativity and 
#     momentary positive affect that can attributed to differences in overall emotional context? 
#     3. What is a reasonable estimate of the association between dispositional negativity and 
#     momentary positive affect that can be attributed to reactivity to recent emotionally salient 
#     events?

# Modeling Notes: 
#   * Exploratory analyses revealed that positive mood ratings were reasonably symmetrical in their
#   distribution, thus a weakly informative normal prior was chosen for these positive mood scores
#   * Missingness was addressed at runtime by taking draws from the posterior predictive 
#   distribution - for both continuous predictors and momentary positive mood scores
#   * To generate a more informative posterior predictive distribution in the missingness models, 
#   we included summary scores of the EMA - (see: S1_PosEvnt_miss and S1_NegEvnt_miss)
#   * POS is mainly an aggregate of cheerful mood - a limitation addressed in Study 2 
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

#---------------------------------------------------------------------------------------------------
user<-ifelse(Sys.getenv("USERPROFILE")=="", "~", Sys.getenv("USERPROFILE"))
wd<-paste0(user, '/Dropbox/UMD/Shackman Lab/EMA_MS')
data.folder<-paste0(wd, '/Data')
study1.out<-paste0(wd, '/Study 1 output')
study1.graphics<-paste0(study1.out, '/Graphics')
study1.model<-paste0(study1.out, '/Model summaries')
stan.code<-paste0(wd, '/Stan_code')
EDA.folder<-paste0(study1.out, '/EDA')
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
#Loading Study 1 Data (from Emotion MS - Shackman et al. 2017)
load(paste0(data.folder, '/Emotion MS environment.RData'))

dat.study1_lv1 <- dat %>% 
  transmute(ID = subid,
            NEG = MAFS_NA, 
            POS = MAFS_PA,
            orig.DN = ZAP_Both,
            NegEvnt = WorstEvent_Neg,
            PosEvnt = BestEvent_Pos,
            SocMotvtn = WantToBeWithPeople, 
            AlnMotvtn = WantToBeAlone, 
            Probe = ProbeNum, 
            Chrfl = PA_Cheerful, 
            Hppy = PA_Happy, 
            Jyfl = PA_Joyful, 
            Nrvs = NA_Nervous, 
            Anxs = NA_Anxious, 
            Unesy = NA_Uneasy)

dat.study1_lv2 <- dat.study1_lv1 %>% 
  group_by(ID) %>% 
  summarize(c.DN = mean(orig.DN),
            m.NegEvnt = mean(NegEvnt, na.rm = TRUE), 
            m.PosEvnt = mean(PosEvnt, na.rm = TRUE), 
            sd.NegEvnt = sd(NegEvnt, na.rm = TRUE), 
            sd.PosEvnt = sd(PosEvnt, na.rm = TRUE),
            m.NEG = mean(NEG, na.rm = TRUE), 
            m.POS = mean(POS, na.rm = TRUE), 
            sd.NEG = sd(NEG, na.rm = TRUE), 
            sd.POS = sd(POS, na.rm = TRUE)) 

dat.study1_model <- merge(dat.study1_lv1, dat.study1_lv2, by = "ID")
dat.study1_model$c.NegEvnt <- dat.study1_model$NegEvnt - dat.study1_model$m.NegEvnt
dat.study1_model$c.PosEvnt <- dat.study1_model$PosEvnt - dat.study1_model$m.PosEvnt

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# Missing data models - incorporating information about individual distributions of EMA data
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
S1_PosEvnt_miss <- bf(
  c.PosEvnt | mi() ~ 1 + m.NegEvnt + m.PosEvnt + sd.NegEvnt + sd.PosEvnt + 
    m.NEG + sd.NEG + m.POS + sd.POS + (1|ID)
) + gaussian()

S1_NegEvnt_miss <- bf(
  c.NegEvnt | mi() ~ 1 + m.NegEvnt + m.PosEvnt + sd.NegEvnt + sd.PosEvnt + 
    m.NEG + sd.NEG + m.POS + sd.POS + (1|ID)
) + gaussian()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# Initial null intercept model - will be necessary to generate final variance estimates
# Note the guassian() prior for the intercept
S1_POS_ucm_form <- bf(
  POS | mi() ~ 1 + (1|ID)
)+gaussian()

# Creating priors for intercept to help constrain final model in response space
mu <- mean(dat.study1_model$POS, na.rm=TRUE)
sigma <- sd(dat.study1_model$POS, na.rm = TRUE)

Int_prior <- set_prior(paste0("normal(", mu, ",", sigma, ")"), 
                       class = "Intercept")

# Running model with priors (see above)
S1_POS_ucm <- brm(S1_POS_ucm_form,
                  data = dat.study1_model, 
                  chains = 3,
                  prior = c(Int_prior),
                  iter = 15000,
                  warmup = 10000, 
                  control = list(adapt_delta = .99, 
                                 max_treedepth = 15))

sink(paste0(study1.model, '/S1_POS_ucm.txt'))
print(summary(S1_POS_ucm), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S1_POS_ucm, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 100)+
  ggtitle("S1_POS_ucm Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_POS_ucm, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S1_POS_ucm Model Posterior Residuals")

# Saving Plots:
png(paste0(study1.graphics, '/S1_POS_ucm_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()

save(file=paste0(study1.out, "/Posteriors/S1_POS_ucm.RData"), 
     list=c("S1_POS_ucm", "S1_POS_ucm_form"))
remove(list=c("S1_POS_ucm", "S1_POS_ucm_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# Model with a level-1 random effect for negative event rating (indiviudally mean centered)
S1_POS_NegEvnt_form <- bf(
  POS | mi() ~ 1 + mi(c.NegEvnt) + (1 + mi(c.NegEvnt)|ID)
)+gaussian()

# Need to update intercept prior here as there are technically now "multiple" responses 
Int_prior <- set_prior(paste0("normal(", mu, ",", sigma, ")"), 
                       class = "Intercept", 
                       resp = "POS")

# Running model with priors (see above)
S1_POS_NegEvnt <- brm(S1_POS_NegEvnt_form + 
                        S1_NegEvnt_miss + 
                        set_rescor(rescor = FALSE),
                      data = dat.study1_model, 
                      chains = 3,
                      prior = c(Int_prior),
                      iter = 15000,
                      warmup = 10000, 
                      control = list(adapt_delta = .99, 
                                     max_treedepth = 15))

sink(paste0(study1.model, '/S1_POS_NegEvnt.txt'))
print(summary(S1_POS_NegEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S1_POS_NegEvnt, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 100, 
           resp = "POS")+
  ggtitle("S1_POS_NegEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_POS_NegEvnt, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 12, 
           type = "error_hist", 
           resp = "POS")+
  ggtitle("S1_POS_NegEvnt Model Posterior Residuals")

# Saving plots:
png(paste0(study1.graphics, '/S1_POS_NegEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()

save(file=paste0(study1.out, "/Posteriors/S1_POS_NegEvnt.RData"), 
     list=c("S1_POS_NegEvnt", "S1_POS_NegEvnt_form"))
remove(list=c("S1_POS_NegEvnt", "S1_POS_NegEvnt_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# Model with a level-1 random effect for positive event rating (individually mean-centered)
S1_POS_PosEvnt_form <- bf(
  POS | mi() ~ 1 + mi(c.PosEvnt) + (1 + mi(c.PosEvnt)|ID)
)+gaussian()

# Running model with priors (see above)
S1_POS_PosEvnt <- brm(S1_POS_PosEvnt_form + 
                        S1_PosEvnt_miss + 
                        set_rescor(rescor = FALSE),
                      data = dat.study1_model, 
                      chains = 3,
                      prior = c(Int_prior),
                      iter = 15000,
                      warmup = 10000, 
                      control = list(adapt_delta = .99, 
                                     max_treedepth = 15))

sink(paste0(study1.model, '/S1_POS_PosEvnt.txt'))
print(summary(S1_POS_PosEvnt), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S1_POS_PosEvnt, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 100, 
           resp = "POS")+
  ggtitle("S1_POS_PosEvnt Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_POS_PosEvnt, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 12, 
           type = "error_hist", 
           resp = "POS")+
  ggtitle("S1_POS_PosEvnt Model Posterior Residuals")

# Saving Plots:
png(paste0(study1.graphics, '/S1_POS_PosEvnt_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()

save(file=paste0(study1.out, "/Posteriors/S1_POS_PosEvnt.RData"), 
     list=c("S1_POS_PosEvnt", "S1_POS_PosEvnt_form"))
remove(list=c("S1_POS_PosEvnt", "S1_POS_PosEvnt_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# Model with level 1 random effects for both negative and positive event ratings 
#   Note this model was not specifically included in the manuscript to simplify variance 
#   decomposition efforts.
S1_POS_lv1_form <- bf(
  POS | mi() ~ 1 + mi(c.NegEvnt) + mi(c.PosEvnt) + 
    (1 + mi(c.NegEvnt) + mi(c.PosEvnt)|ID)
)+gaussian()

# Running model with priors (see above)
S1_POS_lv1 <- brm(S1_POS_lv1_form +
                    S1_NegEvnt_miss + S1_PosEvnt_miss + 
                    set_rescor(rescor = FALSE),
                  data = dat.study1_model, 
                  chains = 3,
                  prior = c(Int_prior),
                  iter = 15000,
                  warmup = 10000, 
                  control = list(adapt_delta = .99, 
                                 max_treedepth = 15))

sink(paste0(study1.model, '/S1_POS_lv1.txt'))
print(summary(S1_POS_lv1), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S1_POS_lv1, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 100, 
           resp = "POS")+
  ggtitle("S1_POS_lv1 Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_POS_lv1, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 12, 
           type = "error_hist", 
           resp = "POS")+
  ggtitle("S1_POS_lv1 Model Posterior Residuals")

# Saving Plots:
png(paste0(study1.graphics, '/S1_POS_lv1_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()

save(file=paste0(study1.out, "/Posteriors/S1_POS_lv1.RData"), 
     list=c("S1_POS_lv1", "S1_POS_lv1_form"))
remove(list=c("S1_POS_lv1", "S1_POS_lv1_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# Model with level 1 random effect for recent negative event ratings and DN 
#   Note that due to the individually mean-centered level-1 effect, this model should capture an 
#   estimate of DN's total between-subjects "effect" on momentary mood. 
S1_POS_NegEvnt_DN_form <- bf(
  POS | mi() ~ 1 + mi(c.NegEvnt) + 
    c.DN + (1 + mi(c.NegEvnt)|ID)
)+gaussian()

# Running model with priors (see above)
S1_POS_NegEvnt_DN <- brm(S1_POS_NegEvnt_DN_form +
                           S1_NegEvnt_miss + 
                           set_rescor(rescor = FALSE),
                         data = dat.study1_model, 
                         chains = 3,
                         prior = c(Int_prior),
                         iter = 15000,
                         warmup = 10000, 
                         control = list(adapt_delta = .99, 
                                        max_treedepth = 15))

sink(paste0(study1.model, '/S1_POS_NegEvnt_DN.txt'))
print(summary(S1_POS_NegEvnt_DN), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S1_POS_NegEvnt_DN, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 100, 
           resp = "POS")+
  ggtitle("S1_POS_NegEvnt_DN Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_POS_NegEvnt_DN, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 12, 
           type = "error_hist", 
           resp = "POS")+
  ggtitle("S1_POS_NegEvnt_DN Model Posterior Residuals")

# Saving Plots:
png(paste0(study1.graphics, '/S1_POS_NegEvnt_DN.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()

save(file=paste0(study1.out, "/Posteriors/S1_POS_NegEvnt_DN.RData"), 
     list=c("S1_POS_NegEvnt_DN", "S1_POS_NegEvnt_DN_form"))
remove(list=c("S1_POS_NegEvnt_DN", "S1_POS_NegEvnt_DN_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# Model with level 1 random effect for recent positive event ratings and DN 
#   Note that due to the individually mean-centered level-1 effect, this model should capture an 
#   estimate of DN's total between-subjects "effect" on momentary mood. Relatedly, accounting for 
#   slight differences related to the Bayesian approach, this effect should be effectively the same
#   as the one observed in the previous model. 
S1_POS_PosEvnt_DN_form <- bf(
  POS | mi() ~ 1 + mi(c.PosEvnt) + 
    c.DN + (1 + mi(c.PosEvnt)|ID)
)+gaussian()

# Running model with priors (see above)
S1_POS_PosEvnt_DN <- brm(S1_POS_PosEvnt_DN_form +
                           S1_PosEvnt_miss + 
                           set_rescor(rescor = FALSE),
                         data = dat.study1_model, 
                         chains = 3,
                         prior = c(Int_prior),
                         iter = 15000,
                         warmup = 10000, 
                         control = list(adapt_delta = .99, 
                                        max_treedepth = 15))

sink(paste0(study1.model, '/S1_POS_PosEvnt_DN.txt'))
print(summary(S1_POS_PosEvnt_DN), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S1_POS_PosEvnt_DN, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 100, 
           resp = "POS")+
  ggtitle("S1_POS_PosEvnt_DN Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_POS_PosEvnt_DN, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 12, 
           type = "error_hist", 
           resp = "POS")+
  ggtitle("S1_POS_PosEvnt_DN Model Posterior Residuals")

# Saving plots:
png(paste0(study1.graphics, '/S1_POS_PosEvnt_DN.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()

save(file=paste0(study1.out, "/Posteriors/S1_POS_PosEvnt_DN.RData"), 
     list=c("S1_POS_PosEvnt_DN", "S1_POS_PosEvnt_DN_form"))
remove(list=c("S1_POS_PosEvnt_DN", "S1_POS_PosEvnt_DN_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# A model in which DN is included as a level-2 predictor alongside both negative and positive event 
# ratings. 
#   Note that this model was not reported in the manuscript to simplify variance extracted 
#   calculations. It is presented here for completeness and the model output is in the repo for 
#   transparency.
S1_POS_lv1_DN_form <- bf(
  POS | mi() ~ 1 + mi(c.NegEvnt) + mi(c.PosEvnt) + 
    c.DN + (1 + mi(c.NegEvnt) + mi(c.PosEvnt)|ID)
)+gaussian()

# Running model with priors (see above)
S1_POS_lv1_DN <- brm(S1_POS_lv1_DN_form + 
                       S1_NegEvnt_miss + 
                       S1_PosEvnt_miss + 
                       set_rescor(rescor = FALSE),
                     data = dat.study1_model, 
                     chains = 3,
                     prior = c(Int_prior),
                     iter = 15000,
                     warmup = 10000, 
                     control = list(adapt_delta = .99, 
                                    max_treedepth = 15))

sink(paste0(study1.model, '/S1_POS_lv1_DN.txt'))
print(summary(S1_POS_lv1_DN), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S1_POS_lv1_DN, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 100, 
           resp = "POS")+
  ggtitle("S1_POS_lv1_DN Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_POS_lv1_DN, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 12, 
           type = "error_hist", 
           resp = "POS")+
  ggtitle("S1_POS_lv1_DN Model Posterior Residuals")

# Saving plots:
png(paste0(study1.graphics, '/S1_POS_lv1_DN_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()

save(file=paste0(study1.out, "/Posteriors/S1_POS_lv1_DN.RData"), 
     list=c("S1_POS_lv1_DN", "S1_POS_lv1_DN_form"))
remove(list=c("S1_POS_lv1_DN", "S1_POS_lv1_DN_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# This model includes *average* negative event ratings as a measure of an individual's overall 
# negative event context. This is our proxy for overall exposure to more intense negative events. 
S1_POS_NegEvnt_Exp_form <- bf(
  POS | mi() ~ 1 + mi(c.NegEvnt) + 
    m.NegEvnt + (1 + mi(c.NegEvnt)|ID)
)+gaussian()

# Running model with priors (see above)
S1_POS_NegEvnt_Exp <- brm(S1_POS_NegEvnt_Exp_form +
                            S1_NegEvnt_miss + 
                            set_rescor(rescor = FALSE),
                          data = dat.study1_model, 
                          chains = 3,
                          prior = c(Int_prior),
                          iter = 15000,
                          warmup = 10000,  
                          control = list(adapt_delta = .99, 
                                         max_treedepth = 15))

sink(paste0(study1.model, '/S1_POS_NegEvnt_Exp.txt'))
print(summary(S1_POS_NegEvnt_Exp), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S1_POS_NegEvnt_Exp, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 100, 
           resp = "POS")+
  ggtitle("S1_POS_NegEvnt_Exp Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_POS_NegEvnt_Exp, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 12, 
           type = "error_hist", 
           resp = "POS")+
  ggtitle("S1_POS_NegEvnt_Exp Model Posterior Residuals")

# Saving plots:
png(paste0(study1.graphics, '/S1_POS_NegEvnt_Exp.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()

save(file=paste0(study1.out, "/Posteriors/S1_POS_NegEvnt_Exp.RData"), 
     list=c("S1_POS_NegEvnt_Exp", "S1_POS_NegEvnt_Exp_form"))
remove(list=c("S1_POS_NegEvnt_Exp", "S1_POS_NegEvnt_Exp_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# Model is the same as the one above, with the exception that the focus is on inclusion of between-
# subjects differences in average positive mood rating. The idea is similar though in that this is 
# our proxy for the overall contextual effect of reporting more intense positive events on average.
S1_POS_PosEvnt_Exp_form <- bf(
  POS | mi() ~ 1 + mi(c.PosEvnt) + 
    m.PosEvnt + (1 + mi(c.PosEvnt)|ID)
)+gaussian()

# Running model with priors (see above)
S1_POS_PosEvnt_Exp <- brm(S1_POS_PosEvnt_Exp_form +
                            S1_PosEvnt_miss + 
                            set_rescor(rescor = FALSE),
                          data = dat.study1_model, 
                          chains = 3,
                          prior = c(Int_prior),
                          iter = 15000,
                          warmup = 10000, 
                          control = list(adapt_delta = .99, 
                                         max_treedepth = 15))

sink(paste0(study1.model, '/S1_POS_PosEvnt_Exp.txt'))
print(summary(S1_POS_PosEvnt_Exp), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S1_POS_PosEvnt_Exp, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 100, 
           resp = "POS")+
  ggtitle("S1_POS_PosEvnt_Exp Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_POS_PosEvnt_Exp, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 12, 
           type = "error_hist", 
           resp = "POS")+
  ggtitle("S1_POS_PosEvnt_Exp Model Posterior Residuals")

# Saving plots:
png(paste0(study1.graphics, '/S1_POS_PosEvnt_Exp.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()

save(file=paste0(study1.out, "/Posteriors/S1_POS_PosEvnt_Exp.RData"), 
     list=c("S1_POS_PosEvnt_Exp", "S1_POS_PosEvnt_Exp_form"))
remove(list=c("S1_POS_PosEvnt_Exp", "S1_POS_PosEvnt_Exp_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# This model includes both level-1 random effects and individual average ratings for both negative
# and positive event ratings as predictors. 
#   Note that as with model models that included both event types as predictors, these models were
#   not included in the manuscript - but are included here along with the relevant outputs elsewhere
#   in this repo. 
S1_POS_lv1_Exp_form <- bf(
  POS | mi() ~ 1 + mi(c.NegEvnt) + mi(c.PosEvnt) + 
    m.NegEvnt + m.PosEvnt + (1 + mi(c.NegEvnt) + mi(c.PosEvnt)|ID)
)+gaussian()

# Running model with priors (see above)
S1_POS_lv1_Exp <- brm(S1_POS_lv1_Exp_form +
                        S1_NegEvnt_miss +
                        S1_PosEvnt_miss + 
                        set_rescor(rescor = FALSE),
                      data = dat.study1_model, 
                      chains = 3,
                      prior = c(Int_prior),
                      iter = 15000,
                      warmup = 10000, 
                      control = list(adapt_delta = .99, 
                                     max_treedepth = 15))

sink(paste0(study1.model, '/S1_POS_lv1_Exp.txt'))
print(summary(S1_POS_lv1_Exp), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S1_POS_lv1_Exp, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 100, 
           resp = "POS")+
  ggtitle("S1_POS_lv1_Exp Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_POS_lv1_Exp, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 12, 
           type = "error_hist", 
           resp = "POS")+
  ggtitle("S1_POS_lv1_Exp Model Posterior Residuals")

# Saving Plots:
png(paste0(study1.graphics, '/S1_POS_lv1_Exp_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()

save(file=paste0(study1.out, "/Posteriors/S1_POS_lv1_Exp.RData"), 
     list=c("S1_POS_lv1_Exp", "S1_POS_lv1_Exp_form"))
remove(list=c("S1_POS_lv1_Exp", "S1_POS_lv1_Exp_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# This model in conjunction with the models including just DN or just average negative event ratings
# will be used to isolate the amount of "unique" variance attributable to DN, and to overall 
# negative contexts
S1_POS_NegEvnt_DN_Exp_form <- bf(
  POS | mi() ~ 1 + mi(c.NegEvnt) + 
    c.DN + m.NegEvnt + (1 + mi(c.NegEvnt)|ID)
)+gaussian()

# Running model with priors (see above)
S1_POS_NegEvnt_DN_Exp <- brm(S1_POS_NegEvnt_DN_Exp_form +
                               S1_NegEvnt_miss + 
                               set_rescor(rescor = FALSE),
                             data = dat.study1_model, 
                             chains = 3,
                             prior = c(Int_prior),
                             iter = 15000,
                             warmup = 10000,  
                             control = list(adapt_delta = .99, 
                                            max_treedepth = 15))

sink(paste0(study1.model, '/S1_POS_NegEvnt_DN_Exp.txt'))
print(summary(S1_POS_NegEvnt_DN_Exp), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S1_POS_NegEvnt_DN_Exp, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 100, 
           resp = "POS")+
  ggtitle("S1_POS_NegEvnt_DN_Exp Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_POS_NegEvnt_DN_Exp, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 12, 
           type = "error_hist", 
           resp = "POS")+
  ggtitle("S1_POS_NegEvnt_DN_Exp Model Posterior Residuals")

# Saving plots:
png(paste0(study1.graphics, '/S1_POS_NegEvnt_DN_Exp.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()

save(file=paste0(study1.out, "/Posteriors/S1_POS_NegEvnt_DN_Exp.RData"), 
     list=c("S1_POS_NegEvnt_DN_Exp", "S1_POS_NegEvnt_DN_Exp_form"))
remove(list=c("S1_POS_NegEvnt_DN_Exp", "S1_POS_NegEvnt_DN_Exp_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# This model in conjunction with the models including just DN or just average positive event ratings
# will be used to isolate the amount of "unique" variance attributable to DN, and to overall 
# positive contexts
S1_POS_PosEvnt_DN_Exp_form <- bf(
  POS | mi() ~ 1 + mi(c.PosEvnt) + 
    c.DN + m.PosEvnt + (1 + mi(c.PosEvnt)|ID)
)+gaussian()

# Running model with priors (see above)
S1_POS_PosEvnt_DN_Exp <- brm(S1_POS_PosEvnt_DN_Exp_form +
                               S1_PosEvnt_miss + 
                               set_rescor(rescor = FALSE),
                             data = dat.study1_model, 
                             chains = 3,
                             prior = c(Int_prior),
                             iter = 15000,
                             warmup = 10000, 
                             control = list(adapt_delta = .99, 
                                            max_treedepth = 15))

sink(paste0(study1.model, '/S1_POS_PosEvnt_DN_Exp.txt'))
print(summary(S1_POS_PosEvnt_DN_Exp), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S1_POS_PosEvnt_DN_Exp, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 100, 
           resp = "POS")+
  ggtitle("S1_POS_PosEvnt_DN_Exp Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_POS_PosEvnt_DN_Exp, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 12, 
           type = "error_hist", 
           resp = "POS")+
  ggtitle("S1_POS_PosEvnt_DN_Exp Model Posterior Residuals")

# Saving plots:
png(paste0(study1.graphics, '/S1_POS_PosEvnt_DN_Exp.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()

save(file=paste0(study1.out, "/Posteriors/S1_POS_PosEvnt_DN_Exp.RData"), 
     list=c("S1_POS_PosEvnt_DN_Exp", "S1_POS_PosEvnt_DN_Exp_form"))
remove(list=c("S1_POS_PosEvnt_DN_Exp", "S1_POS_PosEvnt_DN_Exp_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Agan for completeness we have included a model with all three variables and their "main" effects 
# expressed in a multilevel model. 
#   Note that as with model models that included both event types as predictors, these models were
#   not included in the manuscript - but are included here along with the relevant outputs elsewhere
#   in this repo. 
S1_POS_lv1_DN_Exp_form <- bf(
  POS | mi() ~ 1 + mi(c.NegEvnt) + mi(c.PosEvnt) + 
    c.DN + m.NegEvnt + m.PosEvnt + (1 + mi(c.NegEvnt) + mi(c.PosEvnt)|ID)
)+gaussian()

# Running model with priors (see above)
S1_POS_lv1_DN_Exp <- brm(S1_POS_lv1_DN_Exp_form +
                           S1_NegEvnt_miss +
                           S1_PosEvnt_miss + 
                           set_rescor(rescor = FALSE),
                         data = dat.study1_model, 
                         chains = 3,
                         prior = c(Int_prior),
                         iter = 15000,
                         warmup = 10000, 
                         control = list(adapt_delta = .99, 
                                        max_treedepth = 15))

sink(paste0(study1.model, '/S1_POS_lv1_DN_Exp.txt'))
print(summary(S1_POS_lv1_DN_Exp), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S1_POS_lv1_DN_Exp, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 100, 
           resp = "POS")+
  ggtitle("S1_POS_lv1_DN_Exp Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_POS_lv1_DN_Exp, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 12, 
           type = "error_hist", 
           resp = "POS")+
  ggtitle("S1_POS_lv1_DN_Exp Model Posterior Residuals")

# Saving plots:
png(paste0(study1.graphics, '/S1_POS_lv1_DN_Exp_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()

save(file=paste0(study1.out, "/Posteriors/S1_POS_lv1_DN_Exp.RData"), 
     list=c("S1_POS_lv1_DN_Exp", "S1_POS_lv1_DN_Exp_form"))
remove(list=c("S1_POS_lv1_DN_Exp", "S1_POS_lv1_DN_Exp_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# This model actually represents the "final" model for this modeling tree (DN, negative events, 
# positive mood). All variance components obtained in the present set of analyses stem from 
# isolated R2 differences moving from the unconditional model to this model. 
S1_POS_NegEvnt_Rct_form <- bf(
  POS | mi() ~ 1 + mi(c.NegEvnt)*c.DN + 
    m.NegEvnt + (1 + mi(c.NegEvnt)|ID)
)+gaussian()

# Running model with priors (see above)
S1_POS_NegEvnt_Rct <- brm(S1_POS_NegEvnt_Rct_form +
                            S1_NegEvnt_miss + 
                            set_rescor(rescor = FALSE),
                          data = dat.study1_model, 
                          chains = 3,
                          prior = c(Int_prior),
                          iter = 15000,
                          warmup = 10000, 
                          control = list(adapt_delta = .99, 
                                         max_treedepth = 15))

sink(paste0(study1.model, '/S1_POS_NegEvnt_Rct.txt'))
print(summary(S1_POS_NegEvnt_Rct), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S1_POS_NegEvnt_Rct, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 100, 
           resp = "POS")+
  ggtitle("S1_POS_NegEvnt_Rct Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_POS_NegEvnt_Rct, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 12, 
           type = "error_hist", 
           resp = "POS")+
  ggtitle("S1_POS_NegEvnt_Rct Model Posterior Residuals")

# Saving plots:
png(paste0(study1.graphics, '/S1_POS_NegEvnt_Rct.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()

save(file=paste0(study1.out, "/Posteriors/S1_POS_NegEvnt_Rct.RData"), 
     list=c("S1_POS_NegEvnt_Rct", "S1_POS_NegEvnt_Rct_form"))
remove(list=c("S1_POS_NegEvnt_Rct", "S1_POS_NegEvnt_Rct_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# This model actually represents the "final" model for this modeling tree (DN, positive events, 
# positive mood). All variance components obtained in the present set of analyses stem from 
# isolated R2 differences moving from the unconditional model to this model. 
S1_POS_PosEvnt_Rct_form <- bf(
  POS | mi() ~ 1 + mi(c.PosEvnt)*c.DN + 
    m.PosEvnt + (1 + mi(c.PosEvnt)|ID)
)+gaussian()

# Running model with priors (see above)
S1_POS_PosEvnt_Rct <- brm(S1_POS_PosEvnt_Rct_form +
                            S1_PosEvnt_miss + 
                            set_rescor(rescor = FALSE),
                          data = dat.study1_model, 
                          chains = 3,
                          prior = c(Int_prior),
                          iter = 15000,
                          warmup = 10000, 
                          control = list(adapt_delta = .99, 
                                         max_treedepth = 15))

sink(paste0(study1.model, '/S1_POS_PosEvnt_Rct.txt'))
print(summary(S1_POS_PosEvnt_Rct), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S1_POS_PosEvnt_Rct, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 100, 
           resp = "POS")+
  ggtitle("S1_POS_PosEvnt_Rct Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_POS_PosEvnt_Rct, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 12, 
           type = "error_hist", 
           resp = "POS")+
  ggtitle("S1_POS_PosEvnt_Rct Model Posterior Residuals")

#Saving plots:
png(paste0(study1.graphics, '/S1_POS_PosEvnt_Rct.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()

save(file=paste0(study1.out, "/Posteriors/S1_POS_PosEvnt_Rct.RData"), 
     list=c("S1_POS_PosEvnt_Rct", "S1_POS_PosEvnt_Rct_form"))
remove(list=c("S1_POS_PosEvnt_Rct", "S1_POS_PosEvnt_Rct_form"))
gc()
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# Simultaneous inclustion of all DN reactivity effects for completness 
#   Note not presented in the manuscript (like similar models defined throughout)
S1_POS_lv1_Rct_form <- bf(
  POS | mi() ~ 1 + mi(c.NegEvnt)*c.DN + mi(c.PosEvnt)*c.DN + 
    m.NegEvnt + m.PosEvnt + (1 + mi(c.NegEvnt) + mi(c.PosEvnt)|ID)
)+gaussian()

# Running model with priors (see above)
S1_POS_lv1_Rct <- brm(S1_POS_lv1_Rct_form +
                        S1_NegEvnt_miss +
                        S1_PosEvnt_miss + 
                        set_rescor(rescor = FALSE),
                      data = dat.study1_model, 
                      chains = 3,
                      prior = c(Int_prior),
                      iter = 15000,
                      warmup = 10000, 
                      control = list(adapt_delta = .99, 
                                     max_treedepth = 15))

sink(paste0(study1.model, '/S1_POS_lv1_Rct.txt'))
print(summary(S1_POS_lv1_Rct), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S1_POS_lv1_Rct, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 100, 
           resp = "POS")+
  ggtitle("S1_POS_lv1_Rct Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_POS_lv1_Rct, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 12, 
           type = "error_hist", 
           resp = "POS")+
  ggtitle("S1_POS_lv1_Rct Model Posterior Residuals")

# Saving Plots:
png(paste0(study1.graphics, '/S1_POS_lv1_Rct_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()

save(file=paste0(study1.out, "/Posteriors/S1_POS_lv1_Rct.RData"), 
     list=c("S1_POS_lv1_Rct", "S1_POS_lv1_Rct_form"))
remove(list=c("S1_POS_lv1_Rct", "S1_POS_lv1_Rct_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# The "Flr" part of this model is to address possible floor/ceiling effects could be 
# partially to blame for any DN reactivity effects detected in the models. For instance, one could 
# possibly argue that average ratings of recent negative events could "limit" the range of response
# and only those with less (or more) intense negative events had much room to move in terms of 
# momentary mood. 
S1_POS_NegEvnt_Flr_form <- bf(
  POS | mi() ~ 1 + mi(c.NegEvnt)*c.DN + mi(c.NegEvnt)*m.NegEvnt +
    m.NegEvnt + (1 + mi(c.NegEvnt)|ID)
)+gaussian()

# Running model with priors (see above)
S1_POS_NegEvnt_Flr <- brm(S1_POS_NegEvnt_Flr_form +
                            S1_NegEvnt_miss + 
                            set_rescor(rescor = FALSE),
                          data = dat.study1_model, 
                          chains = 3,
                          prior = c(Int_prior),
                          iter = 15000,
                          warmup = 10000, 
                          control = list(adapt_delta = .99, 
                                         max_treedepth = 15))

sink(paste0(study1.model, '/S1_POS_NegEvnt_Flr.txt'))
print(summary(S1_POS_NegEvnt_Flr), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S1_POS_NegEvnt_Flr, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 100, 
           resp = "POS")+
  ggtitle("S1_POS_NegEvnt_Flr Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_POS_NegEvnt_Flr, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 12, 
           type = "error_hist", 
           resp = "POS")+
  ggtitle("S1_POS_NegEvnt_Flr Model Posterior Residuals")

# Saving Plots:
png(paste0(study1.graphics, '/S1_POS_NegEvnt_Flr.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()

save(file=paste0(study1.out, "/Posteriors/S1_POS_NegEvnt_Flr.RData"), 
     list=c("S1_POS_NegEvnt_Flr", "S1_POS_NegEvnt_Flr_form"))
remove(list=c("S1_POS_NegEvnt_Flr", "S1_POS_NegEvnt_Flr_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# See note above for previous model. This model is essentially the positive event version. 
S1_POS_PosEvnt_Flr_form <- bf(
  POS | mi() ~ 1 + mi(c.PosEvnt)*c.DN + mi(c.PosEvnt)*m.PosEvnt +
    m.PosEvnt + (1 + mi(c.PosEvnt)|ID)
)+gaussian()

# Running model with priors (see above)
S1_POS_PosEvnt_Flr <- brm(S1_POS_PosEvnt_Flr_form +
                            S1_PosEvnt_miss + 
                            set_rescor(rescor = FALSE),
                          data = dat.study1_model, 
                          chains = 3,
                          prior = c(Int_prior),
                          iter = 15000,
                          warmup = 10000, 
                          control = list(adapt_delta = .99, 
                                         max_treedepth = 15))

sink(paste0(study1.model, '/S1_POS_PosEvnt_Flr.txt'))
print(summary(S1_POS_PosEvnt_Flr), digits = 5)
sink()

# Simple model check plotting:
ppc_density <- 
  pp_check(S1_POS_PosEvnt_Flr, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 100, 
           resp = "POS")+
  ggtitle("S1_POS_PosEvnt_Flr Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_POS_PosEvnt_Flr, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 12, 
           type = "error_hist", 
           resp = "POS")+
  ggtitle("S1_POS_PosEvnt_Flr Model Posterior Residuals")

# Saving Plots:
png(paste0(study1.graphics, '/S1_POS_PosEvnt_Flr.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()

save(file=paste0(study1.out, "/Posteriors/S1_POS_PosEvnt_Flr.RData"), 
     list=c("S1_POS_PosEvnt_Flr", "S1_POS_PosEvnt_Flr_form"))
remove(list=c("S1_POS_PosEvnt_Flr", "S1_POS_PosEvnt_Flr_form"))
gc()

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# Really just a complete model with all possible effects of interest, both those specifically 
# highlighted in the research questions posed in the manuscript and tangentially related effects. 
S1_POS_lv1_Flr_form <- bf(
  POS | mi() ~ 1 + mi(c.NegEvnt)*c.DN + mi(c.PosEvnt)*c.DN + 
    mi(c.NegEvnt)*m.NegEvnt + mi(c.PosEvnt)*m.PosEvnt +
    m.NegEvnt + m.PosEvnt + (1 + mi(c.NegEvnt) + mi(c.PosEvnt)|ID)
)+gaussian()

#Running model with priors (see above)
S1_POS_lv1_Flr <- brm(S1_POS_lv1_Flr_form +
                        S1_NegEvnt_miss +
                        S1_PosEvnt_miss + 
                        set_rescor(rescor = FALSE),
                      data = dat.study1_model, 
                      chains = 3,
                      prior = c(Int_prior),
                      iter = 15000,
                      warmup = 10000, 
                      control = list(adapt_delta = .99, 
                                     max_treedepth = 15))

sink(paste0(study1.model, '/S1_POS_lv1_Flr.txt'))
print(summary(S1_POS_lv1_Flr), digits = 5)
sink()

#Simple Model Check plotting:
ppc_density <- 
  pp_check(S1_POS_lv1_Flr, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 100, 
           resp = "POS")+
  ggtitle("S1_POS_lv1_Flr Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_POS_lv1_Flr, 
           newdata = dat.study1_model[!is.na(dat.study1_model$POS),],
           nsamples = 12, 
           type = "error_hist", 
           resp = "POS")+
  ggtitle("S1_POS_lv1_Flr Model Posterior Residuals")

#Saving Plots:
png(paste0(study1.graphics, '/S1_POS_lv1_Flr_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()

save(file=paste0(study1.out, "/Posteriors/S1_POS_lv1_Flr.RData"), 
     list=c("S1_POS_lv1_Flr", "S1_POS_lv1_Flr_form"))
remove(list=c("S1_POS_lv1_Flr", "S1_POS_lv1_Flr_form"))
gc()
