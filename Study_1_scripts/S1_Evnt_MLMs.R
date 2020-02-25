###################################################################################################
# Study 1 Modeling Script

# Negative and Positive event MLMs
# Description: 
#   The analyses below examine a bivariate association between individual dispositional negativity 
#   scores and ratings of recent negative and positive events. Participants provided event ratings 
#   via repeated ecological momentary assessments completed on their smartphones. Dispositional 
#   negativity scores were created by aggregating the same set of items tapping trait anxiety and 
#   negative affect completed at two separate timepoints and in two different contexts. 

# Modeling Notes: 
#   * Exploratory analyses revealed that negative event ratings were positively skewed, thus a 
#   weakly informative lognormal prior was chosen for the negative event model
#   * Missingness was addressed at runtime by taking draws from the posterior predictive 
#   distribution
###################################################################################################

#--------------------------------------------------------------------------------------------------
# Loading relevant packages
library(brms)
library(rstan)
library(rstanarm)
library(bayesplot)
library(tidyverse)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
user<-ifelse(Sys.getenv("USERPROFILE")=="", "~", Sys.getenv("USERPROFILE"))
wd<-paste0(user, '/Dropbox/UMD/Shackman Lab/EMA_MS')
data.folder<-paste0(wd, '/Data')
study1.out<-paste0(wd, '/Study 1 output')
study1.graphics<-paste0(study1.out, '/Graphics')
study1.model<-paste0(study1.out, '/Model summaries')
stan.code<-paste0(wd, '/Stan_code')
EDA.folder<-paste0(study1.out, '/EDA')
#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
# Loading Study 1 Data (from Emotion MS - Shackman et al. 2017)
load(paste0(data.folder, '/Emotion MS environment.RData'))

dat.study1_lv1 <- dat %>% 
  transmute(ID = subid,
            NEG = MAFS_NA, 
            POS = MAFS_PA,
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

dat.study1_lv2 <- dat %>% 
  group_by(subid) %>% 
  summarize(ID = subid[1],
            c.DN = mean(ZAP_Both),
            BFI_A = mean(Battery_BFI_A), 
            BFI_O = mean(Battery_BFI_O), 
            BFI_C = mean(Battery_BFI_C), 
            BFI_N = mean(Battery_BFI_N), 
            m.NegEvnt = mean(WorstEvent_Neg, na.rm = TRUE), 
            m.PosEvnt = mean(BestEvent_Pos, na.rm = TRUE)) 

dat.study1_model <- merge(dat.study1_lv1, dat.study1_lv2, by = "ID")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
# Simple Event Models:

# Substantive Questions: 
#   1. Is dispositional negativity associated with intensity of recent negative events? 
#   2. Is dispositional negativity associated with intensity of recent positive events? 
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
#Negative Event Model: 
NegEvnt_DN <- bf(
  NegEvnt | mi() ~ 1 + c.DN + (1|ID) 
)+lognormal()

#Creating priors for intercept to help constrain final model in response space
mu <- mean(log(dat.study1_model$NegEvnt), na.rm=TRUE)
sigma <- sd(log(dat.study1_model$NegEvnt), na.rm = TRUE)

Int_prior <- set_prior(paste0("normal(", mu, ",", sigma, ")"), 
                       class = "Intercept")

#Running model with priors (see above)
S1_NegEvnt_DN <- brm(NegEvnt_DN,
                     data = dat.study1_model, 
                     chains = 3,
                     prior = c(Int_prior),
                     iter = 15000,
                     warmup = 10000, 
                     control = list(adapt_delta = .99, 
                                    max_treedepth = 15))

save(list = c("S1_NegEvnt_DN", "dat.study1_model"), 
     file = paste0(data.folder, 
                   '/Model_Posteriors', 
                   '/S1_NegEvnt_DN.RData'))

sink(paste0(study1.model, '/S1_NegEvnt_DN.txt'))
print(summary(S1_NegEvnt_DN), digits = 5)
sink()

#Simple Model Check plotting:
ppc_density <- 
pp_check(S1_NegEvnt_DN, 
         newdata = dat.study1_model[!is.na(dat.study1_model$NegEvnt),], 
         nsamples = 100)+
  ggtitle("S1_NegEvnt_DN Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_NegEvnt_DN, 
           newdata = dat.study1_model[!is.na(dat.study1_model$NegEvnt),], 
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S1_NegEvnt_DN Model Posterior Residuals")

png(paste0(study1.graphics, '/S1_NegEvnt_DN_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()
remove("S1_NegEvnt_DN")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-
#Positive Event Model: 
PosEvnt_DN <- bf(
  PosEvnt | mi() ~ 1 + c.DN + (1|ID) 
)

#Creating priors for intercept to help constrain final model in response space
mu <- mean(dat.study1_model$PosEvnt, na.rm=TRUE)
sigma <- sd(dat.study1_model$PosEvnt, na.rm = TRUE)

Int_prior <- set_prior(paste0("normal(", mu, ",", sigma, ")"), 
                       class = "Intercept")

#Running model with priors (see above)
S1_PosEvnt_DN <- brm(PosEvnt_DN,
                     data = dat.study1_model, 
                     chains = 3,
                     prior = c(Int_prior),
                     iter = 15000,
                     warmup = 10000, 
                     control = list(adapt_delta = .99, 
                                    max_treedepth = 15))

save(list = c("S1_PosEvnt_DN", "dat.study1_model"), 
     file = paste0(data.folder, 
                   '/Model_Posteriors', 
                   '/S1_PosEvnt_DN.RData'))

sink(paste0(study1.model, '/S1_PosEvnt_DN.txt'))
print(summary(S1_PosEvnt_DN), digits = 5)
sink()

#Simple Model Check plotting:
ppc_density <- 
  pp_check(S1_PosEvnt_DN, 
           newdata = dat.study1_model[!is.na(dat.study1_model$PosEvnt),], 
           nsamples = 100)+
  ggtitle("S1_PosEvnt_DN Posterior Predictive Distribution")

ppc_hist <- 
  pp_check(S1_PosEvnt_DN, 
           newdata = dat.study1_model[!is.na(dat.study1_model$PosEvnt),], 
           nsamples = 12, 
           type = "error_hist")+
  ggtitle("S1_PosEvnt_DN Model Posterior Residuals")

png(paste0(study1.graphics, '/S1_PosEvnt_DN_ppc.png'), 
    units = "in", 
    height = 5.5, 
    width = 11, 
    res = 900)
cowplot::plot_grid(ppc_hist, 
                   ppc_density)

dev.off()
remove("S1_PosEvnt_DN")
