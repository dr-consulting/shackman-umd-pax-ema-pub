############################################################################################################
#Authors:     Barstead, M.G., DeYoung, K.D., Anderson, A. S., & Shackman, A. J.

#Title:       The moment-to-moment affective experience of dispositionally negative individuals 

#Contact:     barstead@umd.edu

#Contents:    Bayesian regression analyses - Study 2
############################################################################################################
############################################################################################################
library(brms)
library(rstan)
library(rstanarm)
library(bayesplot)
library(pan)
library(mitml)
library(mice)
library(parallel)
library(RColorBrewer)
library(ggridges)
library(riverplot)
library(ggpubr)
library(tidyverse)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
############################################################################################################
user<-ifelse(Sys.getenv("USERPROFILE")=="", "~", Sys.getenv("USERPROFILE"))
wd<-paste0(user, '/Dropbox/UMD/Shackman Lab/EMA_MS')
data.folder<-paste0(wd, '/Data')
study2.out<-paste0(wd, '/Study 2 output')
study2.graphics<-paste0(study2.out, '/Graphics')
study2.model<-paste0(study2.out, '/Model summaries')
stan.code<-paste0(wd, '/Stan_code')
EDA.folder<-paste0(study2.out, '/EDA')

#Loading Study 1 Data (from Emotion MS - Shackman et al. 2017)
load(paste0(data.folder, '/Study2_Model_Data_All_models.RData'))
user<-ifelse(Sys.getenv("USERPROFILE")=="", "~", Sys.getenv("USERPROFILE"))
wd<-paste0(user, '/Dropbox/UMD/Shackman Lab/EMA_MS')
data.folder<-paste0(wd, '/Data')
study2.out<-paste0(wd, '/Study 2 output')
study2.graphics<-paste0(study2.out, '/Graphics')
study2.model<-paste0(study2.out, '/Model summaries')
stan.code<-paste0(wd, '/Stan_code')
EDA.folder<-paste0(study2.out, '/EDA')

############################################################################################################
#Creating a Function to Finalize Imputed Model Output: 
bayes.to.txt<-function(model=NULL, 
                       out.folder=NULL, 
                       file=NULL, 
                       DF=NULL, 
                       tot.pars=3, 
                       impute=TRUE){
  sink(paste0(out.folder, '/', file, '.txt'))
  cat('System Information:')
  cat('\n=====================================================================================')
  cat(paste0('\nProcessor:', '\t\t\t', benchmarkme::get_cpu()$model_name))
  cat(paste0('\nNumber of Threads:', '\t\t', parallel::detectCores(logical=T)))
  cat(paste0('\nRAM:', '\t\t\t\t', paste(round(benchmarkme::get_ram()/1073741824), 'GB')))
  cat('\n=====================================================================================\n')
  
  cat('\n\nModel Information:')
  cat('\n=====================================================================================')
  cat('\n\nFormula (lme4 syntax):')
  cat('\n-------------------------------------------------------------------------------------\n')
  print(model$formula)
  cat('\n\nPriors:')
  cat('\n-------------------------------------------------------------------------------------\n')
  print(prior_summary(model))
  cat('\n\nStan Code:')
  cat('\n-------------------------------------------------------------------------------------\n')
  print(make_stancode(model$formula, data=DF, family=model$family$family))
  cat('\n=====================================================================================\n')
  
  cat('\n\nStan Arguments:')
  cat('\n=====================================================================================')
  cat(paste0('\nAdapt Delta:', '\t\t\t', model$fit@stan_args[[1]]$control$adapt_delta))
  cat(paste0('\nMaximum Tree Depth:', '\t\t', model$fit@stan_args[[1]]$control$max_treedepth))
  cat('\n=====================================================================================\n')
  
  cat('\n\nVariance Explained:')
  cat('\n=====================================================================================')
  cat('\n\nResidual Variance:')
  cat('\n-------------------------------------------------------------------------------------\n')
  print(sjstats::icc(model, posterior = T), prob =.95, digits=5)
  cat('\n\nBayesian R-squared (overall variance explained):')
  cat('\n-------------------------------------------------------------------------------------\n')
  print(bayes_R2(model), digits=5)
  cat('\n=====================================================================================\n')
  
  cat('\n\nModel Summary:')
  cat('\n=====================================================================================\n')
  print(summary(model), digits=5)
  
  if(impute==0){
    cat('WARNING! If using multiply imputed datasets ignore R=hat')
    cat('\nUse the code below to obtain interpretable R=hat values for each data set')
    cat('\nround(modelname$rhats[,1:tot.pars], 3)')  
    cat('\nwhere "tot.pars" is the total number of model parameters to return (usually fixed effects)')
  }
  else{
    cat('\n\nPotential Scale Reduction Factor for Each Imputed Dataset:')
    cat('\n=====================================================================================\n')
    print(round(model$rhats[,1:tot.pars], 3))
    cat('\n=====================================================================================\n')
  }
  sink()
}


#Quick inspection and cleanup
psych::describe(dat.study2)
dat.study2$DN_comb<-as.numeric(scale(dat.study2$DN_comb))

#Multilevel imputation will take an incredibly long time (even when running on eleven cores)

#inspecting missing patterns in the data: 
md.pattern(dat.study2)

#Approximately 20% of the data is missing - will impute using a multilevel approach

#Exploring differences in missing vs. non-missing distributions of scores: 
#Missinginess is almost entirely clustered - will explore using single variable
dat.study2$miss<-ifelse(!is.na(dat.study2$POS), 0, 1)

names(dat.study2)
fit.miss<-lme4::glmer(miss~1+DN_comb+Male+SIG+
                        avg.PosEvnt+avg.NegEvnt+
                        (1+SIG|ID), 
                      data = dat.study2, 
                      control = glmerControl(optimizer = "bobyqa", 
                                             tol = .00001), 
                      family = 'binomial')
summary(fit.miss) #Just a significant linear effect of time (order of survey - will include to reduce noise in imputation models)

#Going to include other gender and centered DN as predictors as well 
dat.study2$NegEvnt_char<-as.factor(dat.study2$NegEvnt)
dat.study2$PosEvnt_char<-as.factor(dat.study2$PosEvnt)
fml<- JOY + CALM + ANX + ANG + TIRED + DEP + PosEvnt_char + NegEvnt_char ~ 
  1 + SIG + DN_comb + Male + (1 + SIG|ID)

#Set number of data sets to impute 
M<-20
imp<-jomoImpute(dat.study2, 
               formula=fml, 
               n.burn=100000, 
               n.iter = 5000, 
               m=M, 
               seed = 021219)

#plot(imp)

#plotting results - assessing convergence on final distributions
dat.imp<-mitmlComplete(imp)

#Inspecting imputation quality properties - especially interested in outcome measures  
dat.long<-data.frame()
for(i in 1:M){
  dat.temp<-dat.imp[[i]]
  dat.temp$IMP<-rep(i, length(dat.temp[,1]))
  dat.long<-rbind(dat.long, dat.temp)
}

dat.study2$IMP<-rep(0, length(dat.study2[,1]))
dat.long<-rbind(dat.study2, dat.long)
dat.long$Orig<-ifelse(dat.long$IMP<1, 'Original', 'Imputed')

#Plotting Negative Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=NEG, group=Orig, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Negative Affect')+
  ylab('Density')+
  xlab('Momentary Negative Affect')
g1

png(paste0(study2.graphics, '/S2_Imputation_NegativeAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Positive Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=POS, group=Orig, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Positive Affect')+
  ylab('Density')+
  xlab('Momentary Positive Affect')
g1

png(paste0(study2.graphics, '/S2_Imputation_PositiveAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#summarizing imputations - checking for convergence
sink(paste0(study2.out, '/S2_Imputation_Summary.txt'))
summary(imp)
psych::describeBy(dat.long, group='IMP')
sink()    #warnings have to do with factor variables for the most part

#Quirky formatting problems - getting around it this way... 
#   The issue is that brm_multiple requires a list of datasets
#   It cannot recognize an mitml.list object though
#   So I made a list out of list

dat.study2.list<-list()
for(i in 1:M){
  dat.study2.list[[i]]<-dat.imp[[i]]
}

#Obtaining individual mean ratings for Best and Worst events
#Individually mean centering predictors to decompose between and within effects

#Checking to make sure events imputed correctly
#table(dat.study2.list[[1]]$NegEvnt_char)
#table(dat.study2.list[[2]]$NegEvnt_char)
#table(dat.study2.list[[3]]$NegEvnt_char)
#table(dat.study2.list[[4]]$NegEvnt_char)
#table(dat.study2.list[[5]]$NegEvnt_char)

#Checking to makes sure positive events imputed correctly
#table(dat.study2.list[[1]]$NegEvnt_char)
#table(dat.study2.list[[2]]$NegEvnt_char)
#table(dat.study2.list[[3]]$NegEvnt_char)
#table(dat.study2.list[[4]]$NegEvnt_char)
#table(dat.study2.list[[5]]$NegEvnt_char)

for(i in 1:M){
  IDs<-unique(dat.study2$ID)
  DF.temp<-dat.study2.list[[i]]
  for(j in 1:length(IDs)){
    DF.temp$prop_NegEvnt[DF.temp$ID == IDs[j]]<-
      mean(as.numeric(as.character(DF.temp$NegEvnt_char[DF.temp$ID == IDs[j]])), 
           na.rm = TRUE)
    
    DF.temp$prop_PosEvnt[DF.temp$ID == IDs[j]]<-
      mean(as.numeric(as.character(DF.temp$PosEvnt_char[DF.temp$ID == IDs[j]])), 
           na.rm = TRUE)
    
  }
  DF.temp$T.NEG<-(DF.temp$DEP*2+DF.temp$ANX*3)/5
  DF.temp$T.NEG[DF.temp$T.NEG<1]<-1
  DF.temp$T.NEG[DF.temp$T.NEG>5]<-5
  
  DF.temp$T.POS<-(DF.temp$JOY + DF.temp$CALM)/2
  DF.temp$T.POS[DF.temp$T.POS<1]<-1
  DF.temp$T.POS[DF.temp$T.POS>5]<-5
  
  DF.temp$T.JOY<-DF.temp$JOY
  DF.temp$T.JOY[DF.temp$T.JOY<1]<-1
  DF.temp$T.JOY[DF.temp$T.JOY>5]<-5
  
  DF.temp$T.CALM<-DF.temp$CALM
  DF.temp$T.CALM[DF.temp$T.CALM<1]<-1
  DF.temp$T.CALM[DF.temp$T.CALM>5]<-5
  
  DF.temp$T.ANX<-DF.temp$ANX
  DF.temp$T.ANX[DF.temp$T.ANX<1]<-1
  DF.temp$T.ANX[DF.temp$T.ANX>5]<-5
  
  DF.temp$T.DEP<-DF.temp$DEP
  DF.temp$T.DEP[DF.temp$T.DEP<1]<-1
  DF.temp$T.DEP[DF.temp$T.DEP>5]<-5
  
  DF.temp$T.ANG<-DF.temp$ANG
  DF.temp$T.ANG[DF.temp$T.ANG<1]<-1
  DF.temp$T.ANG[DF.temp$T.ANG>5]<-5
  
  DF.temp$T.TIRED<-DF.temp$TIRED
  DF.temp$T.TIRED[DF.temp$T.TIRED<1]<-1
  DF.temp$T.TIRED[DF.temp$T.TIRED>5]<-5
  
  DF.temp$c.NegEvnt<-as.numeric(as.character(DF.temp$NegEvnt_char))-DF.temp$prop_NegEvnt
  DF.temp$c.PosEvnt<-as.numeric(as.character(DF.temp$PosEvnt_char))-DF.temp$prop_PosEvnt
  
  DF.temp$NegEvnt_x_DN<-DF.temp$c.NegEvnt*DF.temp$DN_comb
  DF.temp$PosEvnt_x_DN<-DF.temp$c.PosEvnt*DF.temp$DN_comb
  
  dat.study2.list[[i]]<-DF.temp
}

###############################################################################
#Now for Modeling Section - Starting with Negative Mood (attempting to re-create Study 1 models)

#---------------------------------------------------------------------------------------------------------
#NEGATIVE AFFECT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Negative unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
log(mean(dat.study2$NEG, na.rm=TRUE))
hist(dat.study2$NEG)

NEG_ucm<-brms::brm_multiple(T.NEG~1+(1|ID), 
                            data = dat.study2.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                      set_prior('normal(0.3014526,2)', class='Intercept')),
                            warmup = 2000, 
                            iter = 3000, 
                            chains = 3,
                            control = list(adapt_delta=.99, 
                                           max_treedepth=15), 
                            save_model = paste0(stan.code, '/S2_ucm.stan'))
gc()

print(summary(NEG_ucm), digits = 5)
print(NEG_ucm$rhats, digits = 5)
plot(NEG_ucm)
pp_check(NEG_ucm)

#------------------------------------------------------------------------------ 
#Added this modeling section after meeting with AJS & JH on 1/25/19
#Goal is to solely look at negative events and DN intersection in predicting momentary mood variance
NEG_lv1_NDE<-brms::brm_multiple(T.NEG~1+c.NegEvnt+
                                    (1+c.NegEvnt|ID), 
                                  data = dat.study2.list, 
                                  family = 'lognormal',
                                  prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                            set_prior('normal(0.3014526,2)', class='Intercept'), 
                                            set_prior('normal(0,5)', class='b'), 
                                            set_prior('lkj(2)', class='cor')),
                                  warmup = 2000, 
                                  iter = 3000, 
                                  chains = 3,
                                  control = list(adapt_delta=.99, 
                                                 max_treedepth=15), 
                                  save_model = paste0(stan.code, '/S2_lv1NDE.stan'))
gc()

NEG_lv2_NDE_A<-brms::brm_multiple(T.NEG~1+c.NegEvnt+
                                      +prop_NegEvnt+
                                      (1+c.NegEvnt|ID), 
                                    data = dat.study2.list, 
                                    family = 'lognormal',
                                    prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                              set_prior('normal(0.3014526,2)', class='Intercept'), 
                                              set_prior('normal(0,5)', class='b'), 
                                              set_prior('lkj(2)', class='cor')),
                                    warmup = 2000, 
                                    iter = 3000, 
                                    chains = 3,
                                    control = list(adapt_delta=.99, 
                                                   max_treedepth=15), 
                                    save_model = paste0(stan.code, '/S2_lv2NDE_A.stan'))
gc()

NEG_lv2_NDE_B<-brms::brm_multiple(T.NEG~1+c.NegEvnt+
                                      DN_comb+
                                      (1+c.NegEvnt|ID), 
                                    data = dat.study2.list, 
                                    family = 'lognormal',
                                    prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                              set_prior('normal(0.3014526,2)', class='Intercept'), 
                                              set_prior('normal(0,5)', class='b'), 
                                              set_prior('lkj(2)', class='cor')),
                                    warmup = 2000, 
                                    iter = 3000, 
                                    chains = 3,
                                    control = list(adapt_delta=.99, 
                                                   max_treedepth=15), 
                                    save_model = paste0(stan.code, '/S2_lv2NDE_B.stan'))
gc()

NEG_lv2_NDE_C<-brms::brm_multiple(T.NEG~1+c.NegEvnt+
                                      DN_comb+c.NegEvnt+
                                      (1+c.NegEvnt|ID), 
                                    data = dat.study2.list, 
                                    family = 'lognormal',
                                    prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                              set_prior('normal(0.3014526,2)', class='Intercept'), 
                                              set_prior('normal(0,5)', class='b'), 
                                              set_prior('lkj(2)', class='cor')),
                                    warmup = 2000, 
                                    iter = 3000, 
                                    chains = 3,
                                    control = list(adapt_delta=.99, 
                                                   max_treedepth=15), 
                                    save_model = paste0(stan.code, '/S2_lv2NDE_C.stan'))
gc()
pp_check(NEG_NDE_cross, nsamples = 100)


NEG_NDE_cross<-brms::brm_multiple(T.NEG~1+c.NegEvnt+
                                    DN_comb+c.NegEvnt+
                                    c.NegEvnt:DN_comb+
                                    (1+c.NegEvnt|ID), 
                                    data = dat.study2.list, 
                                    family = 'lognormal',
                                    prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                              set_prior('normal(0.3014526,2)', class='Intercept'), 
                                              set_prior('normal(0,5)', class='b'), 
                                              set_prior('lkj(2)', class='cor')),
                                    warmup = 2000, 
                                    iter = 3000, 
                                    chains = 3,
                                    control = list(adapt_delta=.99, 
                                                   max_treedepth=15), 
                                    save_model = paste0(stan.code, '/S2_NDE_cross.stan'))
gc()
pp_check(NEG_NDE_cross, nsamples = 100)

#------------------------------------------------------------------------------
#Positive Daily Event Models
NEG_lv1_PDE<-brms::brm_multiple(T.NEG~1+c.PosEvnt+
                                  (1+c.PosEvnt|ID), 
                                data = dat.study2.list, 
                                family = 'lognormal',
                                prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                          set_prior('normal(0.3014526,2)', class='Intercept'), 
                                          set_prior('normal(0,5)', class='b'), 
                                          set_prior('lkj(2)', class='cor')),
                                warmup = 2000, 
                                iter = 3000, 
                                chains = 3,
                                control = list(adapt_delta=.99, 
                                               max_treedepth=15), 
                                save_model = paste0(stan.code, '/S2_lv1PDE.stan'))
gc()

NEG_lv2_PDE_A<-brms::brm_multiple(T.NEG~1+c.PosEvnt+
                                    +prop_PosEvnt+
                                    (1+c.PosEvnt|ID), 
                                  data = dat.study2.list, 
                                  family = 'lognormal',
                                  prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                            set_prior('normal(0.3014526,2)', class='Intercept'), 
                                            set_prior('normal(0,5)', class='b'), 
                                            set_prior('lkj(2)', class='cor')),
                                  warmup = 2000, 
                                  iter = 3000, 
                                  chains = 3,
                                  control = list(adapt_delta=.99, 
                                                 max_treedepth=15), 
                                  save_model = paste0(stan.code, '/S2_lv2PDE_A.stan'))
gc()

NEG_lv2_PDE_B<-brms::brm_multiple(T.NEG~1+c.PosEvnt+
                                    DN_comb+
                                    (1+c.PosEvnt|ID), 
                                  data = dat.study2.list, 
                                  family = 'lognormal',
                                  prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                            set_prior('normal(0.3014526,2)', class='Intercept'), 
                                            set_prior('normal(0,5)', class='b'), 
                                            set_prior('lkj(2)', class='cor')),
                                  warmup = 2000, 
                                  iter = 3000, 
                                  chains = 3,
                                  control = list(adapt_delta=.99, 
                                                 max_treedepth=15), 
                                  save_model = paste0(stan.code, '/S2_lv2PDE_B.stan'))
gc()

NEG_lv2_PDE_C<-brms::brm_multiple(T.NEG~1+c.PosEvnt+
                                    DN_comb+c.PosEvnt+
                                    (1+c.PosEvnt|ID), 
                                  data = dat.study2.list, 
                                  family = 'lognormal',
                                  prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                            set_prior('normal(0.3014526,2)', class='Intercept'), 
                                            set_prior('normal(0,5)', class='b'), 
                                            set_prior('lkj(2)', class='cor')),
                                  warmup = 2000, 
                                  iter = 3000, 
                                  chains = 3,
                                  control = list(adapt_delta=.99, 
                                                 max_treedepth=15), 
                                  save_model = paste0(stan.code, '/S2_lv2PDE_C.stan'))
gc()
pp_check(NEG_PDE_cross, nsamples = 100)


NEG_PDE_cross<-brms::brm_multiple(T.NEG~1+c.PosEvnt+
                                    DN_comb+c.PosEvnt+
                                    c.PosEvnt:DN_comb+
                                    (1+c.PosEvnt|ID), 
                                  data = dat.study2.list, 
                                  family = 'lognormal',
                                  prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                            set_prior('normal(0.3014526,2)', class='Intercept'), 
                                            set_prior('normal(0,5)', class='b'), 
                                            set_prior('lkj(2)', class='cor')),
                                  warmup = 2000, 
                                  iter = 3000, 
                                  chains = 3,
                                  control = list(adapt_delta=.99, 
                                                 max_treedepth=15), 
                                  save_model = paste0(stan.code, '/S2_PDE_cross.stan'))
gc()
pp_check(NEG_PDE_cross, nsamples = 100)

#Getting posterior predictive distribution (a model check for Negative Mood)
png(paste0(study2.graphics, '/S2_Supplemental_FigureX_final_ppdens_NA.png'), 
    res = 600, 
    units = 'in', 
    height = 6, 
    width = 6)
pp_check(Neg_cross, nsamples=100)+
  ggtitle('Study 1: Posterior Predictive Distribution for Negative Affect')+
  scale_x_continuous(limits=c(0,6))
dev.off()

#---------------------------------------------------------------------------------------------------------
#POSITIVE AFFECT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Positive unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
mean(dat.study2$POS, na.rm=TRUE)
mean(dat.study2$POS, na.rm=TRUE)
hist(dat.study2$POS)

POS_ucm<-brms::brm_multiple(T.POS~1+(1|ID), 
                            data = dat.study2.list, 
                            family = 'normal',
                            prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                      set_prior('normal(02.562117,2)', class='Intercept')),
                            warmup = 2000, 
                            iter = 3000, 
                            chains = 3,
                            control = list(adapt_delta=.99, 
                                           max_treedepth=15), 
                            save_model = paste0(stan.code, '/S2_ucm.stan'))
gc()

print(summary(POS_ucm), digits = 5)
print(POS_ucm$rhats, digits = 5)
plot(POS_ucm)
pp_check(POS_ucm)

#------------------------------------------------------------------------------ 
#Added this modeling section after meeting with AJS & JH on 1/25/19
#Goal is to solely look at negative events and DN intersection in predicting momentary mood variance
POS_lv1_NDE<-brms::brm_multiple(T.POS~1+c.NegEvnt+
                                  (1+c.NegEvnt|ID), 
                                data = dat.study2.list, 
                                family = 'normal',
                                prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                          set_prior('normal(02.562117,2)', class='Intercept'), 
                                          set_prior('normal(0,5)', class='b'), 
                                          set_prior('lkj(2)', class='cor')),
                                warmup = 2000, 
                                iter = 3000, 
                                chains = 3,
                                control = list(adapt_delta=.99, 
                                               max_treedepth=15), 
                                save_model = paste0(stan.code, '/S2_lv1NDE.stan'))
gc()

POS_lv2_NDE_A<-brms::brm_multiple(T.POS~1+c.NegEvnt+
                                    +prop_NegEvnt+
                                    (1+c.NegEvnt|ID), 
                                  data = dat.study2.list, 
                                  family = 'normal',
                                  prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                            set_prior('normal(02.562117,2)', class='Intercept'), 
                                            set_prior('normal(0,5)', class='b'), 
                                            set_prior('lkj(2)', class='cor')),
                                  warmup = 2000, 
                                  iter = 3000, 
                                  chains = 3,
                                  control = list(adapt_delta=.99, 
                                                 max_treedepth=15), 
                                  save_model = paste0(stan.code, '/S2_lv2NDE_A.stan'))
gc()

POS_lv2_NDE_B<-brms::brm_multiple(T.POS~1+c.NegEvnt+
                                    DN_comb+
                                    (1+c.NegEvnt|ID), 
                                  data = dat.study2.list, 
                                  family = 'normal',
                                  prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                            set_prior('normal(02.562117,2)', class='Intercept'), 
                                            set_prior('normal(0,5)', class='b'), 
                                            set_prior('lkj(2)', class='cor')),
                                  warmup = 2000, 
                                  iter = 3000, 
                                  chains = 3,
                                  control = list(adapt_delta=.99, 
                                                 max_treedepth=15), 
                                  save_model = paste0(stan.code, '/S2_lv2NDE_B.stan'))
gc()

POS_lv2_NDE_C<-brms::brm_multiple(T.POS~1+c.NegEvnt+
                                    DN_comb+c.NegEvnt+
                                    (1+c.NegEvnt|ID), 
                                  data = dat.study2.list, 
                                  family = 'normal',
                                  prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                            set_prior('normal(02.562117,2)', class='Intercept'), 
                                            set_prior('normal(0,5)', class='b'), 
                                            set_prior('lkj(2)', class='cor')),
                                  warmup = 2000, 
                                  iter = 3000, 
                                  chains = 3,
                                  control = list(adapt_delta=.99, 
                                                 max_treedepth=15), 
                                  save_model = paste0(stan.code, '/S2_lv2NDE_C.stan'))
gc()
pp_check(POS_NDE_cross, nsamples = 100)


POS_NDE_cross<-brms::brm_multiple(T.POS~1+c.NegEvnt+
                                    DN_comb+c.NegEvnt+
                                    c.NegEvnt:DN_comb+
                                    (1+c.NegEvnt|ID), 
                                  data = dat.study2.list, 
                                  family = 'normal',
                                  prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                            set_prior('normal(02.562117,2)', class='Intercept'), 
                                            set_prior('normal(0,5)', class='b'), 
                                            set_prior('lkj(2)', class='cor')),
                                  warmup = 2000, 
                                  iter = 3000, 
                                  chains = 3,
                                  control = list(adapt_delta=.99, 
                                                 max_treedepth=15), 
                                  save_model = paste0(stan.code, '/S2_NDE_cross.stan'))
gc()
pp_check(POS_NDE_cross, nsamples = 100)

#------------------------------------------------------------------------------
#Positive Daily Event Models
POS_lv1_PDE<-brms::brm_multiple(T.POS~1+c.PosEvnt+
                                  (1+c.PosEvnt|ID), 
                                data = dat.study2.list, 
                                family = 'normal',
                                prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                          set_prior('normal(02.562117,2)', class='Intercept'), 
                                          set_prior('normal(0,5)', class='b'), 
                                          set_prior('lkj(2)', class='cor')),
                                warmup = 2000, 
                                iter = 3000, 
                                chains = 3,
                                control = list(adapt_delta=.99, 
                                               max_treedepth=15), 
                                save_model = paste0(stan.code, '/S2_lv1PDE.stan'))
gc()

POS_lv2_PDE_A<-brms::brm_multiple(T.POS~1+c.PosEvnt+
                                    +prop_PosEvnt+
                                    (1+c.PosEvnt|ID), 
                                  data = dat.study2.list, 
                                  family = 'normal',
                                  prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                            set_prior('normal(02.562117,2)', class='Intercept'), 
                                            set_prior('normal(0,5)', class='b'), 
                                            set_prior('lkj(2)', class='cor')),
                                  warmup = 2000, 
                                  iter = 3000, 
                                  chains = 3,
                                  control = list(adapt_delta=.99, 
                                                 max_treedepth=15), 
                                  save_model = paste0(stan.code, '/S2_lv2PDE_A.stan'))
gc()

POS_lv2_PDE_B<-brms::brm_multiple(T.POS~1+c.PosEvnt+
                                    DN_comb+
                                    (1+c.PosEvnt|ID), 
                                  data = dat.study2.list, 
                                  family = 'normal',
                                  prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                            set_prior('normal(02.562117,2)', class='Intercept'), 
                                            set_prior('normal(0,5)', class='b'), 
                                            set_prior('lkj(2)', class='cor')),
                                  warmup = 2000, 
                                  iter = 3000, 
                                  chains = 3,
                                  control = list(adapt_delta=.99, 
                                                 max_treedepth=15), 
                                  save_model = paste0(stan.code, '/S2_lv2PDE_B.stan'))
gc()

POS_lv2_PDE_C<-brms::brm_multiple(T.POS~1+c.PosEvnt+
                                    DN_comb+c.PosEvnt+
                                    (1+c.PosEvnt|ID), 
                                  data = dat.study2.list, 
                                  family = 'normal',
                                  prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                            set_prior('normal(02.562117,2)', class='Intercept'), 
                                            set_prior('normal(0,5)', class='b'), 
                                            set_prior('lkj(2)', class='cor')),
                                  warmup = 2000, 
                                  iter = 3000, 
                                  chains = 3,
                                  control = list(adapt_delta=.99, 
                                                 max_treedepth=15), 
                                  save_model = paste0(stan.code, '/S2_lv2PDE_C.stan'))
gc()
pp_check(POS_PDE_cross, nsamples = 100)


POS_PDE_cross<-brms::brm_multiple(T.POS~1+c.PosEvnt+
                                    DN_comb+c.PosEvnt+
                                    c.PosEvnt:DN_comb+
                                    (1+c.PosEvnt|ID), 
                                  data = dat.study2.list, 
                                  family = 'normal',
                                  prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                            set_prior('normal(02.562117,2)', class='Intercept'), 
                                            set_prior('normal(0,5)', class='b'), 
                                            set_prior('lkj(2)', class='cor')),
                                  warmup = 2000, 
                                  iter = 3000, 
                                  chains = 3,
                                  control = list(adapt_delta=.99, 
                                                 max_treedepth=15), 
                                  save_model = paste0(stan.code, '/S2_PDE_cross.stan'))
gc()
pp_check(POS_PDE_cross, nsamples = 100)

#Getting posterior predictive distribution (a model check for Positive Mood)
png(paste0(study2.graphics, '/S2_Supplemental_FigureX_final_ppdens_NA.png'), 
    res = 600, 
    units = 'in', 
    height = 6, 
    width = 6)
pp_check(Neg_cross, nsamples=100)+
  ggtitle('Study 1: Posterior Predictive Distribution for Positive Affect')+
  scale_x_continuous(limits=c(0,6))
dev.off()

###############################################################################
#MAIN MODELS FOR POSITIVE AND NEGATIVE MOOD
###############################################################################

#------------------------------------------------------------------------------
#Negative Mood Models: (Final for main publication)
Neg_lv2_DN<-brms::brm_multiple(T.NEG~1+c.NegEvnt+c.PosEvnt+
                                 DN_comb+
                                 (1+c.NegEvnt+c.PosEvnt|ID), 
                               data = dat.study2.list, 
                               family = 'lognormal',
                               prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                         set_prior('normal(0.3014526,2)', class='Intercept'), 
                                         set_prior('normal(0,5)', class='b'), 
                                         set_prior('lkj(2)', class='cor')),
                               warmup = 2000, 
                               iter = 3000, 
                               chains = 3,
                               control = list(adapt_delta=.99, 
                                              max_treedepth=15), 
                               save_model = paste0(stan.code, '/S2_lv2_DN.stan'))
gc()

Neg_lv2_Evnt<-brms::brm_multiple(T.NEG~1+c.NegEvnt+c.PosEvnt+
                                   prop_NegEvnt+prop_PosEvnt+
                                   (1+c.NegEvnt+c.PosEvnt|ID), 
                                 data = dat.study2.list, 
                                 family = 'lognormal',
                                 prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                           set_prior('normal(0.3014526,2)', class='Intercept'), 
                                           set_prior('normal(0,5)', class='b'), 
                                           set_prior('lkj(2)', class='cor')),
                                 warmup = 2000, 
                                 iter = 3000, 
                                 chains = 3,
                                 control = list(adapt_delta=.99, 
                                                max_treedepth=15), 
                                 save_model = paste0(stan.code, '/S2_lv2_Evnt.stan'))
gc()

Neg_lv2_All<-brms::brm_multiple(T.NEG~1+c.NegEvnt+c.PosEvnt+
                                  prop_NegEvnt+prop_PosEvnt+DN_comb+
                                  (1+c.NegEvnt+c.PosEvnt|ID), 
                                data = dat.study2.list, 
                                family = 'lognormal',
                                prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                          set_prior('normal(0.3014526,2)', class='Intercept'), 
                                          set_prior('normal(0,5)', class='b'), 
                                          set_prior('lkj(2)', class='cor')),
                                warmup = 2000, 
                                iter = 3000, 
                                chains = 3,
                                control = list(adapt_delta=.99, 
                                               max_treedepth=15), 
                                save_model = paste0(stan.code, '/S2_lv2_All.stan'))
gc()

Neg_cross<-brms::brm_multiple(T.NEG~1+c.NegEvnt+c.PosEvnt+
                                prop_NegEvnt+prop_PosEvnt+DN_comb+
                                c.PosEvnt:DN_comb+c.NegEvnt:DN_comb+
                                (1+c.NegEvnt+c.PosEvnt|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                        set_prior('normal(0.3014526,2)', class='Intercept'), 
                                        set_prior('normal(0,5)', class='b'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 2000, 
                              iter = 3000, 
                              chains = 3,
                              control = list(adapt_delta=.99, 
                                             max_treedepth=15), 
                              save_model = paste0(stan.code, '/S2_cross.stan'))
gc()
pp_check(Neg_cross)

#---
Pos_lv2_DN<-brms::brm_multiple(T.POS~1+c.NegEvnt+c.PosEvnt+
                                 DN_comb+
                                 (1+c.NegEvnt+c.PosEvnt|ID), 
                               data = dat.study2.list, 
                               family = 'normal',
                               prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                         set_prior('normal(02.562117,2)', class='Intercept'), 
                                         set_prior('normal(0,5)', class='b'), 
                                         set_prior('lkj(2)', class='cor')),
                               warmup = 2000, 
                               iter = 3000, 
                               chains = 3,
                               control = list(adapt_delta=.99, 
                                              max_treedepth=15), 
                               save_model = paste0(stan.code, '/S2_lv2_DN.stan'))
gc()

Pos_lv2_Evnt<-brms::brm_multiple(T.POS~1+c.NegEvnt+c.PosEvnt+
                                   prop_NegEvnt+prop_PosEvnt+
                                   (1+c.NegEvnt+c.PosEvnt|ID), 
                                 data = dat.study2.list, 
                                 family = 'normal',
                                 prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                           set_prior('normal(02.562117,2)', class='Intercept'), 
                                           set_prior('normal(0,5)', class='b'), 
                                           set_prior('lkj(2)', class='cor')),
                                 warmup = 2000, 
                                 iter = 3000, 
                                 chains = 3,
                                 control = list(adapt_delta=.99, 
                                                max_treedepth=15), 
                                 save_model = paste0(stan.code, '/S2_lv2_Evnt.stan'))
gc()

Pos_lv2_All<-brms::brm_multiple(T.POS~1+c.NegEvnt+c.PosEvnt+
                                  prop_NegEvnt+prop_PosEvnt+DN_comb+
                                  (1+c.NegEvnt+c.PosEvnt|ID), 
                                data = dat.study2.list, 
                                family = 'normal',
                                prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                          set_prior('normal(02.562117,2)', class='Intercept'), 
                                          set_prior('normal(0,5)', class='b'), 
                                          set_prior('lkj(2)', class='cor')),
                                warmup = 2000, 
                                iter = 3000, 
                                chains = 3,
                                control = list(adapt_delta=.99, 
                                               max_treedepth=15), 
                                save_model = paste0(stan.code, '/S2_lv2_All.stan'))
gc()

Pos_cross<-brms::brm_multiple(T.POS~1+c.NegEvnt+c.PosEvnt+
                                prop_NegEvnt+prop_PosEvnt+DN_comb+
                                c.PosEvnt:DN_comb+c.NegEvnt:DN_comb+
                                (1+c.NegEvnt+c.PosEvnt|ID), 
                              data = dat.study2.list, 
                              family = 'normal',
                              prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                        set_prior('normal(02.562117,2)', class='Intercept'), 
                                        set_prior('normal(0,5)', class='b'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 2000, 
                              iter = 3000, 
                              chains = 3,
                              control = list(adapt_delta=.99, 
                                             max_treedepth=15), 
                              save_model = paste0(stan.code, '/S2_cross.stan'))
gc()

###############################################################################
#Variance Decomposition area

#1. - The Rights and Sterba Function
#Note I removed the plot so that I could 
r2MLM <- function(data,within_covs,between_covs,random_covs,
                  gamma_w,gamma_b,Tau,sigma2,has_intercept=T,clustermeancentered=T){
  #browser()
  if(has_intercept==T){
    if(length(gamma_b)>1) gamma <- c(1,gamma_w,gamma_b[2:length(gamma_b)])
    if(length(gamma_b)==1) gamma <- c(1,gamma_w)
    if(is.null(within_covs)==T) gamma_w <- 0
  }
  if(has_intercept==F){
    gamma <- c(gamma_w,gamma_b)
    if(is.null(within_covs)==T) gamma_w <- 0
    if(is.null(between_covs)==T) gamma_b <- 0
  }
  if(is.null(gamma)) gamma <- 0
  ##compute phi
  phi <- var(cbind(1,data[,c(within_covs)],data[,c(between_covs)]),na.rm=T)
  if(has_intercept==F) phi <- var(cbind(data[,c(within_covs)],data[,c(between_covs)]),na.rm=T)
  if(is.null(within_covs)==T & is.null(within_covs)==T & has_intercept==F) phi <- 0
  phi_w <- var(data[,within_covs],na.rm=T)
  if(is.null(within_covs)==T) phi_w <- 0
  phi_b <- var(cbind(1,data[,between_covs]),na.rm=T)
  if(is.null(between_covs)==T) phi_b <- 0
  ##compute psi and kappa
  var_randomcovs <- var(cbind(1,data[,c(random_covs)]),na.rm=T)
  if(length(Tau)>1) psi <- matrix(c(diag(Tau)),ncol=1)
  if(length(Tau)==1) psi <- Tau
  if(length(Tau)>1) kappa <- matrix(c(Tau[lower.tri(Tau)==TRUE]),ncol=1)
  if(length(Tau)==1) kappa <- 0
  v <- matrix(c(diag(var_randomcovs)),ncol=1)
  r <- matrix(c(var_randomcovs[lower.tri(var_randomcovs)==TRUE]),ncol=1)
  if(is.null(random_covs)==TRUE){
    v <- 0
    r <- 0
    m <- matrix(1,ncol=1)
  }
  if(length(random_covs)>0) m <- matrix(c(colMeans(cbind(1,data[,c(random_covs)]),na.rm=T)),ncol=1)
  ##total variance
  totalvar_notdecomp <- t(v)%*%psi + 2*(t(r)%*%kappa) + t(gamma)%*%phi%*%gamma + t(m)%*%Tau%*%m + sigma2
  totalwithinvar <- (t(gamma_w)%*%phi_w%*%gamma_w) + (t(v)%*%psi + 2*(t(r)%*%kappa)) + sigma2
  totalbetweenvar <- (t(gamma_b)%*%phi_b%*%gamma_b) + Tau[1]
  totalvar <- totalwithinvar + totalbetweenvar
  ##total decomp
  decomp_fixed_notdecomp <- (t(gamma)%*%phi%*%gamma) / totalvar
  decomp_fixed_within <- (t(gamma_w)%*%phi_w%*%gamma_w) / totalvar
  decomp_fixed_between <- (t(gamma_b)%*%phi_b%*%gamma_b) / totalvar
  decomp_fixed <- decomp_fixed_within + decomp_fixed_between
  decomp_varslopes <- (t(v)%*%psi + 2*(t(r)%*%kappa)) / totalvar
  decomp_varmeans <- (t(m)%*%Tau%*%m) / totalvar
  decomp_sigma <- sigma2/totalvar
  ##within decomp
  decomp_fixed_within_w <- (t(gamma_w)%*%phi_w%*%gamma_w) / totalwithinvar
  decomp_varslopes_w <- (t(v)%*%psi + 2*(t(r)%*%kappa)) / totalwithinvar
  decomp_sigma_w <- sigma2/totalwithinvar
  ##between decomp
  decomp_fixed_between_b <- (t(gamma_b)%*%phi_b%*%gamma_b) / totalbetweenvar
  decomp_varmeans_b <- Tau[1] / totalbetweenvar
  #NEW measures
  if (clustermeancentered==TRUE){
    R2_f <- decomp_fixed
    R2_f1 <- decomp_fixed_within
    R2_f2 <- decomp_fixed_between
    R2_fv <- decomp_fixed + decomp_varslopes
    R2_fvm <- decomp_fixed + decomp_varslopes + decomp_varmeans
    R2_v <- decomp_varslopes
    R2_m <- decomp_varmeans
    R2_f_w <- decomp_fixed_within_w
    R2_f_b <- decomp_fixed_between_b
    R2_fv_w <- decomp_fixed_within_w + decomp_varslopes_w
    R2_v_w <- decomp_varslopes_w
    R2_m_b <- decomp_varmeans_b
  }
  if (clustermeancentered==FALSE){
    R2_f <- decomp_fixed_notdecomp
    R2_fv <- decomp_fixed_notdecomp + decomp_varslopes
    R2_fvm <- decomp_fixed_notdecomp + decomp_varslopes + decomp_varmeans
    R2_v <- decomp_varslopes
    R2_m <- decomp_varmeans
  }
  if(clustermeancentered==TRUE){
    decomp_table <- matrix(c(decomp_fixed_within,decomp_fixed_between,decomp_varslopes,decomp_varmeans,decomp_sigma,
                             decomp_fixed_within_w,"NA",decomp_varslopes_w,"NA",decomp_sigma_w,
                             "NA",decomp_fixed_between_b,"NA",decomp_varmeans_b,"NA"),ncol=3)
    rownames(decomp_table) <- c("fixed, within","fixed, between","slope variation","mean variation","sigma2")
    colnames(decomp_table) <- c("total","within","between")
    R2_table <- matrix(c(R2_f1,R2_f2,R2_v,R2_m,R2_f,R2_fv,R2_fvm,
                         R2_f_w,"NA",R2_v_w,"NA","NA",R2_fv_w,"NA",
                         "NA",R2_f_b,"NA",R2_m_b,"NA","NA","NA")
                       ,ncol=3)
    rownames(R2_table) <- c("f1","f2","v","m","f","fv","fvm")
    colnames(R2_table) <- c("total","within","between")
  }
  ##barchart
  if(clustermeancentered==TRUE){
    contributions_stacked <- matrix(c(decomp_fixed_within,decomp_fixed_between,decomp_varslopes,decomp_varmeans,decomp_sigma,
                                      decomp_fixed_within_w,0,decomp_varslopes_w,0,decomp_sigma_w,
                                      0,decomp_fixed_between_b,0,decomp_varmeans_b,0),5,3)
    colnames(contributions_stacked) <- c("total","within","between")
    rownames(contributions_stacked) <- c("fixed slopes (within)",
                                         "fixed slopes (between)",
                                         "slope variation (within)",
                                         "intercept variation (between)",
                                         "residual (within)")
    #barplot(contributions_stacked, main="Decomposition", horiz=FALSE,
    #        ylim=c(0,1),col=c("darkred","steelblue","darkred","midnightblue","white"),ylab="proportion of variance",
    #        density=c(NA,NA,30,40,NA),angle=c(0,45,0,135,0),xlim=c(0,1),width=c(.3,.3))
    #legend(.30,-.1,legend=rownames(contributions_stacked),fill=c("darkred","steelblue","darkred","midnightblue","white"),
    #       cex=.7, pt.cex = 1,xpd=T,density=c(NA,NA,30,40,NA),angle=c(0,45,0,135,0))
  }
  if(clustermeancentered==FALSE){
    decomp_table <- matrix(c(decomp_fixed_notdecomp,decomp_varslopes,decomp_varmeans,decomp_sigma),ncol=1)
    rownames(decomp_table) <- c("fixed","slope variation","mean variation","sigma2")
    colnames(decomp_table) <- c("total")
    R2_table <- matrix(c(R2_f,R2_v,R2_m,R2_fv,R2_fvm),ncol=1)
    rownames(R2_table) <- c("f","v","m","fv","fvm")
    colnames(R2_table) <- c("total")
    ##barchar
    contributions_stacked <- matrix(c(decomp_fixed_notdecomp,decomp_varslopes,decomp_varmeans,decomp_sigma),4,1)
    colnames(contributions_stacked) <- c("total")
    rownames(contributions_stacked) <- c("fixed slopes",
                                         "slope variation",
                                         "intercept variation",
                                         "residual")
    #barplot(contributions_stacked, main="Decomposition", horiz=FALSE,
    #        ylim=c(0,1),col=c("darkblue","darkblue","darkblue","white"),ylab="proportion of variance",
    #        density=c(NA,30,40,NA),angle=c(0,0,135,0),xlim=c(0,1),width=c(.6))
    #legend(.30,-.1,legend=rownames(contributions_stacked),fill=c("darkblue","darkblue","darkblue","white"),
    #       cex=.7, pt.cex = 1,xpd=TRUE,density=c(NA,30,40,NA),angle=c(0,0,135,0))
  }
  Output <- list(noquote(decomp_table),noquote(R2_table))
  names(Output) <- c("Decompositions","R2s")
  return(Output)
}

#==============================================================================
#NEGATIVE MOOD MODELS - Variance Decomposition
#==============================================================================
#Will need the grand intercept from the null model (on the link scale)
beta_00<-fixef(NEG_ucm)[1,1]
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
#First creating a data set that averages across imputed values:

#Attempting Variance Partionining for DN alone
within_cov<-c(24,25)                    #Columns with group-mean centered predictors
between_cov<-c(22)                       #Columns with between-subject predictors
random_cov<-c(24,25)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Neg_lv2_DN, pars = 'b_Intercept') 
DN = posterior_samples(Neg_lv2_DN, pars = 'b_DN_comb') 
NegEvnt = posterior_samples(Neg_lv2_DN, pars = 'b_c.NegEvnt') 
PosEvnt = posterior_samples(Neg_lv2_DN, pars = 'b_c.PosEvnt') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Neg_lv2_DN, pars = 'sd_ID__Intercept')^2
NegEvnt_var = posterior_samples(Neg_lv2_DN, pars = 'sd_ID__c.NegEvnt')^2 
PosEvnt_var = posterior_samples(Neg_lv2_DN, pars = 'sd_ID__c.PosEvnt')^2

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_NegEvnt = posterior_samples(Neg_lv2_DN, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_lv2_DN, pars = 'sd_ID__c.NegEvnt')*
  posterior_samples(Neg_lv2_DN, pars = 'cor_ID__Intercept__c.NegEvnt')

cov_Int_PosEvnt = posterior_samples(Neg_lv2_DN, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_lv2_DN, pars = 'sd_ID__c.PosEvnt')*
  posterior_samples(Neg_lv2_DN, pars = 'cor_ID__Intercept__c.PosEvnt') 

cov_NegEvnt_PosEvnt = posterior_samples(Neg_lv2_DN, pars = 'sd_ID__c.NegEvnt')*
  posterior_samples(Neg_lv2_DN, pars = 'sd_ID__c.PosEvnt')*
  posterior_samples(Neg_lv2_DN, pars = 'cor_ID__c.NegEvnt__c.PosEvnt')

#Getting level 1 error variance
sigma<-log(1+posterior_samples(Neg_lv2_DN, pars = 'sigma')/beta_00)

post_samples<-data.frame(Intercept,
                         DN, 
                         NegEvnt, 
                         PosEvnt, 
                         Int_var, 
                         NegEvnt_var, 
                         PosEvnt_var, 
                         cov_Int_NegEvnt, 
                         cov_Int_PosEvnt, 
                         cov_NegEvnt_PosEvnt, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'DN', 
                          'NegEvnt', 
                          'PosEvnt', 
                          'Int_var', 
                          'NegEvnt_var', 
                          'PosEvnt_var',
                          'cov_Int_NegEvnt', 
                          'cov_Int_PosEvnt', 
                          'cov_NegEvnt_PosEvnt', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

DN_alone<-list(between_var=vector(), 
               within_var=vector(), 
               between_All_tot=vector(),
               between_All_btw=vector(),
               between_res_btw=vector(),
               within_fix_wthn=vector(),
               within_fix_tot=vector(),
               within_slope_var_wthn=vector(),
               within_res_wthn=vector(), 
               within_unmod_tot=vector())

#Will be used to line up posterior distributions 
#Because  posterior parameter estimates are unique to each imputed data set
#Probably ways to make this more programmatically fluid... 
sampling_list<-list(imp1 = 1:3000, 
                    imp2 = 1:3000 + 3000, 
                    imp3 = 1:3000 + 6000, 
                    imp4 = 1:3000 + 9000,
                    imp5 = 1:3000 + 12000,
                    imp6 = 1:3000 + 15000,
                    imp7 = 1:3000 + 18000,
                    imp8 = 1:3000 + 21000,
                    imp9 = 1:3000 + 24000,
                    imp10 = 1:3000 + 27000,
                    imp11 = 1:3000 + 30000,
                    imp12 = 1:3000 + 33000,
                    imp13 = 1:3000 + 36000,
                    imp14 = 1:3000 + 39000,
                    imp15 = 1:3000 + 42000,
                    imp16 = 1:3000 + 45000,
                    imp17 = 1:3000 + 48000,
                    imp18 = 1:3000 + 51000,
                    imp19 = 1:3000 + 54000,
                    imp20 = 1:3000 + 57000)

for(i in 1:length(dat.study2.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  D<-sample(sampling_list[[i]], size = 1000, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('NegEvnt', 'PosEvnt')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'DN')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_NegEvnt', 'cov_Int_PosEvnt')], 
           post_samples[D[d], c('cov_Int_NegEvnt', 'NegEvnt_var', 'cov_NegEvnt_PosEvnt')], 
           post_samples[D[d], c('cov_Int_PosEvnt', 'cov_NegEvnt_PosEvnt', 'PosEvnt_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=3)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_DN<-r2MLM(data=dat.study2.list[[i]], 
                    within_covs = within_cov, 
                    between_covs = between_cov, 
                    random_covs = random_cov, 
                    gamma_w = Gamma_w, 
                    gamma_b = Gamma_b, 
                    Tau = tau, 
                    sigma2 = Sigma2, 
                    has_intercept = TRUE, 
                    clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    DN_alone$between_var<-c(DN_alone$between_var, as.numeric(r2mlm_DN$Decompositions['fixed, between', 'total'])+
                              as.numeric(r2mlm_DN$Decompositions['mean variation', 'total']))
    DN_alone$within_var<-c(DN_alone$within_var, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                             as.numeric(r2mlm_DN$Decompositions['slope variation', 'total'])+
                             as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    DN_alone$between_All_tot<-c(DN_alone$between_All_tot, as.numeric(r2mlm_DN$R2s['f2', 'total']))
    DN_alone$between_All_btw<-c(DN_alone$between_All_btw, as.numeric(r2mlm_DN$R2s['f2', 'between']))
    DN_alone$between_res_btw<-c(DN_alone$between_res_btw, as.numeric(r2mlm_DN$R2s['m', 'between']))
    DN_alone$within_fix_wthn<-c(DN_alone$within_fix_wthn, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'within']))
    DN_alone$within_fix_tot<-c(DN_alone$within_fix_tot, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    DN_alone$within_slope_var_wthn<-c(DN_alone$within_slope_var_wthn, as.numeric(r2mlm_DN$Decompositions['slope variation', 'within']))
    DN_alone$within_res_wthn<-c(DN_alone$within_res_wthn, as.numeric(r2mlm_DN$Decompositions['sigma2', 'within']))
    DN_alone$within_unmod_tot<-c(DN_alone$within_unmod_tot, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                   as.numeric(r2mlm_DN$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}

hist(DN_alone$between_var)
mean(DN_alone$between_var)
hist(DN_alone$within_var)
mean(DN_alone$within_var)

hist(DN_alone$between_All_tot)
hist(DN_alone$between_All_btw)
hist(DN_alone$within_fix_tot)

summary(Neg_lv2_DN)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=
#Event only model (these are context effects on intercept)
within_cov<-c(24,25)                    #Columns with group-mean centered predictors
between_cov<-c(32,33)                       #Columns with between-subject predictors
random_cov<-c(24,25)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Neg_lv2_Evnt, pars = 'b_Intercept') 
prop_NegEvnt = posterior_samples(Neg_lv2_Evnt, pars = 'b_prop_NegEvnt')
prop_PosEvnt = posterior_samples(Neg_lv2_Evnt, pars = 'b_prop_PosEvnt')
NegEvnt = posterior_samples(Neg_lv2_Evnt, pars = 'b_c.NegEvnt') 
PosEvnt = posterior_samples(Neg_lv2_Evnt, pars = 'b_c.PosEvnt') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Neg_lv2_Evnt, pars = 'sd_ID__Intercept')^2
NegEvnt_var = posterior_samples(Neg_lv2_Evnt, pars = 'sd_ID__c.NegEvnt')^2 
PosEvnt_var = posterior_samples(Neg_lv2_Evnt, pars = 'sd_ID__c.PosEvnt')^2

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_NegEvnt = posterior_samples(Neg_lv2_Evnt, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_lv2_Evnt, pars = 'sd_ID__c.NegEvnt')*
  posterior_samples(Neg_lv2_Evnt, pars = 'cor_ID__Intercept__c.NegEvnt')

cov_Int_PosEvnt = posterior_samples(Neg_lv2_Evnt, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_lv2_Evnt, pars = 'sd_ID__c.PosEvnt')*
  posterior_samples(Neg_lv2_Evnt, pars = 'cor_ID__Intercept__c.PosEvnt') 

cov_NegEvnt_PosEvnt = posterior_samples(Neg_lv2_Evnt, pars = 'sd_ID__c.NegEvnt')*
  posterior_samples(Neg_lv2_Evnt, pars = 'sd_ID__c.PosEvnt')*
  posterior_samples(Neg_lv2_Evnt, pars = 'cor_ID__c.NegEvnt__c.PosEvnt')

#Getting level 1 error variance
sigma<-log(1+posterior_samples(Neg_lv2_Evnt, pars = 'sigma')/beta_00)

post_samples<-data.frame(Intercept,
                         prop_NegEvnt,
                         prop_PosEvnt,
                         NegEvnt, 
                         PosEvnt, 
                         Int_var, 
                         NegEvnt_var, 
                         PosEvnt_var, 
                         cov_Int_NegEvnt, 
                         cov_Int_PosEvnt, 
                         cov_NegEvnt_PosEvnt, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'prop_NegEvnt',
                          'prop_PosEvnt',
                          'NegEvnt', 
                          'PosEvnt', 
                          'Int_var', 
                          'NegEvnt_var', 
                          'PosEvnt_var',
                          'cov_Int_NegEvnt', 
                          'cov_Int_PosEvnt', 
                          'cov_NegEvnt_PosEvnt', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Evnt_alone<-list(between_var=vector(), 
               within_var=vector(), 
               between_All_tot=vector(),
               between_All_btw=vector(),
               between_res_btw=vector(),
               within_fix_wthn=vector(),
               within_fix_tot=vector(),
               within_slope_var_wthn=vector(),
               within_res_wthn=vector(), 
               within_unmod_tot=vector())

for(i in 1:length(dat.study2.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  D<-sample(sampling_list[[i]], size = 1000, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('NegEvnt', 'PosEvnt')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'prop_NegEvnt', 'prop_PosEvnt')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_NegEvnt', 'cov_Int_PosEvnt')], 
           post_samples[D[d], c('cov_Int_NegEvnt', 'NegEvnt_var', 'cov_NegEvnt_PosEvnt')], 
           post_samples[D[d], c('cov_Int_PosEvnt', 'cov_NegEvnt_PosEvnt', 'PosEvnt_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=3)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_DN<-r2MLM(data=dat.study2.list[[i]], 
                    within_covs = within_cov, 
                    between_covs = between_cov, 
                    random_covs = random_cov, 
                    gamma_w = Gamma_w, 
                    gamma_b = Gamma_b, 
                    Tau = tau, 
                    sigma2 = Sigma2, 
                    has_intercept = TRUE, 
                    clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Evnt_alone$between_var<-c(Evnt_alone$between_var, as.numeric(r2mlm_DN$Decompositions['fixed, between', 'total'])+
                              as.numeric(r2mlm_DN$Decompositions['mean variation', 'total']))
    Evnt_alone$within_var<-c(Evnt_alone$within_var, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                             as.numeric(r2mlm_DN$Decompositions['slope variation', 'total'])+
                             as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    Evnt_alone$between_All_tot<-c(Evnt_alone$between_All_tot, as.numeric(r2mlm_DN$R2s['f2', 'total']))
    Evnt_alone$between_All_btw<-c(Evnt_alone$between_All_btw, as.numeric(r2mlm_DN$R2s['f2', 'between']))
    Evnt_alone$between_res_btw<-c(Evnt_alone$between_res_btw, as.numeric(r2mlm_DN$R2s['m', 'between']))
    Evnt_alone$within_fix_wthn<-c(Evnt_alone$within_fix_wthn, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'within']))
    Evnt_alone$within_fix_tot<-c(Evnt_alone$within_fix_tot, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    Evnt_alone$within_slope_var_wthn<-c(Evnt_alone$within_slope_var_wthn, as.numeric(r2mlm_DN$Decompositions['slope variation', 'within']))
    Evnt_alone$within_res_wthn<-c(Evnt_alone$within_res_wthn, as.numeric(r2mlm_DN$Decompositions['sigma2', 'within']))
    Evnt_alone$within_unmod_tot<-c(Evnt_alone$within_unmod_tot, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                   as.numeric(r2mlm_DN$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}
hist(Evnt_alone$between_var)
mean(Evnt_alone$between_var)

hist(Evnt_alone$within_var)
mean(Evnt_alone$within_var)

hist(Evnt_alone$between_All_tot)
hist(Evnt_alone$between_All_btw)
hist(Evnt_alone$within_fix_tot)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=
#Level 2  Combined - Negative Mood 
within_cov<-c(24,25)                    #Columns with group-mean centered predictors
between_cov<-c(22,32,33)                       #Columns with between-subject predictors
random_cov<-c(24,25)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Neg_lv2_All, pars = 'b_Intercept') 
DN = posterior_samples(Neg_lv2_All, pars = 'b_DN_comb') 
prop_NegEvnt = posterior_samples(Neg_lv2_All, pars = 'b_prop_NegEvnt')
prop_PosEvnt = posterior_samples(Neg_lv2_All, pars = 'b_prop_PosEvnt')
NegEvnt = posterior_samples(Neg_lv2_All, pars = 'b_c.NegEvnt') 
PosEvnt = posterior_samples(Neg_lv2_All, pars = 'b_c.PosEvnt') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Neg_lv2_All, pars = 'sd_ID__Intercept')^2
NegEvnt_var = posterior_samples(Neg_lv2_All, pars = 'sd_ID__c.NegEvnt')^2 
PosEvnt_var = posterior_samples(Neg_lv2_All, pars = 'sd_ID__c.PosEvnt')^2

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_NegEvnt = posterior_samples(Neg_lv2_All, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_lv2_All, pars = 'sd_ID__c.NegEvnt')*
  posterior_samples(Neg_lv2_All, pars = 'cor_ID__Intercept__c.NegEvnt')

cov_Int_PosEvnt = posterior_samples(Neg_lv2_All, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_lv2_All, pars = 'sd_ID__c.PosEvnt')*
  posterior_samples(Neg_lv2_All, pars = 'cor_ID__Intercept__c.PosEvnt') 

cov_NegEvnt_PosEvnt = posterior_samples(Neg_lv2_All, pars = 'sd_ID__c.NegEvnt')*
  posterior_samples(Neg_lv2_All, pars = 'sd_ID__c.PosEvnt')*
  posterior_samples(Neg_lv2_All, pars = 'cor_ID__c.NegEvnt__c.PosEvnt')

#Getting level 1 error variance
sigma<-log(1+posterior_samples(Neg_lv2_All, pars = 'sigma')/beta_00)

post_samples<-data.frame(Intercept,
                         DN,
                         prop_NegEvnt,
                         prop_PosEvnt,
                         NegEvnt, 
                         PosEvnt, 
                         Int_var, 
                         NegEvnt_var, 
                         PosEvnt_var, 
                         cov_Int_NegEvnt, 
                         cov_Int_PosEvnt, 
                         cov_NegEvnt_PosEvnt, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          "DN",
                          'prop_NegEvnt',
                          'prop_PosEvnt',
                          'NegEvnt', 
                          'PosEvnt', 
                          'Int_var', 
                          'NegEvnt_var', 
                          'PosEvnt_var',
                          'cov_Int_NegEvnt', 
                          'cov_Int_PosEvnt', 
                          'cov_NegEvnt_PosEvnt', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

Lv2_All<-list(between_var=vector(), 
                 within_var=vector(), 
                 between_All_tot=vector(),
                 between_All_btw=vector(),
                 between_res_btw=vector(),
                 within_fix_wthn=vector(),
                 within_fix_tot=vector(),
                 within_slope_var_wthn=vector(),
                 within_res_wthn=vector(), 
                 within_unmod_tot=vector())

for(i in 1:length(dat.study2.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  D<-sample(sampling_list[[i]], size = 1000, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('NegEvnt', 'PosEvnt')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', "DN", 'prop_NegEvnt', 'prop_PosEvnt')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_NegEvnt', 'cov_Int_PosEvnt')], 
           post_samples[D[d], c('cov_Int_NegEvnt', 'NegEvnt_var', 'cov_NegEvnt_PosEvnt')], 
           post_samples[D[d], c('cov_Int_PosEvnt', 'cov_NegEvnt_PosEvnt', 'PosEvnt_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=3)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_DN<-r2MLM(data=dat.study2.list[[i]], 
                    within_covs = within_cov, 
                    between_covs = between_cov, 
                    random_covs = random_cov, 
                    gamma_w = Gamma_w, 
                    gamma_b = Gamma_b, 
                    Tau = tau, 
                    sigma2 = Sigma2, 
                    has_intercept = TRUE, 
                    clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Lv2_All$between_var<-c(Lv2_All$between_var, as.numeric(r2mlm_DN$Decompositions['fixed, between', 'total'])+
                                as.numeric(r2mlm_DN$Decompositions['mean variation', 'total']))
    Lv2_All$within_var<-c(Lv2_All$within_var, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                               as.numeric(r2mlm_DN$Decompositions['slope variation', 'total'])+
                               as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    Lv2_All$between_All_tot<-c(Lv2_All$between_All_tot, as.numeric(r2mlm_DN$R2s['f2', 'total']))
    Lv2_All$between_All_btw<-c(Lv2_All$between_All_btw, as.numeric(r2mlm_DN$R2s['f2', 'between']))
    Lv2_All$between_res_btw<-c(Lv2_All$between_res_btw, as.numeric(r2mlm_DN$R2s['m', 'between']))
    Lv2_All$within_fix_wthn<-c(Lv2_All$within_fix_wthn, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'within']))
    Lv2_All$within_fix_tot<-c(Lv2_All$within_fix_tot, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    Lv2_All$within_slope_var_wthn<-c(Lv2_All$within_slope_var_wthn, as.numeric(r2mlm_DN$Decompositions['slope variation', 'within']))
    Lv2_All$within_res_wthn<-c(Lv2_All$within_res_wthn, as.numeric(r2mlm_DN$Decompositions['sigma2', 'within']))
    Lv2_All$within_unmod_tot<-c(Lv2_All$within_unmod_tot, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                     as.numeric(r2mlm_DN$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}

hist(Lv2_All$between_var)
mean(Lv2_All$between_var)
hist(Lv2_All$within_var)
mean(Lv2_All$within_var)

hist(Lv2_All$between_All_tot)
hist(Lv2_All$between_All_btw)
hist(Lv2_All$within_fix_tot)

#------------------------------------------------------------------------------
btw_DN_uni_tot<-mean(Lv2_All$between_All_tot)-mean(Evnt_alone$between_All_tot)
btw_shared_tot<-mean(DN_alone$between_All_tot)-btw_DN_uni_tot
btw_context_uni_tot<-mean(Lv2_All$between_All_tot)-mean(DN_alone$between_All_tot)

sum(c(btw_DN_uni_tot, 
      btw_shared_tot,
      btw_context_uni_tot))

#The above total should equal: 
mean(Lv2_All$between_All_tot)

#Now getting the between subjects residual variance from final model
btw_res_tot<-mean(Lv2_All$between_var-Lv2_All$between_All_tot)
btw_res_tot+sum(c(btw_DN_uni_tot, 
                  btw_shared_tot,
                  btw_context_uni_tot))
#Final getting the within subjects variance: 
wthn_best_worst_tot<-mean(Lv2_All$within_fix_wthn)
wthn_res_tot<-mean(Lv2_All$within_var)-wthn_best_worst_tot

#Should be equal to 1
sum(c(wthn_best_worst_tot, 
      wthn_res_tot, 
      mean(Lv2_All$between_var)))

tot_btw<-mean(Lv2_All$between_var)
tot_wthn<-mean(Lv2_All$within_var)  

River_DF<-data.frame(N1 = c('DN',
                            'DN <--> Exposure',
                            'Exposure', 
                            'Unmodeled Between', 
                            'Momentary Events', 
                            'Unmodeled Within', 
                            'Total Between', 
                            'Total Within'), 
                     N2 = c('Total Between', 
                            'Total Between', 
                            'Total Between', 
                            'Total Between', 
                            'Total Within', 
                            'Total Within', 
                            'Total Variance', 
                            'Total Variance'), 
                     Value = c(btw_DN_uni_tot, 
                               btw_shared_tot,
                               btw_context_uni_tot, 
                               btw_res_tot, 
                               wthn_best_worst_tot, 
                               wthn_res_tot, 
                               tot_btw, 
                               tot_wthn), 
                     ID = 1:8)

River_DF$N1<-paste(River_DF$N1, 
                   '\n', 
                   paste0(round(River_DF$Value*100, 
                                digits = 2), '%'
                   )
)
River_DF$N2<-c(rep(River_DF$N1[7], 4),
               rep(River_DF$N1[8], 2), 
               rep('Total Variance', 2)
)

nodes<-data.frame(ID = c(River_DF$N1, 
                         'Total Variance'), 
                  x = c(1,1,1,1,1,1,2,2,3), 
                  y = c(0,1.5,3.25,5.25,8,10.5,2.5,7.5,5))

palette = c(paste0(brewer.pal(9, "Blues"), 90)[c(2, 4, 6, 8)], 
            paste0(brewer.pal(9, "Reds"), 90)[c(6, 8)], 
            paste0(brewer.pal(9, "Blues"), 90)[9], 
            paste0(brewer.pal(9, "Reds"), 90)[9], 
            paste0(brewer.pal(9, 'Purples'), 90)[9])

styles = lapply(nodes$y, function(n) {
  list(col = palette[n], lty = 0, textcol = "black")
})

#Not sure why color is not mapping correctly - trying to force it with this fuction... 
for(i in 1:length(palette)){
  styles[[i]]$col<-palette[i]
}

names(styles) = nodes$ID

riv<-makeRiver(nodes = nodes, 
               edges =  River_DF, 
               node_styles = styles)

#Creating Dataset to add in text - Will have to manipulate the y-values to get text to line up 

png(paste0(study2.graphics, '/S2_NegMood_river.png'), 
    units = 'in', 
    res = 1200, 
    height = 10, 
    width = 10)
riverplot(riv, 
          nodewidth = 3, 
          plot_area = .95)
title(ylab = 'Study 2 - Riverplot of Total Variance Decomposition Estimated from Negative Mood Models')
dev.off()

#==============================================================================
#POSITIVE MOOD MODELS - Variance Decomposition
#==============================================================================
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
#First creating a data set that averages across imputed values:

#Attempting Variance Partionining for DN alone
within_cov<-c(24,25)                    #Columns with group-mean centered predictors
between_cov<-c(22)                       #Columns with between-subject predictors
random_cov<-c(24,25)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Pos_lv2_DN, pars = 'b_Intercept') 
DN = posterior_samples(Pos_lv2_DN, pars = 'b_DN_comb') 
NegEvnt = posterior_samples(Pos_lv2_DN, pars = 'b_c.NegEvnt') 
PosEvnt = posterior_samples(Pos_lv2_DN, pars = 'b_c.PosEvnt') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Pos_lv2_DN, pars = 'sd_ID__Intercept')^2
NegEvnt_var = posterior_samples(Pos_lv2_DN, pars = 'sd_ID__c.NegEvnt')^2 
PosEvnt_var = posterior_samples(Pos_lv2_DN, pars = 'sd_ID__c.PosEvnt')^2

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_NegEvnt = posterior_samples(Pos_lv2_DN, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_lv2_DN, pars = 'sd_ID__c.NegEvnt')*
  posterior_samples(Pos_lv2_DN, pars = 'cor_ID__Intercept__c.NegEvnt')

cov_Int_PosEvnt = posterior_samples(Pos_lv2_DN, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_lv2_DN, pars = 'sd_ID__c.PosEvnt')*
  posterior_samples(Pos_lv2_DN, pars = 'cor_ID__Intercept__c.PosEvnt') 

cov_NegEvnt_PosEvnt = posterior_samples(Pos_lv2_DN, pars = 'sd_ID__c.NegEvnt')*
  posterior_samples(Pos_lv2_DN, pars = 'sd_ID__c.PosEvnt')*
  posterior_samples(Pos_lv2_DN, pars = 'cor_ID__c.NegEvnt__c.PosEvnt')

#Getting level 1 error variance
sigma<-posterior_samples(Pos_lv2_DN, pars = 'sigma')

post_samples<-data.frame(Intercept,
                         DN, 
                         NegEvnt, 
                         PosEvnt, 
                         Int_var, 
                         NegEvnt_var, 
                         PosEvnt_var, 
                         cov_Int_NegEvnt, 
                         cov_Int_PosEvnt, 
                         cov_NegEvnt_PosEvnt, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'DN', 
                          'NegEvnt', 
                          'PosEvnt', 
                          'Int_var', 
                          'NegEvnt_var', 
                          'PosEvnt_var',
                          'cov_Int_NegEvnt', 
                          'cov_Int_PosEvnt', 
                          'cov_NegEvnt_PosEvnt', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

P_DN_alone<-list(between_var=vector(), 
               within_var=vector(), 
               between_All_tot=vector(),
               between_All_btw=vector(),
               between_res_btw=vector(),
               within_fix_wthn=vector(),
               within_fix_tot=vector(),
               within_slope_var_wthn=vector(),
               within_res_wthn=vector(), 
               within_unmod_tot=vector())

for(i in 1:length(dat.study2.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  D<-sample(sampling_list[[i]], size = 1000, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('NegEvnt', 'PosEvnt')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'DN')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_NegEvnt', 'cov_Int_PosEvnt')], 
           post_samples[D[d], c('cov_Int_NegEvnt', 'NegEvnt_var', 'cov_NegEvnt_PosEvnt')], 
           post_samples[D[d], c('cov_Int_PosEvnt', 'cov_NegEvnt_PosEvnt', 'PosEvnt_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=3)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_DN<-r2MLM(data=dat.study2.list[[i]], 
                    within_covs = within_cov, 
                    between_covs = between_cov, 
                    random_covs = random_cov, 
                    gamma_w = Gamma_w, 
                    gamma_b = Gamma_b, 
                    Tau = tau, 
                    sigma2 = Sigma2, 
                    has_intercept = TRUE, 
                    clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    P_DN_alone$between_var<-c(P_DN_alone$between_var, as.numeric(r2mlm_DN$Decompositions['fixed, between', 'total'])+
                              as.numeric(r2mlm_DN$Decompositions['mean variation', 'total']))
    P_DN_alone$within_var<-c(P_DN_alone$within_var, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                             as.numeric(r2mlm_DN$Decompositions['slope variation', 'total'])+
                             as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    P_DN_alone$between_All_tot<-c(P_DN_alone$between_All_tot, as.numeric(r2mlm_DN$R2s['f2', 'total']))
    P_DN_alone$between_All_btw<-c(P_DN_alone$between_All_btw, as.numeric(r2mlm_DN$R2s['f2', 'between']))
    P_DN_alone$between_res_btw<-c(P_DN_alone$between_res_btw, as.numeric(r2mlm_DN$R2s['m', 'between']))
    P_DN_alone$within_fix_wthn<-c(P_DN_alone$within_fix_wthn, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'within']))
    P_DN_alone$within_fix_tot<-c(P_DN_alone$within_fix_tot, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    P_DN_alone$within_slope_var_wthn<-c(P_DN_alone$within_slope_var_wthn, as.numeric(r2mlm_DN$Decompositions['slope variation', 'within']))
    P_DN_alone$within_res_wthn<-c(P_DN_alone$within_res_wthn, as.numeric(r2mlm_DN$Decompositions['sigma2', 'within']))
    P_DN_alone$within_unmod_tot<-c(P_DN_alone$within_unmod_tot, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                   as.numeric(r2mlm_DN$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}

hist(P_DN_alone$between_var)
mean(P_DN_alone$between_var)
hist(P_DN_alone$within_var)
mean(P_DN_alone$within_var)

hist(P_DN_alone$between_All_tot)
hist(P_DN_alone$between_All_btw)
hist(P_DN_alone$within_fix_tot)

summary(Pos_lv2_DN)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=
#Event only model (these are context effects on intercept)
within_cov<-c(24,25)                    #Columns with group-mean centered predictors
between_cov<-c(32,33)                       #Columns with between-subject predictors
random_cov<-c(24,25)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Pos_lv2_Evnt, pars = 'b_Intercept') 
prop_NegEvnt = posterior_samples(Pos_lv2_Evnt, pars = 'b_prop_NegEvnt')
prop_PosEvnt = posterior_samples(Pos_lv2_Evnt, pars = 'b_prop_PosEvnt')
NegEvnt = posterior_samples(Pos_lv2_Evnt, pars = 'b_c.NegEvnt') 
PosEvnt = posterior_samples(Pos_lv2_Evnt, pars = 'b_c.PosEvnt') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Pos_lv2_Evnt, pars = 'sd_ID__Intercept')^2
NegEvnt_var = posterior_samples(Pos_lv2_Evnt, pars = 'sd_ID__c.NegEvnt')^2 
PosEvnt_var = posterior_samples(Pos_lv2_Evnt, pars = 'sd_ID__c.PosEvnt')^2

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_NegEvnt = posterior_samples(Pos_lv2_Evnt, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_lv2_Evnt, pars = 'sd_ID__c.NegEvnt')*
  posterior_samples(Pos_lv2_Evnt, pars = 'cor_ID__Intercept__c.NegEvnt')

cov_Int_PosEvnt = posterior_samples(Pos_lv2_Evnt, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_lv2_Evnt, pars = 'sd_ID__c.PosEvnt')*
  posterior_samples(Pos_lv2_Evnt, pars = 'cor_ID__Intercept__c.PosEvnt') 

cov_NegEvnt_PosEvnt = posterior_samples(Pos_lv2_Evnt, pars = 'sd_ID__c.NegEvnt')*
  posterior_samples(Pos_lv2_Evnt, pars = 'sd_ID__c.PosEvnt')*
  posterior_samples(Pos_lv2_Evnt, pars = 'cor_ID__c.NegEvnt__c.PosEvnt')

#Getting level 1 error variance
sigma<-log(1+posterior_samples(Pos_lv2_Evnt, pars = 'sigma')/beta_00)

post_samples<-data.frame(Intercept,
                         prop_NegEvnt,
                         prop_PosEvnt,
                         NegEvnt, 
                         PosEvnt, 
                         Int_var, 
                         NegEvnt_var, 
                         PosEvnt_var, 
                         cov_Int_NegEvnt, 
                         cov_Int_PosEvnt, 
                         cov_NegEvnt_PosEvnt, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'prop_NegEvnt',
                          'prop_PosEvnt',
                          'NegEvnt', 
                          'PosEvnt', 
                          'Int_var', 
                          'NegEvnt_var', 
                          'PosEvnt_var',
                          'cov_Int_NegEvnt', 
                          'cov_Int_PosEvnt', 
                          'cov_NegEvnt_PosEvnt', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

P_Evnt_alone<-list(between_var=vector(), 
                 within_var=vector(), 
                 between_All_tot=vector(),
                 between_All_btw=vector(),
                 between_res_btw=vector(),
                 within_fix_wthn=vector(),
                 within_fix_tot=vector(),
                 within_slope_var_wthn=vector(),
                 within_res_wthn=vector(), 
                 within_unmod_tot=vector())

for(i in 1:length(dat.study2.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  D<-sample(sampling_list[[i]], size = 1000, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('NegEvnt', 'PosEvnt')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'prop_NegEvnt', 'prop_PosEvnt')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_NegEvnt', 'cov_Int_PosEvnt')], 
           post_samples[D[d], c('cov_Int_NegEvnt', 'NegEvnt_var', 'cov_NegEvnt_PosEvnt')], 
           post_samples[D[d], c('cov_Int_PosEvnt', 'cov_NegEvnt_PosEvnt', 'PosEvnt_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=3)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_DN<-r2MLM(data=dat.study2.list[[i]], 
                    within_covs = within_cov, 
                    between_covs = between_cov, 
                    random_covs = random_cov, 
                    gamma_w = Gamma_w, 
                    gamma_b = Gamma_b, 
                    Tau = tau, 
                    sigma2 = Sigma2, 
                    has_intercept = TRUE, 
                    clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    P_Evnt_alone$between_var<-c(P_Evnt_alone$between_var, as.numeric(r2mlm_DN$Decompositions['fixed, between', 'total'])+
                                as.numeric(r2mlm_DN$Decompositions['mean variation', 'total']))
    P_Evnt_alone$within_var<-c(P_Evnt_alone$within_var, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                               as.numeric(r2mlm_DN$Decompositions['slope variation', 'total'])+
                               as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    P_Evnt_alone$between_All_tot<-c(P_Evnt_alone$between_All_tot, as.numeric(r2mlm_DN$R2s['f2', 'total']))
    P_Evnt_alone$between_All_btw<-c(P_Evnt_alone$between_All_btw, as.numeric(r2mlm_DN$R2s['f2', 'between']))
    P_Evnt_alone$between_res_btw<-c(P_Evnt_alone$between_res_btw, as.numeric(r2mlm_DN$R2s['m', 'between']))
    P_Evnt_alone$within_fix_wthn<-c(P_Evnt_alone$within_fix_wthn, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'within']))
    P_Evnt_alone$within_fix_tot<-c(P_Evnt_alone$within_fix_tot, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    P_Evnt_alone$within_slope_var_wthn<-c(P_Evnt_alone$within_slope_var_wthn, as.numeric(r2mlm_DN$Decompositions['slope variation', 'within']))
    P_Evnt_alone$within_res_wthn<-c(P_Evnt_alone$within_res_wthn, as.numeric(r2mlm_DN$Decompositions['sigma2', 'within']))
    P_Evnt_alone$within_unmod_tot<-c(P_Evnt_alone$within_unmod_tot, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                     as.numeric(r2mlm_DN$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}
hist(P_Evnt_alone$between_var)
mean(P_Evnt_alone$between_var)

hist(P_Evnt_alone$within_var)
mean(P_Evnt_alone$within_var)

hist(P_Evnt_alone$between_All_tot)
hist(P_Evnt_alone$between_All_btw)
hist(P_Evnt_alone$within_fix_tot)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=
#Level 2  Combined - Negative Mood 
within_cov<-c(24,25)                    #Columns with group-mean centered predictors
between_cov<-c(22,32,33)                       #Columns with between-subject predictors
random_cov<-c(24,25)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Pos_lv2_All, pars = 'b_Intercept') 
DN = posterior_samples(Pos_lv2_All, pars = 'b_DN_comb') 
prop_NegEvnt = posterior_samples(Pos_lv2_All, pars = 'b_prop_NegEvnt')
prop_PosEvnt = posterior_samples(Pos_lv2_All, pars = 'b_prop_PosEvnt')
NegEvnt = posterior_samples(Pos_lv2_All, pars = 'b_c.NegEvnt') 
PosEvnt = posterior_samples(Pos_lv2_All, pars = 'b_c.PosEvnt') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Pos_lv2_All, pars = 'sd_ID__Intercept')^2
NegEvnt_var = posterior_samples(Pos_lv2_All, pars = 'sd_ID__c.NegEvnt')^2 
PosEvnt_var = posterior_samples(Pos_lv2_All, pars = 'sd_ID__c.PosEvnt')^2

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_NegEvnt = posterior_samples(Pos_lv2_All, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_lv2_All, pars = 'sd_ID__c.NegEvnt')*
  posterior_samples(Pos_lv2_All, pars = 'cor_ID__Intercept__c.NegEvnt')

cov_Int_PosEvnt = posterior_samples(Pos_lv2_All, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_lv2_All, pars = 'sd_ID__c.PosEvnt')*
  posterior_samples(Pos_lv2_All, pars = 'cor_ID__Intercept__c.PosEvnt') 

cov_NegEvnt_PosEvnt = posterior_samples(Pos_lv2_All, pars = 'sd_ID__c.NegEvnt')*
  posterior_samples(Pos_lv2_All, pars = 'sd_ID__c.PosEvnt')*
  posterior_samples(Pos_lv2_All, pars = 'cor_ID__c.NegEvnt__c.PosEvnt')

#Getting level 1 error variance
sigma<-log(1+posterior_samples(Pos_lv2_All, pars = 'sigma')/beta_00)

post_samples<-data.frame(Intercept,
                         DN,
                         prop_NegEvnt,
                         prop_PosEvnt,
                         NegEvnt, 
                         PosEvnt, 
                         Int_var, 
                         NegEvnt_var, 
                         PosEvnt_var, 
                         cov_Int_NegEvnt, 
                         cov_Int_PosEvnt, 
                         cov_NegEvnt_PosEvnt, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          "DN",
                          'prop_NegEvnt',
                          'prop_PosEvnt',
                          'NegEvnt', 
                          'PosEvnt', 
                          'Int_var', 
                          'NegEvnt_var', 
                          'PosEvnt_var',
                          'cov_Int_NegEvnt', 
                          'cov_Int_PosEvnt', 
                          'cov_NegEvnt_PosEvnt', 
                          'sigma')

#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets

P_Lv2_All<-list(between_var=vector(), 
              within_var=vector(), 
              between_All_tot=vector(),
              between_All_btw=vector(),
              between_res_btw=vector(),
              within_fix_wthn=vector(),
              within_fix_tot=vector(),
              within_slope_var_wthn=vector(),
              within_res_wthn=vector(), 
              within_unmod_tot=vector())

for(i in 1:length(dat.study2.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  D<-sample(sampling_list[[i]], size = 1000, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('NegEvnt', 'PosEvnt')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', "DN", 'prop_NegEvnt', 'prop_PosEvnt')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_NegEvnt', 'cov_Int_PosEvnt')], 
           post_samples[D[d], c('cov_Int_NegEvnt', 'NegEvnt_var', 'cov_NegEvnt_PosEvnt')], 
           post_samples[D[d], c('cov_Int_PosEvnt', 'cov_NegEvnt_PosEvnt', 'PosEvnt_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=3)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_DN<-r2MLM(data=dat.study2.list[[i]], 
                    within_covs = within_cov, 
                    between_covs = between_cov, 
                    random_covs = random_cov, 
                    gamma_w = Gamma_w, 
                    gamma_b = Gamma_b, 
                    Tau = tau, 
                    sigma2 = Sigma2, 
                    has_intercept = TRUE, 
                    clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    P_Lv2_All$between_var<-c(P_Lv2_All$between_var, as.numeric(r2mlm_DN$Decompositions['fixed, between', 'total'])+
                             as.numeric(r2mlm_DN$Decompositions['mean variation', 'total']))
    P_Lv2_All$within_var<-c(P_Lv2_All$within_var, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                            as.numeric(r2mlm_DN$Decompositions['slope variation', 'total'])+
                            as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    P_Lv2_All$between_All_tot<-c(P_Lv2_All$between_All_tot, as.numeric(r2mlm_DN$R2s['f2', 'total']))
    P_Lv2_All$between_All_btw<-c(P_Lv2_All$between_All_btw, as.numeric(r2mlm_DN$R2s['f2', 'between']))
    P_Lv2_All$between_res_btw<-c(P_Lv2_All$between_res_btw, as.numeric(r2mlm_DN$R2s['m', 'between']))
    P_Lv2_All$within_fix_wthn<-c(P_Lv2_All$within_fix_wthn, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'within']))
    P_Lv2_All$within_fix_tot<-c(P_Lv2_All$within_fix_tot, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    P_Lv2_All$within_slope_var_wthn<-c(P_Lv2_All$within_slope_var_wthn, as.numeric(r2mlm_DN$Decompositions['slope variation', 'within']))
    P_Lv2_All$within_res_wthn<-c(P_Lv2_All$within_res_wthn, as.numeric(r2mlm_DN$Decompositions['sigma2', 'within']))
    P_Lv2_All$within_unmod_tot<-c(P_Lv2_All$within_unmod_tot, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                  as.numeric(r2mlm_DN$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}

hist(P_Lv2_All$between_var)
mean(P_Lv2_All$between_var)
hist(P_Lv2_All$within_var)
mean(P_Lv2_All$within_var)

hist(P_Lv2_All$between_All_tot)
hist(P_Lv2_All$between_All_btw)
hist(P_Lv2_All$within_fix_tot)


#------------------------------------------------------------------------------
P_btw_DN_uni_tot<-mean(P_Lv2_All$between_All_tot)-mean(P_Evnt_alone$between_All_tot)
P_btw_shared_tot<-mean(P_DN_alone$between_All_tot)-P_btw_DN_uni_tot
P_btw_context_uni_tot<-mean(P_Lv2_All$between_All_tot)-mean(P_DN_alone$between_All_tot)

sum(c(P_btw_DN_uni_tot, 
      P_btw_shared_tot,
      P_btw_context_uni_tot))

#The above total should equal: 
mean(P_Lv2_All$between_All_tot)

#Now getting the between subjects residual variance from final model
P_btw_res_tot<-mean(P_Lv2_All$between_var-P_Lv2_All$between_All_tot)
P_btw_res_tot+sum(c(P_btw_DN_uni_tot, 
                    P_btw_shared_tot,
                    P_btw_context_uni_tot))
#Final getting the within subjects variance: 
P_wthn_best_worst_tot<-mean(P_Lv2_All$within_fix_wthn)
P_wthn_res_tot<-mean(P_Lv2_All$within_var)-P_wthn_best_worst_tot

#Should be equal to 1
sum(c(P_wthn_best_worst_tot, 
      P_wthn_res_tot, 
      mean(P_Lv2_All$between_var)))

P_tot_btw<-mean(P_Lv2_All$between_var)
P_tot_wthn<-mean(P_Lv2_All$within_var)  

River_DF<-data.frame(N1 = c('DN',
                            'DN <--> Exposure',
                            'Exposure', 
                            'Unmodeled Between', 
                            'Momentary Events', 
                            'Unmodeled Within', 
                            'Total Between', 
                            'Total Within'), 
                     N2 = c('Total Between', 
                            'Total Between', 
                            'Total Between', 
                            'Total Between', 
                            'Total Within', 
                            'Total Within', 
                            'Total Variance', 
                            'Total Variance'), 
                     Value = c(P_btw_DN_uni_tot, 
                               P_btw_shared_tot,
                               P_btw_context_uni_tot, 
                               P_btw_res_tot, 
                               P_wthn_best_worst_tot, 
                               P_wthn_res_tot, 
                               P_tot_btw, 
                               P_tot_wthn), 
                     ID = 1:8)

River_DF$N1<-paste(River_DF$N1, 
                   '\n', 
                   paste0(round(River_DF$Value*100, 
                                digits = 2), '%'
                   )
)
River_DF$N2<-c(rep(River_DF$N1[7], 4),
               rep(River_DF$N1[8], 2), 
               rep('Total Variance', 2)
)

nodes<-data.frame(ID = c(River_DF$N1, 
                         'Total Variance'), 
                  x = c(1,1,1,1,1,1,2,2,3), 
                  y = c(0,1.5,3.25,5.25,8,10.5,2.5,7.5,5))

palette = c(paste0(brewer.pal(9, "Blues"), 90)[c(2, 4, 6, 8)], 
            paste0(brewer.pal(9, "Reds"), 90)[c(6, 8)], 
            paste0(brewer.pal(9, "Blues"), 90)[9], 
            paste0(brewer.pal(9, "Reds"), 90)[9], 
            paste0(brewer.pal(9, 'Purples'), 90)[9])

styles = lapply(nodes$y, function(n) {
  list(col = palette[n], lty = 0, textcol = "black")
})

#Not sure why color is not mapping correctly - trying to force it with this fuction... 
for(i in 1:length(palette)){
  styles[[i]]$col<-palette[i]
}

names(styles) = nodes$ID

riv<-makeRiver(nodes = nodes, 
               edges =  River_DF, 
               node_styles = styles)

#Creating Dataset to add in text - Will have to manipulate the y-values to get text to line up 

png(paste0(study2.graphics, '/S2_PosMood_river.png'), 
    units = 'in', 
    res = 1200, 
    height = 10, 
    width = 10)
riverplot(riv, 
          nodewidth = 3, 
          plot_area = .95)
title(ylab = 'Study 2 - Riverplot of Total Variance Decomposition Estimated from Positive Mood Models')
dev.off()

#------------------------------------------------------------------------------
#Extracting variance from cross-level interactions

#Cross-Level Model 
within_cov<-c(24, 25, 42, 43)           #Columns with group-mean centered predictors & cross-level interactions
between_cov<-c(22, 32, 33)               #Columns with between-subject predictors
#Note I have added a series of cross level interactions
#Variables for these cross-level interactions are calculated in imputed data sets
random_cov<-c(24, 25)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Neg_cross, pars = 'b_Intercept') 
DN = posterior_samples(Neg_cross, pars = 'b_DN_comb')
M.Worst = posterior_samples(Neg_cross, pars = 'b_prop_NegEvnt')
M.Best = posterior_samples(Neg_cross, pars = 'b_prop_PosEvnt')
Worst = posterior_samples(Neg_cross, pars = 'b_c.NegEvnt') 
Best = posterior_samples(Neg_cross, pars = 'b_c.PosEvnt')
DNxWorst = posterior_samples(Neg_cross, pars = 'b_c.NegEvnt:DN_comb')
DNxBest = posterior_samples(Neg_cross, pars = 'b_c.PosEvnt:DN_comb')

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Neg_cross, pars = 'sd_ID__Intercept')^2
Worst_var = posterior_samples(Neg_cross, pars = 'sd_ID__c.NegEvnt')^2 
Best_var = posterior_samples(Neg_cross, pars = 'sd_ID__c.PosEvnt')^2

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Worst = posterior_samples(Neg_cross, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_cross, pars = 'sd_ID__c.NegEvnt')*
  posterior_samples(Neg_cross, pars = 'cor_ID__Intercept__c.NegEvnt')

cov_Int_Best = posterior_samples(Neg_cross, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_cross, pars = 'sd_ID__c.PosEvnt')*
  posterior_samples(Neg_cross, pars = 'cor_ID__Intercept__c.PosEvnt') 

cov_Worst_Best = posterior_samples(Neg_cross, pars = 'sd_ID__c.NegEvnt')*
  posterior_samples(Neg_cross, pars = 'sd_ID__c.PosEvnt')*
  posterior_samples(Neg_cross, pars = 'cor_ID__c.NegEvnt__c.PosEvnt')

#Getting level 1 error variance
sigma<-log(1+posterior_samples(Neg_cross, pars = 'sigma')/beta_00)

post_samples<-data.frame(Intercept, 
                         DN,
                         M.Worst,
                         M.Best,
                         Worst,
                         DNxWorst,
                         Best, 
                         DNxBest,
                         Int_var, 
                         Worst_var, 
                         Best_var,
                         cov_Int_Worst, 
                         cov_Int_Best, 
                         cov_Worst_Best, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'DN',
                          'M.Worst',
                          'M.Best',
                          'Worst',
                          'DNxWorst',
                          'Best', 
                          'DNxBest',
                          'Int_var', 
                          'Worst_var', 
                          'Best_var',
                          'cov_Int_Worst', 
                          'cov_Int_Best', 
                          'cov_Worst_Best', 
                          'sigma')
#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets
Lv2_cross<-list(between_var=vector(), 
                within_var=vector(), 
                between_All_tot=vector(),
                between_All_btw=vector(),
                between_res_btw=vector(),
                within_fix_wthn=vector(),
                within_fix_tot=vector(),
                within_slope_var_wthn=vector(),
                within_res_wthn=vector(), 
                within_unmod_tot=vector())

for(i in 1:length(dat.study2.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  D<-sample(sampling_list[[i]], size = 1000, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    #d<-1
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Worst', 'Best', 'DNxWorst', 'DNxBest')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'DN', 'M.Worst', 'M.Best')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Worst', 'cov_Int_Best')], 
           post_samples[D[d], c('cov_Int_Worst', 'Worst_var', 'cov_Worst_Best')], 
           post_samples[D[d], c('cov_Int_Best', 'cov_Worst_Best', 'Best_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=3)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_All<-r2MLM(data=dat.study2.list[[i]], 
                     within_covs = within_cov, 
                     between_covs = between_cov, 
                     random_covs = random_cov, 
                     gamma_w = Gamma_w, 
                     gamma_b = Gamma_b, 
                     Tau = tau, 
                     sigma2 = Sigma2, 
                     has_intercept = TRUE, 
                     clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    Lv2_cross$between_var<-c(Lv2_cross$between_var, as.numeric(r2mlm_DN$Decompositions['fixed, between', 'total'])+
                             as.numeric(r2mlm_DN$Decompositions['mean variation', 'total']))
    Lv2_cross$within_var<-c(Lv2_cross$within_var, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                            as.numeric(r2mlm_DN$Decompositions['slope variation', 'total'])+
                            as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    Lv2_cross$between_cross_tot<-c(Lv2_cross$between_cross_tot, as.numeric(r2mlm_DN$R2s['f2', 'total']))
    Lv2_cross$between_cross_btw<-c(Lv2_cross$between_cross_btw, as.numeric(r2mlm_DN$R2s['f2', 'between']))
    Lv2_cross$between_res_btw<-c(Lv2_cross$between_res_btw, as.numeric(r2mlm_DN$R2s['m', 'between']))
    Lv2_cross$within_fix_wthn<-c(Lv2_cross$within_fix_wthn, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'within']))
    Lv2_cross$within_fix_tot<-c(Lv2_cross$within_fix_tot, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    Lv2_cross$within_slope_var_wthn<-c(Lv2_cross$within_slope_var_wthn, as.numeric(r2mlm_DN$Decompositions['slope variation', 'within']))
    Lv2_cross$within_res_wthn<-c(Lv2_cross$within_res_wthn, as.numeric(r2mlm_DN$Decompositions['sigma2', 'within']))
    Lv2_cross$within_unmod_tot<-c(Lv2_cross$within_unmod_tot, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                  as.numeric(r2mlm_DN$Decompositions['slope variation', 'total']))
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}

#------------------------------------------------------------------------------
#Reactivity riverplot:
#Obtaining DN-specific plot
tot_var_react<-mean(Lv2_cross$within_fix_tot) - mean(Lv2_All$within_fix_tot) 
DN_total <- btw_DN_uni_tot + btw_shared_tot + tot_var_react 
DN_unique_DN<-btw_DN_uni_tot/DN_total
DN_indirect_DN<-btw_shared_tot/DN_total
DN_react_DN<-tot_var_react/DN_total

sum(c(DN_unique_DN, 
      DN_indirect_DN, 
      DN_react_DN))

River_DF<-data.frame(N1 = c('DN Reactivity',
                            'DN',
                            'DN <--> Exposure'
), 
N2 = c('Total DN Effect', 
       'Total DN Effect', 
       'Total DN Effect'), 
Value = c(DN_react_DN,
          DN_unique_DN, 
          DN_indirect_DN), 
ID = 1:3)

River_DF$N1<-paste(River_DF$N1, 
                   '\n', 
                   paste0(round(River_DF$Value*100, 
                                digits = 2), '%'
                   )
)

nodes<-data.frame(ID = c(River_DF$N1, 
                         'Total DN Effect'), 
                  x = c(1,1,1,2), 
                  y = c(0,1,2,1))

palette = c(paste0(brewer.pal(9, "Reds"),90)[7], 
            paste0(brewer.pal(9, "Blues"),90)[c(7,9)], 
            paste0(brewer.pal(9, "Purples"),90)[9])

styles = lapply(nodes$y, function(n) {
  list(col = palette[n], lty = 0, textcol = "black")
})

#Not sure why color is not mapping correctly - trying to force it with this fuction... 
for(i in 1:length(palette)){
  styles[[i]]$col<-palette[i]
}

names(styles) = nodes$ID

riv<-makeRiver(nodes = nodes, 
               edges =  River_DF, 
               node_styles = styles)

png(paste0(study2.graphics, '/S2_NegMood_DN_River.png'), 
    units = 'in', 
    res = 900, 
    height = 10, 
    width = 10)
riverplot(riv, 
          nodewidth = 3, 
          plot_area = .95)
title(ylab = 'Study 2 - Riverplot of DN Effect Decomposition Derived from Negative Mood Model')
dev.off()
###############################################################################
#Positive Mood River Plot - DN effect decomp
###############################################################################
#Cross-Level Model 
within_cov<-c(24, 25, 42, 43)           #Columns with group-mean centered predictors & cross-level interactions
between_cov<-c(22, 32, 33)               #Columns with between-subject predictors
#Note I have added a series of cross level interactions
#Variables for these cross-level interactions are calculated in imputed data sets
random_cov<-c(24, 25)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Pos_cross, pars = 'b_Intercept') 
DN = posterior_samples(Pos_cross, pars = 'b_DN_comb')
M.Worst = posterior_samples(Pos_cross, pars = 'b_prop_NegEvnt')
M.Best = posterior_samples(Pos_cross, pars = 'b_prop_PosEvnt')
Worst = posterior_samples(Pos_cross, pars = 'b_c.NegEvnt') 
Best = posterior_samples(Pos_cross, pars = 'b_c.PosEvnt')
DNxWorst = posterior_samples(Pos_cross, pars = 'b_c.NegEvnt:DN_comb')
DNxBest = posterior_samples(Pos_cross, pars = 'b_c.PosEvnt:DN_comb')

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Pos_cross, pars = 'sd_ID__Intercept')^2
Worst_var = posterior_samples(Pos_cross, pars = 'sd_ID__c.NegEvnt')^2 
Best_var = posterior_samples(Pos_cross, pars = 'sd_ID__c.PosEvnt')^2

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Worst = posterior_samples(Pos_cross, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_cross, pars = 'sd_ID__c.NegEvnt')*
  posterior_samples(Pos_cross, pars = 'cor_ID__Intercept__c.NegEvnt')

cov_Int_Best = posterior_samples(Pos_cross, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_cross, pars = 'sd_ID__c.PosEvnt')*
  posterior_samples(Pos_cross, pars = 'cor_ID__Intercept__c.PosEvnt') 

cov_Worst_Best = posterior_samples(Pos_cross, pars = 'sd_ID__c.NegEvnt')*
  posterior_samples(Pos_cross, pars = 'sd_ID__c.PosEvnt')*
  posterior_samples(Pos_cross, pars = 'cor_ID__c.NegEvnt__c.PosEvnt')

#Getting level 1 error variance
sigma<-posterior_samples(Pos_cross, pars = 'sigma')

post_samples<-data.frame(Intercept, 
                         DN,
                         M.Worst,
                         M.Best,
                         Worst,
                         DNxWorst,
                         Best, 
                         DNxBest,
                         Int_var, 
                         Worst_var, 
                         Best_var,
                         cov_Int_Worst, 
                         cov_Int_Best, 
                         cov_Worst_Best, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'DN',
                          'M.Worst',
                          'M.Best',
                          'Worst',
                          'DNxWorst',
                          'Best', 
                          'DNxBest',
                          'Int_var', 
                          'Worst_var', 
                          'Best_var',
                          'cov_Int_Worst', 
                          'cov_Int_Best', 
                          'cov_Worst_Best', 
                          'sigma')
#Aggregate across imputed datasets
#Currently going to take 1000 draws from posterior distributions
#Will then apply across all 20 data sets
P_Lv2_cross<-list(between_var=vector(), 
                within_var=vector(), 
                between_All_tot=vector(),
                between_All_btw=vector(),
                between_res_btw=vector(),
                within_fix_wthn=vector(),
                within_fix_tot=vector(),
                within_slope_var_wthn=vector(),
                within_res_wthn=vector(), 
                within_unmod_tot=vector())

for(i in 1:length(dat.study2.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  D<-sample(sampling_list[[i]], size = 1000, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    #d<-1
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Worst', 'Best', 'DNxWorst', 'DNxBest')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'DN', 'M.Worst', 'M.Best')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Worst', 'cov_Int_Best')], 
           post_samples[D[d], c('cov_Int_Worst', 'Worst_var', 'cov_Worst_Best')], 
           post_samples[D[d], c('cov_Int_Best', 'cov_Worst_Best', 'Best_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=3)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_All<-r2MLM(data=dat.study2.list[[i]], 
                     within_covs = within_cov, 
                     between_covs = between_cov, 
                     random_covs = random_cov, 
                     gamma_w = Gamma_w, 
                     gamma_b = Gamma_b, 
                     Tau = tau, 
                     sigma2 = Sigma2, 
                     has_intercept = TRUE, 
                     clustermeancentered = TRUE)
    
    #Extracting relevant values across imputed datasets
    P_Lv2_cross$between_var<-c(P_Lv2_cross$between_var, as.numeric(r2mlm_DN$Decompositions['fixed, between', 'total'])+
                               as.numeric(r2mlm_DN$Decompositions['mean variation', 'total']))
    P_Lv2_cross$within_var<-c(P_Lv2_cross$within_var, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                              as.numeric(r2mlm_DN$Decompositions['slope variation', 'total'])+
                              as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    P_Lv2_cross$between_cross_tot<-c(P_Lv2_cross$between_cross_tot, as.numeric(r2mlm_DN$R2s['f2', 'total']))
    P_Lv2_cross$between_cross_btw<-c(P_Lv2_cross$between_cross_btw, as.numeric(r2mlm_DN$R2s['f2', 'between']))
    P_Lv2_cross$between_res_btw<-c(P_Lv2_cross$between_res_btw, as.numeric(r2mlm_DN$R2s['m', 'between']))
    P_Lv2_cross$within_fix_wthn<-c(P_Lv2_cross$within_fix_wthn, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'within']))
    P_Lv2_cross$within_fix_tot<-c(P_Lv2_cross$within_fix_tot, as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total']))
    P_Lv2_cross$within_slope_var_wthn<-c(P_Lv2_cross$within_slope_var_wthn, as.numeric(r2mlm_DN$Decompositions['slope variation', 'within']))
    P_Lv2_cross$within_res_wthn<-c(P_Lv2_cross$within_res_wthn, as.numeric(r2mlm_DN$Decompositions['sigma2', 'within']))
    P_Lv2_cross$within_unmod_tot<-c(P_Lv2_cross$within_unmod_tot, as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
                                    as.numeric(r2mlm_DN$Decompositions['slope variation', 'total']))
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}

#------------------------------------------------------------------------------
#Reactivity riverplot:
#Obtaining DN-specific plot
P_tot_var_react<-mean(P_Lv2_cross$within_fix_tot) - mean(P_Lv2_All$within_fix_tot) 
P_DN_total <- P_btw_DN_uni_tot + P_btw_shared_tot + P_tot_var_react 
P_DN_unique_DN<-P_btw_DN_uni_tot/P_DN_total
P_DN_indirect_DN<-P_btw_shared_tot/P_DN_total
P_DN_react_DN<-P_tot_var_react/P_DN_total

sum(c(P_DN_unique_DN, 
      P_DN_indirect_DN, 
      P_DN_react_DN))

River_DF<-data.frame(N1 = c('DN Reactivity',
                            'DN',
                            'DN <--> Exposure'
), 
N2 = c('Total DN Effect', 
       'Total DN Effect', 
       'Total DN Effect'), 
Value = c(P_DN_react_DN,
          P_DN_unique_DN, 
          P_DN_indirect_DN), 
ID = 1:3)

River_DF$N1<-paste(River_DF$N1, 
                   '\n', 
                   paste0(round(River_DF$Value*100, 
                                digits = 2), '%'
                   )
)

nodes<-data.frame(ID = c(River_DF$N1, 
                         'Total DN Effect'), 
                  x = c(1,1,1,2), 
                  y = c(0,1,2,1))

palette = c(paste0(brewer.pal(9, "Reds"),90)[7], 
            paste0(brewer.pal(9, "Blues"),90)[c(7,9)], 
            paste0(brewer.pal(9, "Purples"),90)[9])

styles = lapply(nodes$y, function(n) {
  list(col = palette[n], lty = 0, textcol = "black")
})

#Not sure why color is not mapping correctly - trying to force it with this fuction... 
for(i in 1:length(palette)){
  styles[[i]]$col<-palette[i]
}

names(styles) = nodes$ID

riv<-makeRiver(nodes = nodes, 
               edges =  River_DF, 
               node_styles = styles)

png(paste0(study2.graphics, '/S2_PosMood_DN_River.png'), 
    units = 'in', 
    res = 900, 
    height = 10, 
    width = 10)
riverplot(riv, 
          nodewidth = 3, 
          plot_area = .95)
title(ylab = 'Study 2 - Riverplot of DN Effect Decomposition Derived from Positive Mood Model')
dev.off()

###############################################################################
#Simple Exposure Models
#DN --> probability of positive event 
#DN -- probability of negative event
###############################################################################
NegEvnt_DN<-brms::brm_multiple(NegEvnt_char~1+DN_comb+
                                 (1|ID), 
                               data = dat.study2.list, 
                               family = 'bernoulli',
                               prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                         set_prior('normal(0,5)', class='Intercept'), 
                                         set_prior('normal(0,5)', class='b')),
                               warmup = 2000, 
                               iter = 3000, 
                               chains = 3,
                               control = list(adapt_delta=.99, 
                                              max_treedepth=15), 
                               save_model = paste0(stan.code, '/S2_Evnt_DN.stan'))

PosEvnt_DN<-brms::brm_multiple(PosEvnt_char~1+DN_comb+
                                 (1|ID), 
                               data = dat.study2.list, 
                               family = 'bernoulli',
                               prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                         set_prior('normal(0,5)', class='Intercept'), 
                                         set_prior('normal(0,5)', class='b')),
                               warmup = 2000, 
                               iter = 3000, 
                               chains = 3,
                               control = list(adapt_delta=.99, 
                                              max_treedepth=15), 
                               save_model = paste0(stan.code, '/S2_Evnt_DN.stan'))




###############################################################################
#Extracting Model Infomartion to input manually into manuscript:
###############################################################################

#HDIs for model parameters 
hist(DN_alone$between_All_tot)
round(quantile(DN_alone$between_All_tot, c(.025, .975)), digits = 4)*100
round(quantile(P_DN_alone$between_All_tot, c(.025, .975)), digits = 4)*100
round(quantile(Evnt_alone$between_All_tot, c(.025, .975)), digits = 4)*100
round(quantile(P_Evnt_alone$between_All_tot, c(.025, .975)), digits = 4)*100


###############################################################################
#Secondary Mood Variance Decompositions and Analysis: 
###############################################################################
#--
#Anxiety Models
mean(dat.study2$ANX, na.rm=TRUE)
log(mean(dat.study2$ANX, na.rm=TRUE))
hist(dat.study2$ANX)

Anx_lv2_DN<-brms::brm_multiple(T.ANX~1+c.NegEvnt+c.PosEvnt+
                                 DN_comb+
                                 (1+c.NegEvnt+c.PosEvnt|ID), 
                               data = dat.study2.list, 
                               family = 'lognormal',
                               prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                         set_prior('normal(0.3014526,2)', class='Intercept'), 
                                         set_prior('normal(0,5)', class='b'), 
                                         set_prior('lkj(2)', class='cor')),
                               warmup = 2000, 
                               iter = 3000, 
                               chains = 3,
                               control = list(adapt_delta=.99, 
                                              max_treedepth=15), 
                               save_model = paste0(stan.code, '/S2_lv2_DN.stan'))
gc()

Anx_lv2_Evnt<-brms::brm_multiple(T.ANX~1+c.NegEvnt+c.PosEvnt+
                                   prop_NegEvnt+prop_PosEvnt+
                                   (1+c.NegEvnt+c.PosEvnt|ID), 
                                 data = dat.study2.list, 
                                 family = 'lognormal',
                                 prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                           set_prior('normal(0.3014526,2)', class='Intercept'), 
                                           set_prior('normal(0,5)', class='b'), 
                                           set_prior('lkj(2)', class='cor')),
                                 warmup = 2000, 
                                 iter = 3000, 
                                 chains = 3,
                                 control = list(adapt_delta=.99, 
                                                max_treedepth=15), 
                                 save_model = paste0(stan.code, '/S2_lv2_Evnt.stan'))
gc()

Anx_lv2_All<-brms::brm_multiple(T.ANX~1+c.NegEvnt+c.PosEvnt+
                                  prop_NegEvnt+prop_PosEvnt+DN_comb+
                                  (1+c.NegEvnt+c.PosEvnt|ID), 
                                data = dat.study2.list, 
                                family = 'lognormal',
                                prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                          set_prior('normal(0.3014526,2)', class='Intercept'), 
                                          set_prior('normal(0,5)', class='b'), 
                                          set_prior('lkj(2)', class='cor')),
                                warmup = 2000, 
                                iter = 3000, 
                                chains = 3,
                                control = list(adapt_delta=.99, 
                                               max_treedepth=15), 
                                save_model = paste0(stan.code, '/S2_lv2_All.stan'))
gc()

Anx_cross<-brms::brm_multiple(T.ANX~1+c.NegEvnt+c.PosEvnt+
                                prop_NegEvnt+prop_PosEvnt+DN_comb+
                                c.PosEvnt:DN_comb+c.NegEvnt:DN_comb+
                                (1+c.NegEvnt+c.PosEvnt|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                        set_prior('normal(0.3014526,2)', class='Intercept'), 
                                        set_prior('normal(0,5)', class='b'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 2000, 
                              iter = 3000, 
                              chains = 3,
                              control = list(adapt_delta=.99, 
                                             max_treedepth=15), 
                              save_model = paste0(stan.code, '/S2_cross.stan'))
gc()

#-- 
#Depression Models
mean(dat.study2$DEP, na.rm=TRUE)
log(mean(dat.study2$DEP, na.rm=TRUE))
hist(dat.study2$DEP)

Dep_lv2_DN<-brms::brm_multiple(T.NEG~1+c.NegEvnt+c.PosEvnt+
                                 DN_comb+
                                 (1+c.NegEvnt+c.PosEvnt|ID), 
                               data = dat.study2.list, 
                               family = 'lognormal',
                               prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                         set_prior('normal(0.3014526,2)', class='Intercept'), 
                                         set_prior('normal(0,5)', class='b'), 
                                         set_prior('lkj(2)', class='cor')),
                               warmup = 2000, 
                               iter = 3000, 
                               chains = 3,
                               control = list(adapt_delta=.99, 
                                              max_treedepth=15), 
                               save_model = paste0(stan.code, '/S2_lv2_DN.stan'))
gc()

Dep_lv2_Evnt<-brms::brm_multiple(T.NEG~1+c.NegEvnt+c.PosEvnt+
                                   prop_NegEvnt+prop_PosEvnt+
                                   (1+c.NegEvnt+c.PosEvnt|ID), 
                                 data = dat.study2.list, 
                                 family = 'lognormal',
                                 prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                           set_prior('normal(0.3014526,2)', class='Intercept'), 
                                           set_prior('normal(0,5)', class='b'), 
                                           set_prior('lkj(2)', class='cor')),
                                 warmup = 2000, 
                                 iter = 3000, 
                                 chains = 3,
                                 control = list(adapt_delta=.99, 
                                                max_treedepth=15), 
                                 save_model = paste0(stan.code, '/S2_lv2_Evnt.stan'))
gc()

Dep_lv2_All<-brms::brm_multiple(T.NEG~1+c.NegEvnt+c.PosEvnt+
                                  prop_NegEvnt+prop_PosEvnt+DN_comb+
                                  (1+c.NegEvnt+c.PosEvnt|ID), 
                                data = dat.study2.list, 
                                family = 'lognormal',
                                prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                          set_prior('normal(0.3014526,2)', class='Intercept'), 
                                          set_prior('normal(0,5)', class='b'), 
                                          set_prior('lkj(2)', class='cor')),
                                warmup = 2000, 
                                iter = 3000, 
                                chains = 3,
                                control = list(adapt_delta=.99, 
                                               max_treedepth=15), 
                                save_model = paste0(stan.code, '/S2_lv2_All.stan'))
gc()

Dep_cross<-brms::brm_multiple(T.NEG~1+c.NegEvnt+c.PosEvnt+
                                prop_NegEvnt+prop_PosEvnt+DN_comb+
                                c.PosEvnt:DN_comb+c.NegEvnt:DN_comb+
                                (1+c.NegEvnt+c.PosEvnt|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                        set_prior('normal(0.3014526,2)', class='Intercept'), 
                                        set_prior('normal(0,5)', class='b'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 2000, 
                              iter = 3000, 
                              chains = 3,
                              control = list(adapt_delta=.99, 
                                             max_treedepth=15), 
                              save_model = paste0(stan.code, '/S2_cross.stan'))
gc()
pp_check(Dep_cross)

#--
#Cheerful Models
mean(dat.study2$JOY, na.rm=TRUE)
log(mean(dat.study2$JOY, na.rm=TRUE))
hist(dat.study2$JOY)

Joy_lv2_DN<-brms::brm_multiple(T.JOY~1+c.NegEvnt+c.PosEvnt+
                                 DN_comb+
                                 (1+c.NegEvnt+c.PosEvnt|ID), 
                               data = dat.study2.list, 
                               family = 'lognormal',
                               prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                         set_prior('normal(0.3014526,2)', class='Intercept'), 
                                         set_prior('normal(0,5)', class='b'), 
                                         set_prior('lkj(2)', class='cor')),
                               warmup = 2000, 
                               iter = 3000, 
                               chains = 3,
                               control = list(adapt_delta=.99, 
                                              max_treedepth=15), 
                               save_model = paste0(stan.code, '/S2_lv2_DN.stan'))
gc()

Joy_lv2_Evnt<-brms::brm_multiple(T.JOY~1+c.NegEvnt+c.PosEvnt+
                                   prop_NegEvnt+prop_PosEvnt+
                                   (1+c.NegEvnt+c.PosEvnt|ID), 
                                 data = dat.study2.list, 
                                 family = 'lognormal',
                                 prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                           set_prior('normal(0.3014526,2)', class='Intercept'), 
                                           set_prior('normal(0,5)', class='b'), 
                                           set_prior('lkj(2)', class='cor')),
                                 warmup = 2000, 
                                 iter = 3000, 
                                 chains = 3,
                                 control = list(adapt_delta=.99, 
                                                max_treedepth=15), 
                                 save_model = paste0(stan.code, '/S2_lv2_Evnt.stan'))
gc()

Joy_lv2_All<-brms::brm_multiple(T.JOY~1+c.NegEvnt+c.PosEvnt+
                                  prop_NegEvnt+prop_PosEvnt+DN_comb+
                                  (1+c.NegEvnt+c.PosEvnt|ID), 
                                data = dat.study2.list, 
                                family = 'lognormal',
                                prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                          set_prior('normal(0.3014526,2)', class='Intercept'), 
                                          set_prior('normal(0,5)', class='b'), 
                                          set_prior('lkj(2)', class='cor')),
                                warmup = 2000, 
                                iter = 3000, 
                                chains = 3,
                                control = list(adapt_delta=.99, 
                                               max_treedepth=15), 
                                save_model = paste0(stan.code, '/S2_lv2_All.stan'))
gc()

Joy_cross<-brms::brm_multiple(T.JOY~1+c.NegEvnt+c.PosEvnt+
                                prop_NegEvnt+prop_PosEvnt+DN_comb+
                                c.PosEvnt:DN_comb+c.NegEvnt:DN_comb+
                                (1+c.NegEvnt+c.PosEvnt|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                        set_prior('normal(0.3014526,2)', class='Intercept'), 
                                        set_prior('normal(0,5)', class='b'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 2000, 
                              iter = 3000, 
                              chains = 3,
                              control = list(adapt_delta=.99, 
                                             max_treedepth=15), 
                              save_model = paste0(stan.code, '/S2_cross.stan'))
gc()

#-- 
#Calm Models
#Anxiety Models
mean(dat.study2$CALM, na.rm=TRUE)
log(mean(dat.study2$CALM, na.rm=TRUE))
hist(dat.study2$CALM)

Calm_lv2_DN<-brms::brm_multiple(T.CALM~1+c.NegEvnt+c.PosEvnt+
                                 DN_comb+
                                 (1+c.NegEvnt+c.PosEvnt|ID), 
                               data = dat.study2.list, 
                               family = 'lognormal',
                               prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                         set_prior('normal(0.3014526,2)', class='Intercept'), 
                                         set_prior('normal(0,5)', class='b'), 
                                         set_prior('lkj(2)', class='cor')),
                               warmup = 2000, 
                               iter = 3000, 
                               chains = 3,
                               control = list(adapt_delta=.99, 
                                              max_treedepth=15), 
                               save_model = paste0(stan.code, '/S2_lv2_DN.stan'))
gc()

Calm_lv2_Evnt<-brms::brm_multiple(T.CALM~1+c.NegEvnt+c.PosEvnt+
                                   prop_NegEvnt+prop_PosEvnt+
                                   (1+c.NegEvnt+c.PosEvnt|ID), 
                                 data = dat.study2.list, 
                                 family = 'lognormal',
                                 prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                           set_prior('normal(0.3014526,2)', class='Intercept'), 
                                           set_prior('normal(0,5)', class='b'), 
                                           set_prior('lkj(2)', class='cor')),
                                 warmup = 2000, 
                                 iter = 3000, 
                                 chains = 3,
                                 control = list(adapt_delta=.99, 
                                                max_treedepth=15), 
                                 save_model = paste0(stan.code, '/S2_lv2_Evnt.stan'))
gc()

Calm_lv2_All<-brms::brm_multiple(T.CALM~1+c.NegEvnt+c.PosEvnt+
                                  prop_NegEvnt+prop_PosEvnt+DN_comb+
                                  (1+c.NegEvnt+c.PosEvnt|ID), 
                                data = dat.study2.list, 
                                family = 'lognormal',
                                prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                          set_prior('normal(0.3014526,2)', class='Intercept'), 
                                          set_prior('normal(0,5)', class='b'), 
                                          set_prior('lkj(2)', class='cor')),
                                warmup = 2000, 
                                iter = 3000, 
                                chains = 3,
                                control = list(adapt_delta=.99, 
                                               max_treedepth=15), 
                                save_model = paste0(stan.code, '/S2_lv2_All.stan'))
gc()

Calm_cross<-brms::brm_multiple(T.CALM~1+c.NegEvnt+c.PosEvnt+
                                prop_NegEvnt+prop_PosEvnt+DN_comb+
                                c.PosEvnt:DN_comb+c.NegEvnt:DN_comb+
                                (1+c.NegEvnt+c.PosEvnt|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                        set_prior('normal(0.3014526,2)', class='Intercept'), 
                                        set_prior('normal(0,5)', class='b'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 2000, 
                              iter = 3000, 
                              chains = 3,
                              control = list(adapt_delta=.99, 
                                             max_treedepth=15), 
                              save_model = paste0(stan.code, '/S2_cross.stan'))
gc()

###############################################################################
#Quick Model Variance Extraction
###############################################################################

Anx_lv2_DN<-brms::brm_multiple(T.ANX~1+c.NegEvnt+c.PosEvnt+
                                 DN_comb+
                                 (1+c.NegEvnt+c.PosEvnt|ID), 
                               data = dat.study2.list, 
                               family = 'lognormal',
                               prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                         set_prior('normal(0.3014526,2)', class='Intercept'), 
                                         set_prior('normal(0,5)', class='b'), 
                                         set_prior('lkj(2)', class='cor')),
                               warmup = 2000, 
                               iter = 3000, 
                               chains = 3,
                               control = list(adapt_delta=.99, 
                                              max_treedepth=15), 
                               save_model = paste0(stan.code, '/S2_lv2_DN.stan'))
gc()

Anx_lv2_Evnt<-brms::brm_multiple(T.ANX~1+c.NegEvnt+c.PosEvnt+
                                   prop_NegEvnt+prop_PosEvnt+
                                   (1+c.NegEvnt+c.PosEvnt|ID), 
                                 data = dat.study2.list, 
                                 family = 'lognormal',
                                 prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                           set_prior('normal(0.3014526,2)', class='Intercept'), 
                                           set_prior('normal(0,5)', class='b'), 
                                           set_prior('lkj(2)', class='cor')),
                                 warmup = 2000, 
                                 iter = 3000, 
                                 chains = 3,
                                 control = list(adapt_delta=.99, 
                                                max_treedepth=15), 
                                 save_model = paste0(stan.code, '/S2_lv2_Evnt.stan'))
gc()

Anx_lv2_All<-brms::brm_multiple(T.ANX~1+c.NegEvnt+c.PosEvnt+
                                  prop_NegEvnt+prop_PosEvnt+DN_comb+
                                  (1+c.NegEvnt+c.PosEvnt|ID), 
                                data = dat.study2.list, 
                                family = 'lognormal',
                                prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                          set_prior('normal(0.3014526,2)', class='Intercept'), 
                                          set_prior('normal(0,5)', class='b'), 
                                          set_prior('lkj(2)', class='cor')),
                                warmup = 2000, 
                                iter = 3000, 
                                chains = 3,
                                control = list(adapt_delta=.99, 
                                               max_treedepth=15), 
                                save_model = paste0(stan.code, '/S2_lv2_All.stan'))
gc()

Anx_cross<-brms::brm_multiple(T.ANX~1+c.NegEvnt+c.PosEvnt+
                                prop_NegEvnt+prop_PosEvnt+DN_comb+
                                c.PosEvnt:DN_comb+c.NegEvnt:DN_comb+
                                (1+c.NegEvnt+c.PosEvnt|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                        set_prior('normal(0.3014526,2)', class='Intercept'), 
                                        set_prior('normal(0,5)', class='b'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 2000, 
                              iter = 3000, 
                              chains = 3,
                              control = list(adapt_delta=.99, 
                                             max_treedepth=15), 
                              save_model = paste0(stan.code, '/S2_cross.stan'))
gc()

#-- 
#Depression Models
mean(dat.study2$DEP, na.rm=TRUE)
log(mean(dat.study2$DEP, na.rm=TRUE))
hist(dat.study2$DEP)

Dep_lv2_DN<-brms::brm_multiple(T.NEG~1+c.NegEvnt+c.PosEvnt+
                                 DN_comb+
                                 (1+c.NegEvnt+c.PosEvnt|ID), 
                               data = dat.study2.list, 
                               family = 'lognormal',
                               prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                         set_prior('normal(0.3014526,2)', class='Intercept'), 
                                         set_prior('normal(0,5)', class='b'), 
                                         set_prior('lkj(2)', class='cor')),
                               warmup = 2000, 
                               iter = 3000, 
                               chains = 3,
                               control = list(adapt_delta=.99, 
                                              max_treedepth=15), 
                               save_model = paste0(stan.code, '/S2_lv2_DN.stan'))
gc()

Dep_lv2_Evnt<-brms::brm_multiple(T.NEG~1+c.NegEvnt+c.PosEvnt+
                                   prop_NegEvnt+prop_PosEvnt+
                                   (1+c.NegEvnt+c.PosEvnt|ID), 
                                 data = dat.study2.list, 
                                 family = 'lognormal',
                                 prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                           set_prior('normal(0.3014526,2)', class='Intercept'), 
                                           set_prior('normal(0,5)', class='b'), 
                                           set_prior('lkj(2)', class='cor')),
                                 warmup = 2000, 
                                 iter = 3000, 
                                 chains = 3,
                                 control = list(adapt_delta=.99, 
                                                max_treedepth=15), 
                                 save_model = paste0(stan.code, '/S2_lv2_Evnt.stan'))
gc()

Dep_lv2_All<-brms::brm_multiple(T.NEG~1+c.NegEvnt+c.PosEvnt+
                                  prop_NegEvnt+prop_PosEvnt+DN_comb+
                                  (1+c.NegEvnt+c.PosEvnt|ID), 
                                data = dat.study2.list, 
                                family = 'lognormal',
                                prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                          set_prior('normal(0.3014526,2)', class='Intercept'), 
                                          set_prior('normal(0,5)', class='b'), 
                                          set_prior('lkj(2)', class='cor')),
                                warmup = 2000, 
                                iter = 3000, 
                                chains = 3,
                                control = list(adapt_delta=.99, 
                                               max_treedepth=15), 
                                save_model = paste0(stan.code, '/S2_lv2_All.stan'))
gc()

Dep_cross<-brms::brm_multiple(T.NEG~1+c.NegEvnt+c.PosEvnt+
                                prop_NegEvnt+prop_PosEvnt+DN_comb+
                                c.PosEvnt:DN_comb+c.NegEvnt:DN_comb+
                                (1+c.NegEvnt+c.PosEvnt|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                        set_prior('normal(0.3014526,2)', class='Intercept'), 
                                        set_prior('normal(0,5)', class='b'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 2000, 
                              iter = 3000, 
                              chains = 3,
                              control = list(adapt_delta=.99, 
                                             max_treedepth=15), 
                              save_model = paste0(stan.code, '/S2_cross.stan'))
gc()
pp_check(Dep_cross)

#--
#Cheerful Models
mean(dat.study2$JOY, na.rm=TRUE)
log(mean(dat.study2$JOY, na.rm=TRUE))
hist(dat.study2$JOY)

Joy_lv2_DN<-brms::brm_multiple(T.JOY~1+c.NegEvnt+c.PosEvnt+
                                 DN_comb+
                                 (1+c.NegEvnt+c.PosEvnt|ID), 
                               data = dat.study2.list, 
                               family = 'lognormal',
                               prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                         set_prior('normal(0.3014526,2)', class='Intercept'), 
                                         set_prior('normal(0,5)', class='b'), 
                                         set_prior('lkj(2)', class='cor')),
                               warmup = 2000, 
                               iter = 3000, 
                               chains = 3,
                               control = list(adapt_delta=.99, 
                                              max_treedepth=15), 
                               save_model = paste0(stan.code, '/S2_lv2_DN.stan'))
gc()

Joy_lv2_Evnt<-brms::brm_multiple(T.JOY~1+c.NegEvnt+c.PosEvnt+
                                   prop_NegEvnt+prop_PosEvnt+
                                   (1+c.NegEvnt+c.PosEvnt|ID), 
                                 data = dat.study2.list, 
                                 family = 'lognormal',
                                 prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                           set_prior('normal(0.3014526,2)', class='Intercept'), 
                                           set_prior('normal(0,5)', class='b'), 
                                           set_prior('lkj(2)', class='cor')),
                                 warmup = 2000, 
                                 iter = 3000, 
                                 chains = 3,
                                 control = list(adapt_delta=.99, 
                                                max_treedepth=15), 
                                 save_model = paste0(stan.code, '/S2_lv2_Evnt.stan'))
gc()

Joy_lv2_All<-brms::brm_multiple(T.JOY~1+c.NegEvnt+c.PosEvnt+
                                  prop_NegEvnt+prop_PosEvnt+DN_comb+
                                  (1+c.NegEvnt+c.PosEvnt|ID), 
                                data = dat.study2.list, 
                                family = 'lognormal',
                                prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                          set_prior('normal(0.3014526,2)', class='Intercept'), 
                                          set_prior('normal(0,5)', class='b'), 
                                          set_prior('lkj(2)', class='cor')),
                                warmup = 2000, 
                                iter = 3000, 
                                chains = 3,
                                control = list(adapt_delta=.99, 
                                               max_treedepth=15), 
                                save_model = paste0(stan.code, '/S2_lv2_All.stan'))
gc()

Joy_cross<-brms::brm_multiple(T.JOY~1+c.NegEvnt+c.PosEvnt+
                                prop_NegEvnt+prop_PosEvnt+DN_comb+
                                c.PosEvnt:DN_comb+c.NegEvnt:DN_comb+
                                (1+c.NegEvnt+c.PosEvnt|ID), 
                              data = dat.study2.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                        set_prior('normal(0.3014526,2)', class='Intercept'), 
                                        set_prior('normal(0,5)', class='b'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 2000, 
                              iter = 3000, 
                              chains = 3,
                              control = list(adapt_delta=.99, 
                                             max_treedepth=15), 
                              save_model = paste0(stan.code, '/S2_cross.stan'))
gc()

#-- 
#Calm Models
#Anxiety Models
mean(dat.study2$CALM, na.rm=TRUE)
log(mean(dat.study2$CALM, na.rm=TRUE))
hist(dat.study2$CALM)

Calm_lv2_DN<-brms::brm_multiple(T.CALM~1+c.NegEvnt+c.PosEvnt+
                                  DN_comb+
                                  (1+c.NegEvnt+c.PosEvnt|ID), 
                                data = dat.study2.list, 
                                family = 'lognormal',
                                prior = c(set_prior('student_t(3,0,10)', class='sd'), 
                                          set_prior('normal(0.3014526,2)', class='Intercept'), 
                                          set_prior('normal(0,5)', class='b'), 
                                          set_prior('lkj(2)', class='cor')),
                                warmup = 2000, 
                                iter = 3000, 
                                chains = 3,
                                control = list(adapt_delta=.99, 
                                               max_treedepth=15), 
                                save_model = paste0(stan.code, '/S2_lv2_DN.stan'))
gc()


#------------------------------------------------------------------------------
#Calm Model 1 - DN predicting Calm Mood
#------------------------------------------------------------------------------
Calm_lv2_DN<-lme4::lmer(CALM~1+c.NegEvnt+c.PosEvnt+
                          DN_comb+
                           (1+c.NegEvnt+c.PosEvnt|ID), 
                         data = dat.study2)
summary(Calm_lv2_DN)

#******************************************************************************
within_cov<-c(22,23)
between_cov<-c(20)
random_cov<-c(22,23)
Gamma_w<-lme4::fixef(Calm_lv2_DN)[2:3]
Gamma_b<-lme4::fixef(Calm_lv2_DN)[c(1,4)]
tau<-as.matrix(vcov(Calm_lv2_DN)[1:3, 1:3])
Sigma2<-sigma(Calm_lv2_DN)

#Matrix rows/columns ordered starting with tau.00 (intercept variance)
#then add columns/rows for random effects of slopes in order of within_cov

r2mlm_DN<-r2MLM(data=dat.study2, 
                 within_covs = within_cov, 
                 between_covs = between_cov, 
                 random_covs = random_cov, 
                 gamma_w = Gamma_w, 
                 gamma_b = Gamma_b, 
                 Tau = tau, 
                 sigma2 = Sigma2, 
                 has_intercept = TRUE, 
                 clustermeancentered = TRUE)
#******************************************************************************

#Extracting relevant values across imputed datasets
DN_between_var<-as.numeric(r2mlm_DN$Decompositions['fixed, between', 'total'])+
  as.numeric(r2mlm_DN$Decompositions['mean variation', 'total'])
DN_within_var<-as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
  as.numeric(r2mlm_DN$Decompositions['slope variation', 'total'])+
  as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total'])
DN_between_DN_tot<-as.numeric(r2mlm_DN$R2s['f2', 'total'])
DN_between_DN_btw<-as.numeric(r2mlm_DN$R2s['f2', 'between'])
DN_between_res_btw<-as.numeric(r2mlm_DN$R2s['m', 'between'])
DN_within_fix_wthn<-as.numeric(r2mlm_DN$Decompositions['fixed, within', 'within'])
DN_within_fix_tot<-as.numeric(r2mlm_DN$Decompositions['fixed, within', 'total'])
DN_within_slope_var_wthn<-as.numeric(r2mlm_DN$Decompositions['slope variation', 'within'])
DN_within_res_wthn<-as.numeric(r2mlm_DN$Decompositions['sigma2', 'within'])
DN_within_unmod_tot<-as.numeric(r2mlm_DN$Decompositions['sigma2', 'total'])+
  as.numeric(r2mlm_DN$Decompositions['slope variation', 'total'])

#------------------------------------------------------------------------------
#Calm Model 2 - DN & average exposure predicting calm mood. 
#------------------------------------------------------------------------------
Calm_lv2_All<-lme4::lmer(CALM~1+c.NegEvnt+c.PosEvnt+
                           avg.NegEvnt+avg.PosEvnt+DN_comb+
                           (1+c.NegEvnt+c.PosEvnt|ID), 
                         data = dat.study2)
summary(Calm_lv2_All)

#******************************************************************************
within_cov<-c(22,23)
between_cov<-c(20,24,25)
random_cov<-c(22,23)
Gamma_w<-lme4::fixef(Calm_lv2_All)[2:3]
Gamma_b<-lme4::fixef(Calm_lv2_All)[c(1,4:6)]
tau<-as.matrix(vcov(Calm_lv2_All)[1:3, 1:3])
Sigma2<-sigma(Calm_lv2_All)

#Matrix rows/columns ordered starting with tau.00 (intercept variance)
#then add columns/rows for random effects of slopes in order of within_cov

r2mlm_All<-r2MLM(data=dat.study2, 
                within_covs = within_cov, 
                between_covs = between_cov, 
                random_covs = random_cov, 
                gamma_w = Gamma_w, 
                gamma_b = Gamma_b, 
                Tau = tau, 
                sigma2 = Sigma2, 
                has_intercept = TRUE, 
                clustermeancentered = TRUE)
#******************************************************************************

#Extracting relevant values across imputed datasets
All_between_var<-as.numeric(r2mlm_All$Decompositions['fixed, between', 'total'])+
                           as.numeric(r2mlm_All$Decompositions['mean variation', 'total'])
All_within_var<-as.numeric(r2mlm_All$Decompositions['sigma2', 'total'])+
                          as.numeric(r2mlm_All$Decompositions['slope variation', 'total'])+
                          as.numeric(r2mlm_All$Decompositions['fixed, within', 'total'])
All_between_All_tot<-as.numeric(r2mlm_All$R2s['f2', 'total'])
All_between_All_btw<-as.numeric(r2mlm_All$R2s['f2', 'between'])
All_between_res_btw<-as.numeric(r2mlm_All$R2s['m', 'between'])
All_within_fix_wthn<-as.numeric(r2mlm_All$Decompositions['fixed, within', 'within'])
All_within_fix_tot<-as.numeric(r2mlm_All$Decompositions['fixed, within', 'total'])
All_within_slope_var_wthn<-as.numeric(r2mlm_All$Decompositions['slope variation', 'within'])
All_within_res_wthn<-as.numeric(r2mlm_All$Decompositions['sigma2', 'within'])
All_within_unmod_tot<-as.numeric(r2mlm_All$Decompositions['sigma2', 'total'])+
                                as.numeric(r2mlm_All$Decompositions['slope variation', 'total'])

#******************************************************************************
gc()

Calm_cross<-lme4::lmer(CALM~1+c.NegEvnt+c.PosEvnt+
                         avg.NegEvnt+avg.PosEvnt+DN_comb+
                         c.PosEvnt:DN_comb+c.NegEvnt:DN_comb+
                         (1+c.NegEvnt+c.PosEvnt|ID), 
                       data = dat.study2)

summary(Calm_cross)
gc()
