############################################################################################################
#Authors:     Barstead, M.G., DeYoung, K.D., Anderson, A. S., & Shackman, A. J.

#Title:       The moment-to-moment affective experience of dispositionally negative individuals 

#Contact:     barstead@umd.edu

#Contents:    Bayesian regression analyses - Study 1
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
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
############################################################################################################
wd<-'/media/matthew/HDD1/EMA_MS_drafts'
data.folder<-paste0(wd, '/Model_Data')
study1.out<-paste0(wd, '/Study1')
study1.graphics<-paste0(study1.out, '/Graphics')
study1.model<-paste0(study1.out, '/Model_Summaries')
stan.code<-paste0(wd, '/Stan_code')

#Loading Study 1 Data (from Emotion MS - Shackman et al. 2017)
load(paste0(data.folder, '/Emotion MS environment.RData'))

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

############################################################################################################
#Getting vector of IDs (required for group mean centering values)
IDs<-unique(dat$subid)

ID<-vector()
Worst<-vector()
Best<-vector()
c.DN<-vector()
NegAff<-vector()
PosAff<-vector()

#Selecting additional dispositional variables to consider in the imputation. 
Gender<-vector()
BFI_E<-vector()
BFI_A<-vector()
BFI_C<-vector()
BFI_O<-vector()
GenDep<-vector()

#Bringing in covariates for imputation - 
#Want to include gender, BFI_E, BFI_C, BFI_O, BFI_A, GenDep

for(i in 1:length(IDs)){
  dat.temp<-dat[dat$subid==IDs[i],]
  for(j in 1:length(dat.temp[,1])){
    ID<-c(ID, IDs[i])
    Worst<-c(Worst, dat.temp$WorstEvent_Neg[j])
    Best<-c(Best, dat.temp$BestEvent_Pos[j])
    c.DN<-c(c.DN, dat.temp$ZAP_Both[j])
    NegAff<-c(NegAff, dat.temp$MAFS_NA[j])
    PosAff<-c(PosAff, dat.temp$MAFS_PA[j])
    Gender<-c(Gender, dat.temp$Gender[j])
    BFI_E<-c(BFI_E, dat.temp$Battery_BFI_E[j])
    BFI_A<-c(BFI_A, dat.temp$Battery_BFI_A[j])
    BFI_C<-c(BFI_C, dat.temp$Battery_BFI_C[j])
    BFI_O<-c(BFI_O, dat.temp$Battery_BFI_O[j])
    GenDep<-c(GenDep, dat.temp$Battery_IDASII64_GenDep[j])
  }
}

dat.study1<-data.frame(ID, 
                       Worst, 
                       Best, 
                       c.DN, 
                       NegAff,
                       PosAff,
                       Gender, 
                       BFI_A, 
                       BFI_C, 
                       BFI_E, 
                       BFI_O, 
                       GenDep)                
psych::describe(dat.study1)

dat.study1$Gender<-dat.study1$Gender-1  #Putting this in 0-1 coding
#Males = 1, Females = 0
dat.study1$c.DN<-as.numeric(scale(dat.study1$c.DN)) #Re-standardizing to give variable unit variance in models

##########################################################################################################
##########################################################################################################
#Creating a brms model for Study 1 - using as a discovery sample
##########################################################################################################
##########################################################################################################

#Multilevel imputation will take an incredibly long time (even when running on eleven cores)

#inspecting missing patterns in the data: 
md.pattern(dat.study1)

#Approximately 20% of the data is missing - will impute using a multilevel approach

#Exploring differences in missing vs. non-missing distributions of scores: 
#Missinginess is almost entirely clustered - will explore using single variable
dat.study1$miss<-ifelse(!is.na(dat.study1$PosAff), 'Complete', 'Missing')

g1<-ggplot(data=dat.study1, 
           aes(x=BFI_O, y=miss, fill=miss))+
  geom_density_ridges(alpha=.5)+
  ylab('')+
  theme_bw()+
  xlab('BFI Openess Scores (1-5)')+
  guides(fill = guide_legend(""))+
  scale_fill_brewer(palette = 'Set1')+
  geom_vline(xintercept = mean(dat.study1$BFI_O[dat.study1$miss=='Missing']), color=brewer.pal(3,'Set1')[2])+
  geom_vline(xintercept = mean(dat.study1$BFI_O[dat.study1$miss=='Complete']), color=brewer.pal(3,'Set1')[1])+
  annotate(geom = 'text', label=paste('Mean Complete =', 
                                      round(mean(dat.study1$BFI_O[dat.study1$miss=='Complete']), digits = 2)), 
           color=brewer.pal(3,'Set1')[1], 
           x=2.9, 
           y=.5)+
  annotate(geom = 'text', label=paste('Mean Missing =', 
                                      round(mean(dat.study1$BFI_O[dat.study1$miss=='Missing']), digits = 2)), 
           color=brewer.pal(3,'Set1')[2], 
           x=4.25, 
           y=.5)

g2<-ggplot(data=dat.study1, 
           aes(x=BFI_C, y=miss, fill=miss))+
  geom_density_ridges(alpha=.5)+
  ylab('')+
  theme_bw()+
  xlab('BFI Conscientiousness Scores (1-5)')+
  guides(fill = guide_legend(""))+
  scale_fill_brewer(palette = 'Set1')+
  geom_vline(xintercept = mean(dat.study1$BFI_C[dat.study1$miss=='Missing']), color=brewer.pal(3,'Set1')[2])+
  geom_vline(xintercept = mean(dat.study1$BFI_C[dat.study1$miss=='Complete']), color=brewer.pal(3,'Set1')[1])+
  annotate(geom = 'text', label=paste('Mean Complete =', 
                                      round(mean(dat.study1$BFI_C[dat.study1$miss=='Complete']), digits = 2)), 
           color=brewer.pal(3,'Set1')[1], 
           x=4.35, 
           y=.5)+
  annotate(geom = 'text', label=paste('Mean Missing =', 
                                      round(mean(dat.study1$BFI_C[dat.study1$miss=='Missing']), digits = 2)), 
           color=brewer.pal(3,'Set1')[2], 
           x=2.75, 
           y=.5)

g3<-ggplot(data=dat.study1, 
           aes(x=BFI_E, y=miss, fill=miss))+
  geom_density_ridges(alpha=.5)+
  ylab('')+
  theme_bw()+
  xlab('BFI Extraversion Scores (1-5)')+
  guides(fill = guide_legend(""))+
  scale_fill_brewer(palette = 'Set1')+
  geom_vline(xintercept = mean(dat.study1$BFI_E[dat.study1$miss=='Missing']), color=brewer.pal(3,'Set1')[2])+
  geom_vline(xintercept = mean(dat.study1$BFI_E[dat.study1$miss=='Complete']), color=brewer.pal(3,'Set1')[1])+
  annotate(geom = 'text', label=paste('Mean Complete =', 
                                      round(mean(dat.study1$BFI_E[dat.study1$miss=='Complete']), digits = 2)), 
           color=brewer.pal(3,'Set1')[1], 
           x=2.2, 
           y=.5)+
  annotate(geom = 'text', label=paste('Mean Missing =', 
                                      round(mean(dat.study1$BFI_E[dat.study1$miss=='Missing']), digits = 2)), 
           color=brewer.pal(3,'Set1')[2], 
           x=4, 
           y=.5)


g4<-ggplot(data=dat.study1, 
           aes(x=BFI_A, y=miss, fill=miss))+
  geom_density_ridges(alpha=.5)+
  ylab('')+
  theme_bw()+
  xlab('BFI Agreeableness Scores (1-5)')+
  guides(fill = guide_legend(""))+
  scale_fill_brewer(palette = 'Set1')+
  geom_vline(xintercept = mean(dat.study1$BFI_A[dat.study1$miss=='Missing']), color=brewer.pal(3,'Set1')[2])+
  geom_vline(xintercept = mean(dat.study1$BFI_A[dat.study1$miss=='Complete']), color=brewer.pal(3,'Set1')[1])+
  annotate(geom = 'text', label=paste('Mean Complete =', 
                                      round(mean(dat.study1$BFI_A[dat.study1$miss=='Complete']), digits = 2)), 
           color=brewer.pal(3,'Set1')[1], 
           x=4.5, 
           y=.5)+
  annotate(geom = 'text', label=paste('Mean Missing =', 
                                      round(mean(dat.study1$BFI_A[dat.study1$miss=='Missing']), digits = 2)), 
           color=brewer.pal(3,'Set1')[2], 
           x=3.15, 
           y=.5)

g5<-ggplot(data=dat.study1, 
           aes(x=c.DN, y=miss, fill=miss))+
  geom_density_ridges(alpha=.5)+
  ylab('')+
  theme_bw()+
  xlab('Centered DN scores')+
  guides(fill = guide_legend(""))+
  scale_fill_brewer(palette = 'Set1')+
  geom_vline(xintercept = mean(dat.study1$c.DN[dat.study1$miss=='Missing']), color=brewer.pal(3,'Set1')[2])+
  geom_vline(xintercept = mean(dat.study1$c.DN[dat.study1$miss=='Complete']), color=brewer.pal(3,'Set1')[1])+
  annotate(geom = 'text', label=paste('Mean Complete =', 
                                      round(mean(dat.study1$c.DN[dat.study1$miss=='Complete']), digits = 2)), 
           color=brewer.pal(3,'Set1')[1], 
           x=-.95, 
           y=.5)+
  annotate(geom = 'text', label=paste('Mean Missing =', 
                                      round(mean(dat.study1$c.DN[dat.study1$miss=='Missing']), digits = 2)), 
           color=brewer.pal(3,'Set1')[2], 
           x=.95, 
           y=.5)

g6<-ggplot(data=dat.study1, 
           aes(x=GenDep, y=miss, fill=miss))+
  geom_density_ridges(alpha=.5)+
  ylab('')+
  theme_bw()+
  xlab('IDAS Depression Scores')+
  guides(fill = guide_legend(""))+
  scale_fill_brewer(palette = 'Set1')+
  geom_vline(xintercept = mean(dat.study1$GenDep[dat.study1$miss=='Missing']), color=brewer.pal(3,'Set1')[2])+
  geom_vline(xintercept = mean(dat.study1$GenDep[dat.study1$miss=='Complete']), color=brewer.pal(3,'Set1')[1])+
  annotate(geom = 'text', label=paste('Mean Complete =', 
                                      round(mean(dat.study1$GenDep[dat.study1$miss=='Complete']), digits = 2)), 
           color=brewer.pal(3,'Set1')[1], 
           x=30, 
           y=.5)+
  annotate(geom = 'text', label=paste('Mean Missing =', 
                                      round(mean(dat.study1$GenDep[dat.study1$miss=='Missing']), digits = 2)), 
           color=brewer.pal(3,'Set1')[2], 
           x=55, 
           y=.5)
png(paste0(study1.graphics, '/Missing_Data_Exploration.png'), 
    res=600, 
    units='in', 
    height=10, 
    width =16)
cowplot::plot_grid(g1,g2,g3,g4,g5,g6, nrow = 3)
dev.off()

#Nothing obvious in this relative unsophisticated way of looking at the data
#Will Model missingness using a random effects approach. 
dat.study1$miss.dich<-ifelse(dat.study1$miss=='Missing', 1, 0)
dat.study1$zGenDep<-as.numeric(scale(dat.study1$GenDep))

fit.miss<-lme4::glmer(miss.dich~1+c.DN+Gender+zGenDep+BFI_A+BFI_C+BFI_E+BFI_O+
                        (1|ID), 
                      data=dat.study1, 
                      family='binomial',
                      control=glmerControl(optimizer="bobyqa",
                                           optCtrl=list(maxfun=2e5)))
summary(fit.miss)
#Suggests that the following variables are related to missingness: 
#   1. BFI Agreeableness (less likely to have missing data)
#   2. BFI Conscientiouness (less likely to have missing data)
#   3. BFI Openness (more likely to have missing data... interesting)
#=========================================================================================================
#Create initial matrix for convenience
#Prediciting missing values using variables that will go in the model (c.DN)
#Also including variables noted in the missing data analysis above as being related to missingness

fml<- Worst + Best + NegAff + PosAff  ~ 
  1 + c.DN + BFI_A + BFI_C + BFI_O + (1|ID)

#Set number of data sets to impute 
M<-20
imp<-panImpute(dat.study1, 
               formula=fml, 
               n.burn=10000, 
               n.iter = 5000, 
               m=M, 
               seed = 0716)

#plotting results - assessing convergence on final distributions
dat.imp<-mitmlComplete(imp)

#Inspecting imputation quality properties - especially interested in outcome measures  
dat.long<-data.frame()
for(i in 1:M){
  dat.temp<-dat.imp[[i]]
  dat.temp$IMP<-rep(i, length(dat.temp[,1]))
  dat.long<-rbind(dat.long, dat.temp)
}

dat.study1$IMP<-rep(0, length(dat.study1[,1]))
dat.long<-rbind(dat.study1, dat.long)
dat.long$Orig<-ifelse(dat.long$IMP<1, 'Original', 'Imputed')

#Plotting Negative Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=NegAff, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Negative Affect')+
  ylab('Density')+
  xlab('Momentary Negative Affect')
g1

png(paste0(study1.graphics, '/S1_Imputation_NegativeAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Positive Affect Imputation
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=PosAff, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Momentary Positive Affect')+
  ylab('Density')+
  xlab('Momentary Positive Affect')
g1

png(paste0(study1.graphics, '/S1_Imputation_PositiveAffect.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Deviations from Mean Worst Ratings
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=Worst, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Worst Event')+
  ylab('Density')+
  xlab('Worst Event Ratings')
g1

png(paste0(study1.graphics, '/S1_Imputation_WorstEvent.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#Plotting Deviations from Mean Best Ratings
g1<-ggplot()+
  geom_freqpoly(data = dat.long, aes(x=Best, group=IMP, color=Orig), stat = 'density')+
  theme(legend.title = element_blank())+
  ggtitle('Imputed Distributions of Best Event')+
  ylab('Density')+
  xlab('Best Event Ratings')
g1

png(paste0(study1.graphics, '/S1_Imputation_BestEvent.png'), res=300, 
    units='in', height=5, width = 8)
g1
dev.off()

#summarizing imputations - checking for convergence
sink(paste0(study1.out, '/Imputation_Summary.txt'))
summary(imp)
psych::describeBy(dat.long, group='IMP')
sink()    #warnings have to do with factor variables for the most part


#Looping through and constraining values to be 1 <= x <= 5.  
#Note that 20 is the number of imputed data sets:
for(i in 1:M){
  dat.imp[[i]]$T.NegAff<-ifelse(dat.imp[[i]]$NegAff>5, 5, dat.imp[[i]]$NegAff)
  dat.imp[[i]]$T.NegAff<-ifelse(dat.imp[[i]]$T.NegAff<1, 1, dat.imp[[i]]$T.NegAff)
  dat.imp[[i]]$T.PosAff<-ifelse(dat.imp[[i]]$PosAff>5, 5, dat.imp[[i]]$PosAff)
  dat.imp[[i]]$T.PosAff<-ifelse(dat.imp[[i]]$T.PosAff<1, 1, dat.imp[[i]]$T.PosAff)
  dat.imp[[i]]$T.Worst<-ifelse(dat.imp[[i]]$Worst>5, 5, dat.imp[[i]]$Worst)
  dat.imp[[i]]$T.Worst<-ifelse(dat.imp[[i]]$T.Worst<1, 1, dat.imp[[i]]$T.Worst)
  dat.imp[[i]]$T.Best<-ifelse(dat.imp[[i]]$Best>5, 5, dat.imp[[i]]$Best)
  dat.imp[[i]]$T.Best<-ifelse(dat.imp[[i]]$T.Best<1, 1, dat.imp[[i]]$T.Best)
}

#Quirky formatting problems - getting around it this way... 
#   The issue is that brm_multiple requires a list of datasets
#   It cannot recognize an mitml.list object though
#   So I made a list out of list
dat.study1.list<-list()

for(i in 1:M){
  dat.study1.list[[i]]<-dat.imp[[i]]
}

#Obtaining individual mean ratings for Best and Worst events
#Individually mean centering predictors to decompose between and within effects

for(i in 1:M){
  IDs<-unique(dat.study1$ID)
  DF.temp<-dat.study1.list[[i]]
  for(j in 1:length(IDs)){
    DF.temp$c.Worst[DF.temp$ID==IDs[j]]<-
      DF.temp$Worst[DF.temp$ID==IDs[j]] - 
      mean(DF.temp$Worst[DF.temp$ID==IDs[j]], na.rm=T)
    
    DF.temp$mean.Worst[DF.temp$ID==IDs[j]]<-rep(mean(DF.temp$Worst[DF.temp$ID==IDs[j]]))
    
    DF.temp$c.Best[DF.temp$ID==IDs[j]]<-
      DF.temp$Best[DF.temp$ID==IDs[j]] - 
      mean(DF.temp$Best[DF.temp$ID==IDs[j]], na.rm=T)
    
    DF.temp$mean.Best[DF.temp$ID==IDs[j]]<-rep(mean(DF.temp$Best[DF.temp$ID==IDs[j]]))
  }
  dat.study1.list[[i]]<-DF.temp
}

#---------------------------------------------------------------------------------------------------------
#NEGATIVE AFFECT MODEL
#---------------------------------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Negative unconditional model (i.e., random intercepts only)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Neg_ucm<-brms::brm_multiple(T.NegAff~1+(1|ID), 
                            data = dat.study1.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,2)', class='Intercept')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.95, 
                                           max_treedepth=15), 
                            save_model = paste0(stan.code, '/S1_ucm.stan'))
print(summary(Neg_ucm), digits=3)
plot(Neg_ucm)
print(sjstats::icc(Neg_ucm, posterior = T, prob=.90), digits=4)
bayes.to.txt(model = Neg_ucm,
             out.folder = study1.out, 
             file = 'S1_NegMood_UCM', 
             DF = dat.study1.list[[2]],
             tot.pars = 3, 
             impute = TRUE)

#--
Neg_lv1<-brms::brm_multiple(T.NegAff~1+c.Worst+c.Best+(1+c.Worst+c.Best|ID), 
                            data = dat.study1.list, 
                            family = 'lognormal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(0,5)', class='Intercept'), 
                                      set_prior('normal(0,5)', class='b'), 
                                      set_prior('lkj(2)', class='cor')),
                            warmup = 3000, 
                            iter = 6000,
                            thin = 3, 
                            chains = 4,
                            control = list(adapt_delta=.95, 
                                           max_treedepth=15), 
                            save_model = paste0(stan.code, '/S1_lv1.stan'))

print(summary(Neg_lv1), digits=3)
plot(Neg_lv1, N=10)
print(sjstats::icc(Neg_lv1, posterior = T, prob=.90), digits=4)

#Calculating within-subjects variance accounted for
wthn_var_Best_Worst<-(.1221-.1077)/.1221
wthn_var_Best_Worst

bayes.to.txt(model = Neg_lv1,
             out.folder = study1.out, 
             file = 'S1_NegMood_LV1', 
             DF = dat.study1.list[[2]],
             tot.pars = 10, 
             impute = TRUE)

#--
Neg_lv2_Best<-brms::brm_multiple(T.NegAff~1+c.Worst+c.Best+
                                   mean.Best+
                                   (1+c.Worst+c.Best|ID), 
                                 data = dat.study1.list, 
                                 family = 'lognormal',
                                 prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                           set_prior('normal(0,5)', class='Intercept'), 
                                           set_prior('normal(0,5)', class='b'), 
                                           set_prior('lkj(2)', class='cor')),
                                 warmup = 3000, 
                                 iter = 6000,
                                 thin = 3, 
                                 chains = 4,
                                 control = list(adapt_delta=.95, 
                                                max_treedepth=15), 
                                 save_model = paste0(stan.code, '/S1_lv2a.stan'))
print(summary(Neg_lv2_Best), digits = 3)
print(sjstats::icc(Neg_lv2_Best, posterior = T, prob=.90), digits=4)

#Looking at variance accounted for by the contextual effect here - and it is not significant
#Comparing the between residuals for the Neg_lv1 and Neg_lv2_Best - no difference mean's 0 variability
btw_var_mean.Best<-0

bayes.to.txt(model = Neg_lv2_Best,
             out.folder = study1.out, 
             file = 'S1_NegMood_LV2_Best', 
             DF = dat.study1.list[[2]],
             tot.pars = 11, 
             impute = TRUE)

Neg_lv2_Worst<-brms::brm_multiple(T.NegAff~1+c.Worst+c.Best+
                                    mean.Worst+
                                    (1+c.Worst+c.Best|ID), 
                                  data = dat.study1.list, 
                                  family = 'lognormal',
                                  prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                            set_prior('normal(0,5)', class='Intercept'), 
                                            set_prior('normal(0,5)', class='b'), 
                                            set_prior('lkj(2)', class='cor')),
                                  warmup = 3000, 
                                  iter = 6000,
                                  thin = 3, 
                                  chains = 4,
                                  control = list(adapt_delta=.95, 
                                                 max_treedepth=15), 
                                  save_model = paste0(stan.code, '/S1_lv2a.stan'))

print(summary(Neg_lv2_Worst), digits = 3)
print(sjstats::icc(Neg_lv1, posterior = T, prob=.90), digits=4)
print(sjstats::icc(Neg_lv2_Worst, posterior = T, prob=.90), digits=4)

#Between-Subjects - Total Variance explained by average negative events
btw_var_mean.Best<-(.0901-.0670)/.0901

bayes.to.txt(model = Neg_lv2_Worst,
             out.folder = study1.out, 
             file = 'S1_NegMood_LV2_Worst', 
             DF = dat.study1.list[[2]],
             tot.pars = 11, 
             impute = TRUE)

Neg_lv2_DN<-brms::brm_multiple(T.NegAff~1+c.Worst+c.Best+
                                 c.DN+
                                 (1+c.Worst+c.Best|ID), 
                               data = dat.study1.list, 
                               family = 'lognormal',
                               prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                         set_prior('normal(0,5)', class='Intercept'), 
                                         set_prior('normal(0,5)', class='b'), 
                                         set_prior('lkj(2)', class='cor')),
                               warmup = 3000, 
                               iter = 6000,
                               thin = 3, 
                               chains = 4,
                               control = list(adapt_delta=.95, 
                                              max_treedepth=15), 
                               save_model = paste0(stan.code, '/S1_lv2a.stan'))
print(summary(Neg_lv2_DN), digits = 3)
print(sjstats::icc(Neg_lv1, posterior = T, prob=.90), digits=4)
print(sjstats::icc(Neg_lv2_DN, posterior = T, prob=.90), digits=4)

#Overall DN - (not controlling for anything...)
btw_var_DN<-(.0901-.0720)/.0901

bayes.to.txt(model = Neg_lv2_DN,
             out.folder = study1.out, 
             file = 'S1_NegMood_LV2_DN', 
             DF = dat.study1.list[[2]],
             tot.pars = 11, 
             impute = TRUE)

Neg_lv2_DN_Best<-brms::brm_multiple(T.NegAff~1+c.Worst+c.Best+
                                      c.DN+mean.Best+
                                      (1+c.Worst+c.Best|ID), 
                                    data = dat.study1.list, 
                                    family = 'lognormal',
                                    prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                              set_prior('normal(0,5)', class='Intercept'), 
                                              set_prior('normal(0,5)', class='b'), 
                                              set_prior('lkj(2)', class='cor')),
                                    warmup = 3000, 
                                    iter = 6000,
                                    thin = 3, 
                                    chains = 4,
                                    control = list(adapt_delta=.95, 
                                                   max_treedepth=15), 
                                    save_model = paste0(stan.code, '/S1_lv2a.stan'))

bayes.to.txt(model = Neg_lv2_DN_Best,
             out.folder = study1.out, 
             file = 'S1_NegMood_LV2_DN_Best', 
             DF = dat.study1.list[[2]],
             tot.pars = 11, 
             impute = TRUE)

Neg_lv2_DN_Worst<-brms::brm_multiple(T.NegAff~1+c.Worst+c.Best+
                                       c.DN+mean.Worst+
                                       (1+c.Worst+c.Best|ID), 
                                     data = dat.study1.list, 
                                     family = 'lognormal',
                                     prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                               set_prior('normal(0,5)', class='Intercept'), 
                                               set_prior('normal(0,5)', class='b'), 
                                               set_prior('lkj(2)', class='cor')),
                                     warmup = 3000, 
                                     iter = 6000,
                                     thin = 3, 
                                     chains = 4,
                                     control = list(adapt_delta=.95, 
                                                    max_treedepth=15), 
                                     save_model = paste0(stan.code, '/S1_lv2a.stan'))

print(summary(Neg_lv2_DN_Worst), digits = 3)
print(sjstats::icc(Neg_lv1, posterior = T, prob=.90), digits=4)
print(sjstats::icc(Neg_lv2_DN_Worst, posterior = T, prob=.90), digits=4)
btw_var_un_DN_Worst<-(.0670-.0589)/.0901

bayes.to.txt(model = Neg_lv2_DN_Worst,
             out.folder = study1.out, 
             file = 'S1_NegMood_LV2_DN_Worst', 
             DF = dat.study1.list[[2]],
             tot.pars = 11, 
             impute = TRUE)

Neg_lv2_Evnt<-brms::brm_multiple(T.NegAff~1+c.Worst+c.Best+
                                   mean.Worst+mean.Best+
                                   (1+c.Worst+c.Best|ID), 
                                 data = dat.study1.list, 
                                 family = 'lognormal',
                                 prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                           set_prior('normal(0,5)', class='Intercept'), 
                                           set_prior('normal(0,5)', class='b'), 
                                           set_prior('lkj(2)', class='cor')),
                                 warmup = 3000, 
                                 iter = 6000,
                                 thin = 3, 
                                 chains = 4,
                                 control = list(adapt_delta=.95, 
                                                max_treedepth=15), 
                                 save_model = paste0(stan.code, '/S1_lv2a.stan'))
print(summary(Neg_lv2_Evnt), digits=3)
print(sjstats::icc(Neg_lv2_Evnt, posterior = T, prob=.90), digits=4)

#Calculating amount of between-subjects variance accounted for by average event ratings
btw_var_mean.Best_mean.Worst<-(.0901-.0574)/.0901

bayes.to.txt(model = Neg_lv2_Evnt,
             out.folder = study1.out, 
             file = 'S1_NegMood_LV2_Evnt', 
             DF = dat.study1.list[[2]],
             tot.pars = 12, 
             impute = TRUE)

Neg_lv2b<-brms::brm_multiple(T.NegAff~1+c.Worst+c.Best+
                               mean.Worst+mean.Best+c.DN+
                               (1+c.Worst+c.Best|ID), 
                             data = dat.study1.list, 
                             family = 'lognormal',
                             prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                       set_prior('normal(0,5)', class='Intercept'), 
                                       set_prior('normal(0,5)', class='b'), 
                                       set_prior('lkj(2)', class='cor')),
                             warmup = 3000, 
                             iter = 6000,
                             thin = 3, 
                             chains = 4,
                             control = list(adapt_delta=.95, 
                                            max_treedepth=15), 
                             save_model = paste0(stan.code, '/S1_lv2b.stan'))
print(summary(Neg_lv2b), digits=6)
print(bayes_R2(Neg_lv2b), digits = 7)
print(sjstats::icc(Neg_lv2b, posterior = T, prob=.90), digits=8)

#Unique Variance accounted for in between subjects mood differences (controlling for average event)
btw_var_DN_un<-(.0574-.0511)/.0901

bayes.to.txt(model = Neg_lv2b,
             out.folder = study1.out, 
             file = 'S1_NegMood_LV2_ALL', 
             DF = dat.study1.list[[2]],
             tot.pars = 13, 
             impute = TRUE)

Neg_cross_DN_Best<-brms::brm_multiple(T.NegAff~1+c.Worst+c.Best+
                                        mean.Worst+mean.Best+c.DN+
                                        c.Best:c.DN+
                                        (1+c.Worst+c.Best|ID), 
                                      data = dat.study1.list, 
                                      family = 'lognormal',
                                      prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                                set_prior('normal(0,5)', class='Intercept'), 
                                                set_prior('normal(0,5)', class='b'), 
                                                set_prior('lkj(2)', class='cor')),
                                      warmup = 3000, 
                                      iter = 6000,
                                      thin = 3, 
                                      chains = 4,
                                      control = list(adapt_delta=.95, 
                                                     max_treedepth=15), 
                                      save_model = paste0(stan.code, '/S1_cross.stan'))
print(summary(Neg_cross_DN_Best), digits=5)
print(summary(Neg_lv1), digits=5)

bayes.to.txt(model = Neg_cross_DN_Best,
             out.folder = study1.out, 
             file = 'S1_NegMood_CROSS_DN_Best', 
             DF = dat.study1.list[[2]],
             tot.pars = 14, 
             impute = TRUE)

#Negative (Anxious) Mood Reactivity in response to positive events accounted for by DN
sd_Best_null<-.03454 
sd_Best_DN<-.03241

#Reactivity variance 
vcov(Neg_cross_DN_Best)
sum(vcov(Neg_cross_DN_Best))

print(sjstats::icc(Neg_lv1, posterior = T, prob=.90), digits=4)
print(sjstats::icc(Neg_cross_DN_Best, posterior = T, prob=.90), digits=4)

Neg_cross_DN_Worst<-brms::brm_multiple(T.NegAff~1+c.Worst+c.Best+
                                         mean.Worst+mean.Best+c.DN+
                                         c.Worst:c.DN+
                                         (1+c.Worst+c.Best|ID), 
                                       data = dat.study1.list, 
                                       family = 'lognormal',
                                       prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                                 set_prior('normal(0,5)', class='Intercept'), 
                                                 set_prior('normal(0,5)', class='b'), 
                                                 set_prior('lkj(2)', class='cor')),
                                       warmup = 3000, 
                                       iter = 6000,
                                       thin = 3, 
                                       chains = 4,
                                       control = list(adapt_delta=.95, 
                                                      max_treedepth=15), 
                                       save_model = paste0(stan.code, '/S1_cross.stan'))

bayes.to.txt(model = Neg_cross_DN_Worst,
             out.folder = study1.out, 
             file = 'S1_NegMood_CROSS_DN_Worst', 
             DF = dat.study1.list[[2]],
             tot.pars = 14, 
             impute = TRUE)


Neg_cross<-brms::brm_multiple(T.NegAff~1+c.Worst+c.Best+
                                mean.Worst+mean.Best+c.DN+
                                c.Best:c.DN+c.Worst:c.DN+
                                (1+c.Worst+c.Best|ID), 
                              data = dat.study1.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                        set_prior('normal(0,5)', class='Intercept'), 
                                        set_prior('normal(0,5)', class='b'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 3000, 
                              iter = 6000,
                              thin = 3, 
                              chains = 4,
                              control = list(adapt_delta=.95, 
                                             max_treedepth=15), 
                              save_model = paste0(stan.code, '/S1_cross.stan'))
print(bayes_R2(Neg_cross), digits=7)
print(sjstats::icc(Neg_cross), digits=8)
bayes.to.txt(model = Neg_cross,
             out.folder = study1.out, 
             file = 'S1_NegMood_CROSS_ALL', 
             DF = dat.study1.list[[2]],
             tot.pars = 15, 
             impute = TRUE)

#Getting posterior predictive distribution (a model check for Negative Mood)
png(paste0(study1.graphics, '/S1_Supplemental_FigureX_final_ppdens_NA.png'), 
    res = 600, 
    units = 'in', 
    height = 6, 
    width = 6)
pp_check(Neg_cross, nsamples=100)+
  ggtitle('Study 1: Posterior Predictive Distribution for Negative Affect')+
  scale_x_continuous(limits=c(0,6))
dev.off()


#========================================================================================
#NEED TO COME BACK THROUGH AND CLEAN UP PLOTTING CODE BELOW 
#========================================================================================
#Reactivity Variance from model R2's (may want to build a function for this eventually - hard code for now to test out)
tot_var_react<-.4750223-.4739584 

#Going to build out data frame manually for graphing partialed out variance from model 
#The intraclass correlation from UCM - removing cross-level interaction variance (very smally ~.1%) 
#The idea here is to decompose "pure" between subjects variance, pure within-subjects variance, and the cross-level interaction variance
tot_btw<-.42006       
tot_wthn<-1-tot_btw      #"Pure" within-subjects variance is what remains here 

#Within-subjects variables (i.e., events) accounting for total variance
wthn_best_worst_tot<-(.12207-.10770)/.12207*tot_wthn  #Based on reduction of within-subjects residual variance
wthn_res_tot<-(1-(.12207-.10770)/.12207)*tot_wthn     #Suggests that approximately 50% of total variance is based on within-subjects residuals

wthn_best_worst_wthn<-wthn_best_worst_tot/tot_wthn    #Converting back to within-subjects variability only
wthn_res_wthn<-1-wthn_best_worst_wthn

#Between-subjects variables
btw_context_tot<-(.09007-.05740)/.09007*tot_btw   #Overall effect of living in varying contexts (i.e., high negative events, low positive events)
#.09007 comes from level 1 model - after controlling for within-subjects differences
#Interpretation is that ~15% of total variance is attributable to between-subjects differencs in context

btw_DN_tot<-(.09007-.07205)/.09007*tot_btw        #Overall effect of DN in predicting momentary mood ~8% of overall variance (how much unique?)

btw_DN_uni_tot<-(.05740-.05110)/.09007*tot_btw    #when combined with btw_context_tot = total between-sujects variance explained (~3%)
btw_context_uni_tot<-(.07205-.05110)/.09007*tot_btw   #when combined with btw_DN_tot = total between-subjects variance
btw_shared_tot<-btw_DN_tot-btw_DN_uni_tot         #can arrive at this value multiple ways 
btw_all_tot<-(.09007-.05110)/.09007*tot_btw
btw_res_tot<-(1-(.09007-.05110)/.09007)*tot_btw

btw_context_btw<-btw_context_tot/tot_btw              #Converting back to between subjects variability
btw_DN_btw<-btw_DN_tot/tot_btw                        #Overall effect of DN in predicting momentary mood ~8% of overall variance (how much unique?)
btw_DN_uni_btw<-btw_DN_uni_tot/tot_btw                #when combined with btw_context_tot = total between-sujects variance explained (~3%)
btw_context_uni_btw<-btw_context_uni_tot/tot_btw      #when combined with btw_DN_tot = total between-subjects variance
btw_shared_btw<-(btw_DN_tot-btw_DN_uni_tot)/tot_btw   #can arrive at this value multiple ways 
btw_all_btw<-(.09007-.05110)/.09007                   #should equal btw_DN_uni_btw+btw_context_uni_btw+btw_shared_btw
btw_res_btw<-1-(.09007-.05110)/.09007

#Variation in momentary reactivity - Need to calculate (will include in text - no graph needed for this value)

#How much does DN account for variation in momentary reactions to events... (which definitely exist)

#Total contribution of DN = unique DN effect + reactivity effect (which is very small) + shared contextual variance (potential indirect pathway...)
#Can think of this as the pathways of DN to negative mood
DN_total<-btw_DN_tot+tot_var_react
DN_unique_DN<-btw_DN_uni_tot/DN_total
DN_indirect_DN<-btw_shared_tot/DN_total
DN_react_DN<-tot_var_react/DN_total

DF_tot_var<-data.frame(Value = c(btw_context_uni_tot,
                                 btw_shared_tot, 
                                 btw_DN_uni_tot, 
                                 btw_res_tot,
                                 wthn_best_worst_tot,
                                 wthn_res_tot, 
                                 wthn_best_worst_wthn, 
                                 wthn_res_wthn, 
                                 btw_context_uni_btw,
                                 btw_shared_btw, 
                                 btw_DN_uni_btw, 
                                 btw_res_btw
), 
Component = c(rep('Total', 6), 
              rep('Within-Subjects', 2), 
              rep('Between-Subjects', 4)), 
Term = c('Event Context Variance', 
         'DN-Event Context Shared Variance', 
         'DN Variance', 
         'Unmodeled Between-Subjects Variance',
         'Momentary Event Variance', 
         'Unmodeled Momentary Variance',
         'Momentary Event Variance', 
         'Unmodeled Momentary Variance',
         'Event Context Variance', 
         'DN-Event Context Shared Variance', 
         'DN Variance', 
         'Unmodeled Between-Subjects Variance'
)
)

#Note that for the pie-charts to line up as expected I have to use the reverse factor order 
#Same function, same values, different order compared to below 
DF_tot_var$Term<-forcats::fct_relevel(DF_tot_var$Term, 
                                      c('Unmodeled Momentary Variance',
                                        'Momentary Event Variance', 
                                        'Unmodeled Between-Subjects Variance',
                                        'DN Variance',
                                        'DN-Event Context Shared Variance',
                                        'Event Context Variance'))

sum(DF_tot_var$Value)   #Should sum to 3 (1 for each variance component: total, within, and between)

#Extracting distribution from models (will try a violin plot for posterior)

y_hat<-posterior_predict(Neg_cross)
y_hat<-colMeans(y_hat)                #Getting point estimates for each observation - don't need all 80,000 samples for that
hist(y_hat)

#Secondary plotting DF 
DF_plot_temp<-data.frame(Term = c('Event Context Variance', 
                                  'DN-Event Context Shared Variance', 
                                  'DN Variance', 
                                  'Unmodeled Between-Subjects Variance',
                                  'Momentary Event Variance', 
                                  'Unmodeled Momentary Variance'),
                         ymin = c(min(y_hat), 
                                  as.numeric(quantile(y_hat, probs = DF_tot_var$Value[1])),
                                  as.numeric(quantile(y_hat, probs = sum(DF_tot_var$Value[1:2]))),
                                  as.numeric(quantile(y_hat, probs = sum(DF_tot_var$Value[1:3]))),
                                  as.numeric(quantile(y_hat, probs = sum(DF_tot_var$Value[1:4]))), 
                                  as.numeric(quantile(y_hat, probs = sum(DF_tot_var$Value[1:5])))
                         ), 
                         ymax = c(as.numeric(quantile(y_hat, probs = DF_tot_var$Value[1])), 
                                  as.numeric(quantile(y_hat, probs = sum(DF_tot_var$Value[1:2]))),
                                  as.numeric(quantile(y_hat, probs = sum(DF_tot_var$Value[1:3]))),
                                  as.numeric(quantile(y_hat, probs = sum(DF_tot_var$Value[1:4]))),
                                  as.numeric(quantile(y_hat, probs = sum(DF_tot_var$Value[1:5]))),
                                  as.numeric(quantile(y_hat, probs = sum(DF_tot_var$Value[1:6])))
                         ), 
                         stringsAsFactors = FALSE
)

DF_plot_temp$y_pos<-rowMeans(DF_plot_temp[,2:3])
DF_plot_temp$x_pos<-rep(1.75)
DF_plot_temp$Value <- DF_tot_var$Value[1:6]

DF_plot_temp$Term<-forcats::fct_relevel(DF_plot_temp$Term, 
                                        c('Event Context Variance', 
                                          'DN-Event Context Shared Variance', 
                                          'DN Variance', 
                                          'Unmodeled Between-Subjects Variance', 
                                          'Momentary Event Variance', 
                                          'Unmodeled Momentary Variance'))

#Need to extract the violin shape here... 
g1<-ggplot()+
  geom_violin(aes(y=y_hat, x=1),
              bw = .25)

p_build<-ggplot_build(g1)$data[[1]]

p_build<-transform(p_build,
                   xminv = x - violinwidth * (x - xmin),
                   xmaxv = x + violinwidth * (xmax - x))

p_build<-rbind(plyr::arrange(transform(p_build, x = xminv), y),
               plyr::arrange(transform(p_build, x = xmaxv), -y))

p_build$fill_group<-rep(NA)
p_build$fill_group[p_build$y>=DF_plot_temp$ymin[1] & p_build$y<=DF_plot_temp$ymax[1]]<-DF_plot_temp$Term[1]
p_build$fill_group[p_build$y>DF_plot_temp$ymin[2] & p_build$y<=DF_plot_temp$ymax[2]]<-DF_plot_temp$Term[2]
p_build$fill_group[p_build$y>DF_plot_temp$ymin[3] & p_build$y<=DF_plot_temp$ymax[3]]<-DF_plot_temp$Term[3]
p_build$fill_group[p_build$y>DF_plot_temp$ymin[4] & p_build$y<=DF_plot_temp$ymax[4]]<-DF_plot_temp$Term[4]
p_build$fill_group[p_build$y>DF_plot_temp$ymin[5] & p_build$y<=DF_plot_temp$ymax[5]]<-DF_plot_temp$Term[5]
p_build$fill_group[p_build$y>DF_plot_temp$ymin[6] & p_build$y<=DF_plot_temp$ymax[6]]<-DF_plot_temp$Term[6]
table(p_build$fill_group)

p_build$group1 <- with(p_build,interaction(factor(group),factor(fill_group)))

table(p_build$group1)

p_build$fill_group<-forcats::fct_relevel(p_build$fill_group, 
                                         c('Event Context Variance', 
                                           'DN-Event Context Shared Variance', 
                                           'DN Variance', 
                                           'Unmodeled Between-Subjects Variance', 
                                           'Momentary Event Variance', 
                                           'Unmodeled Momentary Variance'))
palette_1 = c(paste0(brewer.pal(11, "Spectral"), "90")[2:5], 
              paste0(brewer.pal(11, "Spectral"), "90")[10:9])

p_fill <- ggplot() + 
  geom_violin(aes(y=y_hat, x=1),
              bw = .25, 
              alpha = .75, 
              fill='gray80')+
  geom_rect(aes(xmin = 1.70, 
                xmax = 2.30, 
                ymin = min(DF_plot_temp$ymin), 
                ymax = max(DF_plot_temp$ymax)), 
            fill = 'gray80', 
            alpha = .75)+
  geom_polygon(data = p_build,
               aes(x = x,y = y,group = fill_group, fill = fill_group, color=fill_group),
               alpha=.85)+
  geom_label(data = DF_plot_temp, 
             aes(x = x_pos, 
                 y = y_pos, 
                 label = paste0(round(Value*100, digits = 2), '%'), 
                 fill = Term), 
             show.legend = FALSE)+
  theme_classic()+
  scale_fill_manual(name = 'Variance Component', 
                    values = palette_1)+
  scale_color_manual(name = 'Variance Component', 
                     values = palette_1)+
  annotate(geom = 'text', x = 2, y = 4.5, label = 'Total Variance Explained', 
           size = 4.25)+
  ylab('Posterior Distribution of Momentary Negative Mood')+
  xlab('')+
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 14), 
        title = element_text(size = 16))+
  ggtitle('Momentary Negative Mood Variance Decomposition')

p_fill

#Note used edit function to tinker with values for aesthetic purposes
DF_plot_temp<-edit(DF_plot_temp)

png(paste0(study1.graphics, '/S1_NegMood_TotVar_Violin.png'), 
    res=900, 
    units = 'in', 
    height = 14, 
    width = 10.5)
p_fill
dev.off()

#Have to create a new variable to get this plot to work
DF_tot_var<- DF_tot_var %>% 
  group_by(Component) %>%
  arrange(Component) %>%
  mutate(pos = cumsum(Value) - Value/2)

#Creating Pie Charts for Within- and Between-Subjects Variance
pie_between<-ggplot(data = DF_tot_var[DF_tot_var$Component=='Between-Subjects',])+
  geom_bar(aes(x=Component, y=Value, fill = Term), stat = 'identity', show.legend = FALSE)+
  coord_polar('y', start=0)+
  theme_minimal()+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        panel.border = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        title = element_text(size = 16))+
  scale_fill_manual(values = palette_1[1:4])+
  geom_text(data = DF_tot_var[DF_tot_var$Component=='Between-Subjects',], 
            aes(label = paste0(round(Value*100, digits = 2), '%'), 
                y = pos, 
                x = Component))+
  ggtitle('Between-Subjects Variance Decomposition')

pie_between

png(paste0(study1.graphics, '/S1_NegMood_between_pie.png'), 
    res = 900, 
    units = 'in', 
    height = 8, 
    width = 8)
pie_between
dev.off()

#Within Subjects Version
pie_within<-ggplot(data = DF_tot_var[DF_tot_var$Component=='Within-Subjects',])+
  geom_bar(aes(x=Component, y=Value, fill = Term), stat = 'identity', show.legend = FALSE)+
  coord_polar('y', start=0)+
  theme_minimal()+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        panel.border = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        title = element_text(size = 16))+
  scale_fill_manual(values = palette_1[5:6])+
  geom_text(data = DF_tot_var[DF_tot_var$Component=='Within-Subjects',], 
            aes(label = paste0(round(Value*100, digits = 2), '%'), 
                y = pos, 
                x = Component))+
  ggtitle('Within-Subjects Variance Decomposition')

pie_within

png(paste0(study1.graphics, '/S1_NegMood_within_pie.png'), 
    res = 900, 
    units = 'in', 
    height = 8, 
    width = 8)
pie_within
dev.off()

#Combining plots
png(paste0(study1.graphics, '/Test_comb_vardecomp.png'), 
    res=900, 
    units = 'in', 
    height = 12, 
    width = 20)
ggarrange(p_fill,
          labels = 'A',
          ncol=2,
          ggarrange(pie_between, 
                    pie_within,
                    nrow = 2, 
                    labels = c('B', 'C')))
dev.off()

#And now for how DN contributes to negative mood... (three routes - general, shared variance with exposure, reactivity)
#Attemping a modification of a Sankey Plot: 
#Will look to make three total Nodes 
#First are the different Components linking to between and within-subject variance 
#Then Between & Within Subject Variance linking to Total Variance - trying just total DN Effect... 
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

palette = c(paste0(brewer.pal(9, "Blues"), 90)[5:8], 
            paste0(brewer.pal(9, "Reds"), 90)[7:8], 
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

png(paste0(study1.graphics, '/test_river.png'), 
    units = 'in', 
    res = 900, 
    height = 10, 
    width = 10)
riverplot(riv, 
          nodewidth = 3, 
          plot_area = .95)
title(ylab = 'Study 1 - Riverplot of Total Variance Decomposition Derived from Negative Mood Model')
dev.off()


#Getting separate DN only River plot decomposing the ways in which DN contributes to mood
DN_total<-btw_DN_tot+tot_var_react
DN_unique_DN<-btw_DN_uni_tot/DN_total
DN_indirect_DN<-btw_shared_tot/DN_total
DN_react_DN<-tot_var_react/DN_total

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

palette = c(paste0(brewer.pal(11, "Spectral"))[2:3], 
            paste0(brewer.pal(11, "Spectral"))[6], 
            paste0(brewer.pal(11, "Spectral"))[4])

styles = lapply(nodes$y, function(n) {
  list(col = palette[n], lty = 0, textcol = "black")
})

names(styles) = nodes$ID

riv<-makeRiver(nodes = nodes, 
               edges =  River_DF, 
               node_styles = styles)

png(paste0(study1.graphics, '/test_riverDN.png'), 
    units = 'in', 
    res = 900, 
    height = 10, 
    width = 10)
riverplot(riv, 
          nodewidth = 3, 
          plot_area = .95)
title(ylab = 'Study 1 - Riverplot of DN Effect Decomposition Derived from Negative Mood Model')
dev.off()

postscript(paste0(study1.graphics, '/test_riverDN.eps'), 
           width = 10, 
           height = 10, 
           horizontal = FALSE, 
           onefile = FALSE)
riverplot(riv, 
          nodewidth = 3, 
          plot_area = .95)
title(ylab = 'Study 1 - Riverplot of DN Effect Decomposition Derived from Negative Mood Model')
dev.off()

#---------------------------------------------------------------------------------------------------------
#POSITIVE AFFECT MODELS
#---------------------------------------------------------------------------------------------------------
mean(dat.study1$PosAff, na.rm=TRUE)
#using mean of original scale to have as a location prior for grand intercept - goal is to aid convergence 
#Standard deviation of the prior is larger than observed so should allow ample, reasonable parameter space
Pos_ucm<-brms::brm_multiple(T.PosAff~1+(1|ID), 
                            data = dat.study1.list, 
                            family = 'normal',
                            prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                      set_prior('normal(3.104,2)', class='Intercept')),
                            warmup = 2000, 
                            iter = 3000, 
                            chains = 3,
                            control = list(adapt_delta=.99, 
                                           max_treedepth=15), 
                            save_model = paste0(stan.code, '/S1_ucm.stan'))
gc()

Pos_lv2_DN<-brms::brm_multiple(T.PosAff~1+c.Worst+c.Best+
                                 c.DN+
                                 (1+c.Worst+c.Best|ID), 
                               data = dat.study1.list, 
                               family = 'normal',
                               prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                         set_prior('normal(3.104,2)', class='Intercept'), 
                                         set_prior('normal(0,5)', class='b'), 
                                         set_prior('lkj(2)', class='cor')),
                               warmup = 2000, 
                               iter = 3000, 
                               chains = 3,
                               control = list(adapt_delta=.99, 
                                              max_treedepth=15), 
                               save_model = paste0(stan.code, '/S1_lv2a.stan'))
gc()

Pos_lv2_Evnt<-brms::brm_multiple(T.PosAff~1+c.Worst+c.Best+
                                   mean.Worst+mean.Best+
                                   (1+c.Worst+c.Best|ID), 
                                 data = dat.study1.list, 
                                 family = 'normal',
                                 prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                           set_prior('normal(3.104,2)', class='Intercept'), 
                                           set_prior('normal(0,5)', class='b'), 
                                           set_prior('lkj(2)', class='cor')),
                                 warmup = 2000, 
                                 iter = 3000, 
                                 chains = 3,
                                 control = list(adapt_delta=.99, 
                                                max_treedepth=15), 
                                 save_model = paste0(stan.code, '/S1_lv2a.stan'))
gc()

Pos_lv2_All<-brms::brm_multiple(T.PosAff~1+c.Worst+c.Best+
                                  mean.Worst+mean.Best+c.DN+
                                  (1+c.Worst+c.Best|ID), 
                                data = dat.study1.list, 
                                family = 'normal',
                                prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                          set_prior('normal(3.104,2)', class='Intercept'), 
                                          set_prior('normal(0,5)', class='b'), 
                                          set_prior('lkj(2)', class='cor')),
                                warmup = 2000, 
                                iter = 3000, 
                                chains = 3,
                                control = list(adapt_delta=.99, 
                                               max_treedepth=15), 
                                save_model = paste0(stan.code, '/S1_lv2b.stan'))
gc()

Pos_cross<-brms::brm_multiple(T.PosAff~1+c.Worst+c.Best+
                                mean.Worst+mean.Best+c.DN+
                                c.Best:c.DN+c.Worst:c.DN+
                                (1+c.Worst+c.Best|ID), 
                              data = dat.study1.list, 
                              family = 'normal',
                              prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                        set_prior('normal(3.104,2)', class='Intercept'), 
                                        set_prior('normal(0,5)', class='b'), 
                                        set_prior('lkj(2)', class='cor')),
                              warmup = 2000, 
                              iter = 3000, 
                              chains = 3,
                              control = list(adapt_delta=.99, 
                                             max_treedepth=15), 
                              save_model = paste0(stan.code, '/S1_cross.stan'))
gc()

#Getting posterior predictive distribution (a model check for Positive Mood)
png(paste0(study1.graphics, '/S1_Supplemental_FigureX_final_ppdens_PA.png'), 
    res = 600, 
    units = 'in', 
    height = 6, 
    width = 6)
pp_check(Pos_cross, nsamples=100)+
  ggtitle('Study 1: Posterior Predictive Distribution for Positive Affect')
dev.off()

#------------------------------------------------------------------------------
#Calculating variance terms of interest:
#Reactivity Variance from model R2's (may want to build a function for this eventually - hard code for now to test out)
tot_var_react<-.62169-.62157

#Going to build out data frame manually for graphing partialed out variance from model 
#The intraclass correlation from UCM - removing cross-level interaction variance (very smally ~.1%) 
#The idea here is to decompose "pure" between subjects variance, pure within-subjects variance, and the cross-level interaction variance
tot_btw<-.51139    
tot_wthn<-1-tot_btw      #"Pure" within-subjects variance is what remains here 

#Within-subjects variables (i.e., events) accounting for total variance
wthn_best_worst_tot<-(.68426-.49353)/.68426*tot_wthn  #Based on reduction of within-subjects residual variance
wthn_res_tot<-(1-(.68426-.49353)/.68426)*tot_wthn     #Suggests that approximately 50% of total variance is based on within-subjects residuals

wthn_best_worst_wthn<-wthn_best_worst_tot/tot_wthn    #Converting back to within-subjects variability only
wthn_res_wthn<-1-wthn_best_worst_wthn

#Between-subjects variables
btw_context_tot<-(.69463-.44546)/.69463*tot_btw   #Overall effect of living in varying contexts (i.e., high negative events, low positive events)
#.09007 comes from level 1 model - after controlling for within-subjects differences
#Interpretation is that ~15% of total variance is attributable to between-subjects differencs in context

btw_DN_tot<-(.69463-.60130)/.69463*tot_btw        #Overall effect of DN in predicting momentary mood ~8% of overall variance (how much unique?)
btw_DN_uni_tot<-(.44546-.40244)/.69463*tot_btw    
btw_context_uni_tot<-(.60130-.40244)/.69463*tot_btw   #when combined with btw_DN_tot = total between-subjects variance
btw_shared_tot<-btw_DN_tot-btw_DN_uni_tot         #can arrive at this value multiple ways 
btw_all_tot<-(.69463-.40244)/.69463*tot_btw
btw_res_tot<-(1-(.69463-.40244)/.69463)*tot_btw

btw_context_btw<-btw_context_tot/tot_btw              #Converting back to between subjects variability
btw_DN_btw<-btw_DN_tot/tot_btw                        #Overall effect of DN in predicting momentary mood ~8% of overall variance (how much unique?)
btw_DN_uni_btw<-btw_DN_uni_tot/tot_btw                #when combined with btw_context_tot = total between-sujects variance explained (~3%)
btw_context_uni_btw<-btw_context_uni_tot/tot_btw      #when combined with btw_DN_tot = total between-subjects variance
btw_shared_btw<-(btw_DN_tot-btw_DN_uni_tot)/tot_btw   #can arrive at this value multiple ways 
btw_all_btw<-(.69463-.40244)/.69463                  #should equal btw_DN_uni_btw+btw_context_uni_btw+btw_shared_btw
btw_res_btw<-(1-(.69463-.40244)/.69463)

#Variation in momentary reactivity - Need to calculate (will include in text - no graph needed for this value)

#How much does DN account for variation in momentary reactions to events... (which definitely exist)

#Total contribution of DN = unique DN effect + reactivity effect (which is very small) + shared contextual variance (potential indirect pathway...)
#Can think of this as the pathways of DN to negative mood
DN_total<-btw_DN_tot+tot_var_react
DN_unique_DN<-btw_DN_uni_tot/DN_total
DN_indirect_DN<-btw_shared_tot/DN_total
DN_react_DN<-tot_var_react/DN_total

DF_tot_var<-data.frame(Value = c(btw_context_uni_tot,
                                 btw_shared_tot, 
                                 btw_DN_uni_tot, 
                                 btw_res_tot,
                                 wthn_best_worst_tot,
                                 wthn_res_tot, 
                                 wthn_best_worst_wthn, 
                                 wthn_res_wthn, 
                                 btw_context_uni_btw,
                                 btw_shared_btw, 
                                 btw_DN_uni_btw, 
                                 btw_res_btw
), 
Component = c(rep('Total', 6), 
              rep('Within-Subjects', 2), 
              rep('Between-Subjects', 4)), 
Term = c('Event Context Variance', 
         'DN-Event Context Shared Variance', 
         'DN Variance', 
         'Unmodeled Between-Subjects Variance',
         'Momentary Event Variance', 
         'Unmodeled Momentary Variance',
         'Momentary Event Variance', 
         'Unmodeled Momentary Variance',
         'Event Context Variance', 
         'DN-Event Context Shared Variance', 
         'DN Variance', 
         'Unmodeled Between-Subjects Variance'
)
)

#------------------------------------------------------------------------------
#Making riverplots (Sankey plots) for Positive Mood Model
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

palette = c(paste0(brewer.pal(9, "Blues"), 90)[5:8], 
            paste0(brewer.pal(9, "Reds"), 90)[7:8], 
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

png(paste0(study1.graphics, '/S1_PosMood_river.png'), 
    units = 'in', 
    res = 900, 
    height = 10, 
    width = 10)
riverplot(riv, 
          nodewidth = 3, 
          plot_area = .95)
title(ylab = 'Study 1 - Riverplot of Total Variance Decomposition Derived from Positive Mood Model')
dev.off()


#Getting separate DN only River plot decomposing the ways in which DN contributes to mood
DN_total<-btw_DN_tot+tot_var_react
DN_unique_DN<-btw_DN_uni_tot/DN_total
DN_indirect_DN<-btw_shared_tot/DN_total
DN_react_DN<-tot_var_react/DN_total

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

palette = c(paste0(brewer.pal(9, "Blues"), 90)[5:6], 
            paste0(brewer.pal(9, "Blues"), 90)[4], 
            paste0(brewer.pal(9, "Blues"), 90)[7])

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

png(paste0(study1.graphics, '/S1_PosMood_DN_river.png'), 
    units = 'in', 
    res = 900, 
    height = 10, 
    width = 10)
riverplot(riv, 
          nodewidth = 3, 
          plot_area = .95)
title(ylab = 'Study 1 - Riverplot of DN Effect Decomposition Derived from Positive Mood Model')
dev.off()

#=========================================================================================================
#Descriptive models - determining direct associations between DN and event ratings
#=========================================================================================================

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Worst Magnitude
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
log(mean(dat.study1$Worst, na.rm=TRUE))
Worst_lv2<-brms::brm_multiple(T.Worst~1+c.DN+
                                (1|ID), 
                              data = dat.study1.list, 
                              family = 'lognormal',
                              prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                        set_prior('normal(0.8605,2)', class='Intercept'), 
                                        set_prior('normal(0,2)', class='b')),
                              warmup = 3000, 
                              iter = 6000,
                              thin = 3, 
                              chains = 4,
                              control = list(adapt_delta=.90, 
                                             max_treedepth=10), 
                              save_model = paste0(stan.code, '/S1_worst_comb.stan'))

summary(Worst_lv2)

round(Worst_lv2$rhats[,1:3], 3)

png(paste0(study1.graphics, '/S1_Worst_lv2_traceplots.png'), 
    units='in', width = 8, height=14, res=300)
plot(Worst_lv2, pars = c("^b_", 'sigma', '^sd_'), N=4)
dev.off()

#plot recovery of original distribution
png(paste0(study1.graphics, '/S1_Worst_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(Worst_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = Worst_lv2, 
             out.folder = study1.model,
             file = paste0('S1_Worst_lv2_cross_Summary_', Sys.Date()),
             DF = dat.study1.list[[1]], 
             tot.pars = 13, 
             impute = TRUE)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Best Magnitude
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
mean(dat.study1$Best, na.rm=TRUE)
Best_lv2<-brms::brm_multiple(T.Best~1+c.DN+
                               (1|ID), 
                             data = dat.study1.list, 
                             family = 'normal',
                             prior = c(set_prior('student_t(3,0,1)', class='sd'), 
                                       set_prior('normal(3.44598,2)', class='Intercept'), 
                                       set_prior('normal(0,2)', class='b')),
                             warmup = 3000, 
                             iter = 6000,
                             thin = 3, 
                             chains = 4,
                             control = list(adapt_delta=.90, 
                                            max_treedepth=10), 
                             save_model = paste0(stan.code, '/S1_worst_comb.stan'))

summary(Best_lv2)

round(Best_lv2$rhats[,1:3], 3)

png(paste0(study1.graphics, '/S1_Best_lv2_traceplots.png'), 
    units='in', width = 8, height=14, res=300)
plot(Best_lv2, pars = c("^b_", 'sigma', '^sd_'), N=4)
dev.off()

#plot recovery of original distribution
png(paste0(study1.graphics, '/S1_Best_lv2_pp_dens.png'), 
    units='in', width = 8, height=8, res=300)
pp_check(Best_lv2, nsamples=100)
dev.off()

#Storing most important aspects of Summary
bayes.to.txt(model = Best_lv2, 
             out.folder = study1.model,
             file = paste0('S1_Best_lv2_Summary_', Sys.Date()),
             DF = dat.study1.list[[1]], 
             tot.pars = 13, 
             impute = TRUE)

#####################################################################################################
#####################################################################################################
#Generating Additional Figures for Manuscript:
#####################################################################################################
#####################################################################################################

DN.mean<-vector()
for(i in 1:length(dat.id[,1])){
  DN.mean[i]<-mean(dat.study1$c.DN[dat.study1$ID==dat.id$`unique(dat$subid)`[i]], na.rm=T)
}

DF.plot<-cbind(dat.id, DN.mean)

g1<-ggplot(data=DF.plot)+
  geom_histogram(aes(x=DN.mean, y=..density..), 
                 color=brewer.pal(5,'Dark2')[1], 
                 fill=brewer.pal(5,'Accent')[1], 
                 bins = 17)+
  stat_density(aes(x=DN.mean),
               geom='line', 
               color=brewer.pal(5,'Accent')[5], 
               lwd=1.5)+
  scale_x_continuous(limits = c(-2.5, 2.5))+
  xlab(expression(paste('DN composite Score (', italic('N'), ' = 127)')))+
  ylab('Density')+
  ggtitle('Study 1 - Dispositional Negativity Score Distribution')

png(paste0(study1.graphics, '/Supplemental_Figure1_DN_distribution.png'), 
    res=300, 
    units='in', 
    height=6, 
    width=6)
g1
dev.off()

png(paste0(study1.graphics, '/Supplemental_Figure2_Exploratory_Mood_Figure.png'), 
    res=300, 
    units='in', 
    width=6, 
    height = 6)
psych::pairs.panels(dat.study1[c('NegAff', 'PosAff', 'Worst', 'Best')])
dev.off()

#Testing area
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
beta_00<-fixef(Neg_ucm)[1,1]
#=-=-=-=-=-=-=-=-=-=-=-=-=-=
#First creating a data set that averages across imputed values:

#Attempting Variance Partionining for DN alone
within_cov<-c(21,23)                    #Columns with group-mean centered predictors
between_cov<-c(6)                       #Columns with between-subject predictors
random_cov<-c(21,23)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Neg_lv2_DN, pars = 'b_Intercept') 
DN = posterior_samples(Neg_lv2_DN, pars = 'b_c.DN') 
Worst = posterior_samples(Neg_lv2_DN, pars = 'b_c.Worst') 
Best = posterior_samples(Neg_lv2_DN, pars = 'b_c.Best') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Neg_lv2_DN, pars = 'sd_ID__Intercept')^2
Worst_var = posterior_samples(Neg_lv2_DN, pars = 'sd_ID__c.Worst')^2 
Best_var = posterior_samples(Neg_lv2_DN, pars = 'sd_ID__c.Best')^2

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Worst = posterior_samples(Neg_lv2_DN, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_lv2_DN, pars = 'sd_ID__c.Worst')*
  posterior_samples(Neg_lv2_DN, pars = 'cor_ID__Intercept__c.Worst')

cov_Int_Best = posterior_samples(Neg_lv2_DN, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_lv2_DN, pars = 'sd_ID__c.Best')*
  posterior_samples(Neg_lv2_DN, pars = 'cor_ID__Intercept__c.Best') 

cov_Worst_Best = posterior_samples(Neg_lv2_DN, pars = 'sd_ID__c.Worst')*
  posterior_samples(Neg_lv2_DN, pars = 'sd_ID__c.Best')*
  posterior_samples(Neg_lv2_DN, pars = 'cor_ID__c.Worst__c.Best')

#Getting level 1 error variance
sigma<-log(1+posterior_samples(Neg_lv2_DN, pars = 'sigma')/beta_00)

post_samples<-data.frame(Intercept,
                         DN, 
                         Worst, 
                         Best, 
                         Int_var, 
                         Worst_var, 
                         Best_var, 
                         cov_Int_Worst, 
                         cov_Int_Best, 
                         cov_Worst_Best, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'DN', 
                          'Worst', 
                          'Best', 
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

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  D<-sample(sampling_list[[i]], size = 1000, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Worst', 'Best')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'DN')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Worst', 'cov_Int_Best')], 
           post_samples[D[d], c('cov_Int_Worst', 'Worst_var', 'cov_Worst_Best')], 
           post_samples[D[d], c('cov_Int_Best', 'cov_Worst_Best', 'Best_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=3)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_DN<-r2MLM(data=dat.study1.list[[i]], 
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

#Statement should be true if correct values are extracted from data/model 
mean(DN_alone$between_var)+mean(DN_alone$within_var)==1

hist(DN_alone$between_All_tot)
mean(DN_alone$between_All_tot)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=
#Event only model (these are context effects on intercept)
#Attempting Variance Partionining for DN alone
within_cov<-c(21,23)                    #Columns with group-mean centered predictors
between_cov<-c(22, 24)                       #Columns with between-subject predictors
random_cov<-c(21,23)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Neg_lv2_Evnt, pars = 'b_Intercept') 
M.Worst = posterior_samples(Neg_lv2_Evnt, pars = 'b_mean.Worst')
M.Best = posterior_samples(Neg_lv2_Evnt, pars = 'b_mean.Best')
Worst = posterior_samples(Neg_lv2_Evnt, pars = 'b_c.Worst') 
Best = posterior_samples(Neg_lv2_Evnt, pars = 'b_c.Best') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Neg_lv2_Evnt, pars = 'sd_ID__Intercept')^2
Worst_var = posterior_samples(Neg_lv2_Evnt, pars = 'sd_ID__c.Worst')^2 
Best_var = posterior_samples(Neg_lv2_Evnt, pars = 'sd_ID__c.Best')^2

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Worst = posterior_samples(Neg_lv2_Evnt, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_lv2_Evnt, pars = 'sd_ID__c.Worst')*
  posterior_samples(Neg_lv2_Evnt, pars = 'cor_ID__Intercept__c.Worst')

cov_Int_Best = posterior_samples(Neg_lv2_Evnt, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_lv2_Evnt, pars = 'sd_ID__c.Best')*
  posterior_samples(Neg_lv2_Evnt, pars = 'cor_ID__Intercept__c.Best') 

cov_Worst_Best = posterior_samples(Neg_lv2_Evnt, pars = 'sd_ID__c.Worst')*
  posterior_samples(Neg_lv2_Evnt, pars = 'sd_ID__c.Best')*
  posterior_samples(Neg_lv2_Evnt, pars = 'cor_ID__c.Worst__c.Best')

#Getting level 1 error variance
sigma<-log(1+posterior_samples(Neg_lv2_Evnt, pars = 'sigma')/beta_00)
post_samples<-data.frame(Intercept,
                         M.Worst,
                         M.Best,
                         Worst, 
                         Best, 
                         Int_var, 
                         Worst_var, 
                         Best_var, 
                         cov_Int_Worst, 
                         cov_Int_Best, 
                         cov_Worst_Best, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'M.Worst',
                          'M.Best',
                          'Worst', 
                          'Best', 
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

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  D<-sample(sampling_list[[i]], size = 1000, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Worst', 'Best')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'M.Worst', 'M.Best')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Worst', 'cov_Int_Best')], 
           post_samples[D[d], c('cov_Int_Worst', 'Worst_var', 'cov_Worst_Best')], 
           post_samples[D[d], c('cov_Int_Best', 'cov_Worst_Best', 'Best_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=3)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_Evnt<-r2MLM(data=dat.study1.list[[i]], 
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
    Evnt_alone$between_var<-c(Evnt_alone$between_var, as.numeric(r2mlm_Evnt$Decompositions['fixed, between', 'total'])+
                              as.numeric(r2mlm_Evnt$Decompositions['mean variation', 'total']))
    Evnt_alone$within_var<-c(Evnt_alone$within_var, as.numeric(r2mlm_Evnt$Decompositions['sigma2', 'total'])+
                             as.numeric(r2mlm_Evnt$Decompositions['slope variation', 'total'])+
                             as.numeric(r2mlm_Evnt$Decompositions['fixed, within', 'total']))
    Evnt_alone$between_All_tot<-c(Evnt_alone$between_All_tot, as.numeric(r2mlm_Evnt$R2s['f2', 'total']))
    Evnt_alone$between_All_btw<-c(Evnt_alone$between_All_btw, as.numeric(r2mlm_Evnt$R2s['f2', 'between']))
    Evnt_alone$between_res_btw<-c(Evnt_alone$between_res_btw, as.numeric(r2mlm_Evnt$R2s['m', 'between']))
    Evnt_alone$within_fix_wthn<-c(Evnt_alone$within_fix_wthn, as.numeric(r2mlm_Evnt$Decompositions['fixed, within', 'within']))
    Evnt_alone$within_fix_tot<-c(Evnt_alone$within_fix_tot, as.numeric(r2mlm_Evnt$Decompositions['fixed, within', 'total']))
    Evnt_alone$within_slope_var_wthn<-c(Evnt_alone$within_slope_var_wthn, as.numeric(r2mlm_Evnt$Decompositions['slope variation', 'within']))
    Evnt_alone$within_res_wthn<-c(Evnt_alone$within_res_wthn, as.numeric(r2mlm_Evnt$Decompositions['sigma2', 'within']))
    Evnt_alone$within_unmod_tot<-c(Evnt_alone$within_unmod_tot, as.numeric(r2mlm_Evnt$Decompositions['sigma2', 'total'])+
                                   as.numeric(r2mlm_Evnt$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}
hist(Evnt_alone$between_var)
mean(Evnt_alone$between_var)
hist(Evnt_alone$within_var)
mean(Evnt_alone$within_var)

#Should be TRUE
mean(Evnt_alone$between_var)+mean(Evnt_alone$within_var)==1

#Total variance explained by level 2 effects in model:
mean(Evnt_alone$between_All_tot)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=
#Combined Level 2 Model (these are context effects on intercept)
#Attempting Variance Partionining for DN alone
within_cov<-c(21,23)                    #Columns with group-mean centered predictors
between_cov<-c(6, 22, 24)                       #Columns with between-subject predictors
random_cov<-c(21,23)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Neg_lv2_All, pars = 'b_Intercept') 
DN = posterior_samples(Neg_lv2_All, pars = 'b_c.DN')
M.Worst = posterior_samples(Neg_lv2_All, pars = 'b_mean.Worst')
M.Best = posterior_samples(Neg_lv2_All, pars = 'b_mean.Best')
Worst = posterior_samples(Neg_lv2_All, pars = 'b_c.Worst') 
Best = posterior_samples(Neg_lv2_All, pars = 'b_c.Best') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Neg_lv2_All, pars = 'sd_ID__Intercept')^2
Worst_var = posterior_samples(Neg_lv2_All, pars = 'sd_ID__c.Worst')^2 
Best_var = posterior_samples(Neg_lv2_All, pars = 'sd_ID__c.Best')^2

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Worst = posterior_samples(Neg_lv2_All, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_lv2_All, pars = 'sd_ID__c.Worst')*
  posterior_samples(Neg_lv2_All, pars = 'cor_ID__Intercept__c.Worst')

cov_Int_Best = posterior_samples(Neg_lv2_All, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_lv2_All, pars = 'sd_ID__c.Best')*
  posterior_samples(Neg_lv2_All, pars = 'cor_ID__Intercept__c.Best') 

cov_Worst_Best = posterior_samples(Neg_lv2_All, pars = 'sd_ID__c.Worst')*
  posterior_samples(Neg_lv2_All, pars = 'sd_ID__c.Best')*
  posterior_samples(Neg_lv2_All, pars = 'cor_ID__c.Worst__c.Best')

#Getting level 1 error variance
sigma<-log(1+posterior_samples(Neg_lv2_All, pars = 'sigma')/beta_00)
post_samples<-data.frame(Intercept,
                         DN,
                         M.Worst,
                         M.Best,
                         Worst, 
                         Best, 
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
                          'Best', 
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
Lv2_comb<-list(between_var=vector(), 
               within_var=vector(), 
               between_All_tot=vector(),
               between_All_btw=vector(),
               between_res_btw=vector(),
               within_fix_wthn=vector(),
               within_fix_tot=vector(),
               within_slope_var_wthn=vector(),
               within_res_wthn=vector(), 
               within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  D<-sample(sampling_list[[i]], size = 1000, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Worst', 'Best')]))       #Make sure the effects line up - in order of within_cov
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
    
    r2mlm_All<-r2MLM(data=dat.study1.list[[i]], 
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
    Lv2_comb$between_var<-c(Lv2_comb$between_var, as.numeric(r2mlm_All$Decompositions['fixed, between', 'total'])+
                                as.numeric(r2mlm_All$Decompositions['mean variation', 'total']))
    Lv2_comb$within_var<-c(Lv2_comb$within_var, as.numeric(r2mlm_All$Decompositions['sigma2', 'total'])+
                               as.numeric(r2mlm_All$Decompositions['slope variation', 'total'])+
                               as.numeric(r2mlm_All$Decompositions['fixed, within', 'total']))
    Lv2_comb$between_All_tot<-c(Lv2_comb$between_All_tot, as.numeric(r2mlm_All$R2s['f2', 'total']))
    Lv2_comb$between_All_btw<-c(Lv2_comb$between_All_btw, as.numeric(r2mlm_All$R2s['f2', 'between']))
    Lv2_comb$between_res_btw<-c(Lv2_comb$between_res_btw, as.numeric(r2mlm_All$R2s['m', 'between']))
    Lv2_comb$within_fix_wthn<-c(Lv2_comb$within_fix_wthn, as.numeric(r2mlm_All$Decompositions['fixed, within', 'within']))
    Lv2_comb$within_fix_tot<-c(Lv2_comb$within_fix_tot, as.numeric(r2mlm_All$Decompositions['fixed, within', 'total']))
    Lv2_comb$within_slope_var_wthn<-c(Lv2_comb$within_slope_var_wthn, as.numeric(r2mlm_All$Decompositions['slope variation', 'within']))
    Lv2_comb$within_res_wthn<-c(Lv2_comb$within_res_wthn, as.numeric(r2mlm_All$Decompositions['sigma2', 'within']))
    Lv2_comb$within_unmod_tot<-c(Lv2_comb$within_unmod_tot, as.numeric(r2mlm_All$Decompositions['sigma2', 'total'])+
                                     as.numeric(r2mlm_All$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}

hist(Lv2_comb$between_var)
mean(Lv2_comb$between_var)
hist(Lv2_comb$within_var)
mean(Lv2_comb$within_var)

#Should be TRUE
mean(Lv2_comb$between_var)+mean(Lv2_comb$within_var)==1

#Total variance explained by level 2 effects in model:
mean(Lv2_comb$between_All_tot)

#=-=-=-=-=-=-=-=-=-=
#Cross-Level Model 
within_cov<-c(21, 23, 25, 26)           #Columns with group-mean centered predictors & cross-level interactions
between_cov<-c(6, 22, 24)               #Columns with between-subject predictors
                                        #Note I have added a series of cross level interactions
                                        #Variables for these cross-level interactions are calculated in imputed data sets
random_cov<-c(21,23)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Neg_cross, pars = 'b_Intercept') 
DN = posterior_samples(Neg_cross, pars = 'b_c.DN')
M.Worst = posterior_samples(Neg_cross, pars = 'b_mean.Worst')
M.Best = posterior_samples(Neg_cross, pars = 'b_mean.Best')
Worst = posterior_samples(Neg_cross, pars = 'b_c.Worst') 
Best = posterior_samples(Neg_cross, pars = 'b_c.Best')
DNxWorst = posterior_samples(Neg_cross, pars = 'b_c.Worst:c.DN')
DNxBest = posterior_samples(Neg_cross, pars = 'b_c.Best:c.DN')

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Neg_cross, pars = 'sd_ID__Intercept')^2
Worst_var = posterior_samples(Neg_cross, pars = 'sd_ID__c.Worst')^2 
Best_var = posterior_samples(Neg_cross, pars = 'sd_ID__c.Best')^2

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Worst = posterior_samples(Neg_cross, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_cross, pars = 'sd_ID__c.Worst')*
  posterior_samples(Neg_cross, pars = 'cor_ID__Intercept__c.Worst')

cov_Int_Best = posterior_samples(Neg_cross, pars = 'sd_ID__Intercept')*
  posterior_samples(Neg_cross, pars = 'sd_ID__c.Best')*
  posterior_samples(Neg_cross, pars = 'cor_ID__Intercept__c.Best') 

cov_Worst_Best = posterior_samples(Neg_cross, pars = 'sd_ID__c.Worst')*
  posterior_samples(Neg_cross, pars = 'sd_ID__c.Best')*
  posterior_samples(Neg_cross, pars = 'cor_ID__c.Worst__c.Best')

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

for(i in 1:length(dat.study1.list)){
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
    
    dat.study1.list[[i]]$DNxWorst<-dat.study1.list[[i]]$c.DN*dat.study1.list[[i]]$c.Worst
    dat.study1.list[[i]]$DNxBest<-dat.study1.list[[i]]$c.DN*dat.study1.list[[i]]$c.Best
    
    r2mlm_All<-r2MLM(data=dat.study1.list[[i]], 
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
    Lv2_cross$between_var<-c(Lv2_cross$between_var, as.numeric(r2mlm_All$Decompositions['fixed, between', 'total'])+
                              as.numeric(r2mlm_All$Decompositions['mean variation', 'total']))
    Lv2_cross$within_var<-c(Lv2_cross$within_var, as.numeric(r2mlm_All$Decompositions['sigma2', 'total'])+
                             as.numeric(r2mlm_All$Decompositions['slope variation', 'total'])+
                             as.numeric(r2mlm_All$Decompositions['fixed, within', 'total']))
    Lv2_cross$between_All_tot<-c(Lv2_cross$between_All_tot, as.numeric(r2mlm_All$R2s['f2', 'total']))
    Lv2_cross$between_All_btw<-c(Lv2_cross$between_All_btw, as.numeric(r2mlm_All$R2s['f2', 'between']))
    Lv2_cross$between_res_btw<-c(Lv2_cross$between_res_btw, as.numeric(r2mlm_All$R2s['m', 'between']))
    Lv2_cross$within_fix_wthn<-c(Lv2_cross$within_fix_wthn, as.numeric(r2mlm_All$Decompositions['fixed, within', 'within']))
    Lv2_cross$within_fix_tot<-c(Lv2_cross$within_fix_tot, as.numeric(r2mlm_All$Decompositions['fixed, within', 'total']))
    Lv2_cross$within_slope_var_wthn<-c(Lv2_cross$within_slope_var_wthn, as.numeric(r2mlm_All$Decompositions['slope variation', 'within']))
    Lv2_cross$within_res_wthn<-c(Lv2_cross$within_res_wthn, as.numeric(r2mlm_All$Decompositions['sigma2', 'within']))
    Lv2_cross$within_unmod_tot<-c(Lv2_cross$within_unmod_tot, as.numeric(r2mlm_All$Decompositions['sigma2', 'total'])+
                                   as.numeric(r2mlm_All$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}

#To View Within and Between Variation in Event model (for comparison)
hist(DN_alone$between_var)
mean(DN_alone$between_var)
hist(DN_alone$within_var)
mean(DN_alone$within_var)

#To View Within and Between Variation in Event model (for comparison)
hist(Evnt_alone$between_var)
mean(Evnt_alone$between_var)
hist(Evnt_alone$within_var)
mean(Evnt_alone$within_var)

#To View Within and Between Variation in Combined model (for comparison)
hist(Lv2_comb$between_var)
mean(Lv2_comb$between_var)
hist(Lv2_comb$within_var)
mean(Lv2_comb$within_var)

#Examining Within and Between Variation in Final Cross Model
hist(Lv2_cross$between_var)
mean(Lv2_cross$between_var)
hist(Lv2_cross$within_var)
mean(Lv2_cross$within_var)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=
#Using Combined model as reference point for variance 
#Slight differences in final estimates can be achieved with different model combinations 
#This is because there are slight, non-meaningful differences in the total between & within variances from model to model
btw_DN_uni_tot<-mean(Lv2_comb$between_All_tot)-mean(Evnt_alone$between_All_tot)
btw_shared_tot<-mean(DN_alone$between_All_tot)-btw_DN_uni_tot
btw_context_uni_tot<-mean(Lv2_comb$between_All_tot)-mean(DN_alone$between_All_tot)

sum(c(btw_DN_uni_tot, 
      btw_shared_tot,
      btw_context_uni_tot))

#The above total should equal: 
mean(Lv2_comb$between_All_tot)

#Now getting the between subjects residual variance from final model
btw_res_tot<-mean(Lv2_comb$between_var-Lv2_comb$between_All_tot)
btw_res_tot+sum(c(btw_DN_uni_tot, 
                  btw_shared_tot,
                  btw_context_uni_tot))
#Final getting the within subjects variance: 
wthn_best_worst_tot<-mean(Lv2_comb$within_fix_wthn)
wthn_res_tot<-mean(Lv2_comb$within_var)-wthn_best_worst_tot

#Should be equal to 1
sum(c(wthn_best_worst_tot, 
      wthn_res_tot, 
      mean(Lv2_comb$between_var)))

tot_btw<-mean(Lv2_comb$between_var)
tot_wthn<-mean(Lv2_comb$within_var)  

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

png(paste0(study1.graphics, '/S1_NegMood_river.png'), 
    units = 'in', 
    res = 1200, 
    height = 10, 
    width = 10)
riverplot(riv, 
          nodewidth = 3, 
          plot_area = .95)
title(ylab = 'Study 1 - Riverplot of Total Variance Decomposition Estimated from Negative Mood Models')
dev.off()

#-#-#-#-#-#-#-#-#-#
#Obtaining DN-specific plot
tot_var_react<-mean(Lv2_cross$within_fix_tot) - mean(Lv2_comb$within_fix_tot) 
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

png(paste0(study1.graphics, '/S1_NegMood_DN_River.png'), 
    units = 'in', 
    res = 900, 
    height = 10, 
    width = 10)
riverplot(riv, 
          nodewidth = 3, 
          plot_area = .95)
title(ylab = 'Study 1 - Riverplot of DN Effect Decomposition Derived from Negative Mood Model')
dev.off()

postscript(paste0(study1.graphics, '/test_riverDN.eps'), 
           width = 10, 
           height = 10, 
           horizontal = FALSE, 
           onefile = FALSE)
riverplot(riv, 
          nodewidth = 3, 
          plot_area = .95)
title(ylab = 'Study 1 - Riverplot of DN Effect Decomposition Derived from Negative Mood Model')
dev.off()

#==============================================================================
#POSITIVE MOOD MODELS - Variance Decomposition
#==============================================================================

#Attempting Variance Partionining for DN alone
within_cov<-c(21,23)                    #Columns with group-mean centered predictors
between_cov<-c(6)                       #Columns with between-subject predictors
random_cov<-c(21,23)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Pos_lv2_DN, pars = 'b_Intercept') 
DN = posterior_samples(Pos_lv2_DN, pars = 'b_c.DN') 
Worst = posterior_samples(Pos_lv2_DN, pars = 'b_c.Worst') 
Best = posterior_samples(Pos_lv2_DN, pars = 'b_c.Best') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Pos_lv2_DN, pars = 'sd_ID__Intercept')^2
Worst_var = posterior_samples(Pos_lv2_DN, pars = 'sd_ID__c.Worst')^2 
Best_var = posterior_samples(Pos_lv2_DN, pars = 'sd_ID__c.Best')^2

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Worst = posterior_samples(Pos_lv2_DN, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_lv2_DN, pars = 'sd_ID__c.Worst')*
  posterior_samples(Pos_lv2_DN, pars = 'cor_ID__Intercept__c.Worst')

cov_Int_Best = posterior_samples(Pos_lv2_DN, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_lv2_DN, pars = 'sd_ID__c.Best')*
  posterior_samples(Pos_lv2_DN, pars = 'cor_ID__Intercept__c.Best') 

cov_Worst_Best = posterior_samples(Pos_lv2_DN, pars = 'sd_ID__c.Worst')*
  posterior_samples(Pos_lv2_DN, pars = 'sd_ID__c.Best')*
  posterior_samples(Pos_lv2_DN, pars = 'cor_ID__c.Worst__c.Best')

#Getting level 1 error variance
sigma<-posterior_samples(Pos_lv2_DN, pars = 'sigma')

post_samples<-data.frame(Intercept,
                         DN, 
                         Worst, 
                         Best, 
                         Int_var, 
                         Worst_var, 
                         Best_var, 
                         cov_Int_Worst, 
                         cov_Int_Best, 
                         cov_Worst_Best, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'DN', 
                          'Worst', 
                          'Best', 
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

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  D<-sample(sampling_list[[i]], size = 1000, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Worst', 'Best')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'DN')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Worst', 'cov_Int_Best')], 
           post_samples[D[d], c('cov_Int_Worst', 'Worst_var', 'cov_Worst_Best')], 
           post_samples[D[d], c('cov_Int_Best', 'cov_Worst_Best', 'Best_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=3)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_DN<-r2MLM(data=dat.study1.list[[i]], 
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
    P_DN_alone$within_slope_var_wthn<-c(P_DN_alone$within_slope_wthn, as.numeric(r2mlm_DN$Decompositions['slope variation', 'within']))
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

#Statement should be true if correct values are extracted from data/model 
mean(P_DN_alone$between_var)+mean(P_DN_alone$within_var)==1

hist(P_DN_alone$between_All_tot)
mean(P_DN_alone$between_All_tot)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=
#Event only model (these are context effects on intercept)
#Attempting Variance Partionining for DN alone
within_cov<-c(21,23)                    #Columns with group-mean centered predictors
between_cov<-c(22, 24)                       #Columns with between-subject predictors
random_cov<-c(21,23)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Pos_lv2_Evnt, pars = 'b_Intercept') 
M.Worst = posterior_samples(Pos_lv2_Evnt, pars = 'b_mean.Worst')
M.Best = posterior_samples(Pos_lv2_Evnt, pars = 'b_mean.Best')
Worst = posterior_samples(Pos_lv2_Evnt, pars = 'b_c.Worst') 
Best = posterior_samples(Pos_lv2_Evnt, pars = 'b_c.Best') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Pos_lv2_Evnt, pars = 'sd_ID__Intercept')^2
Worst_var = posterior_samples(Pos_lv2_Evnt, pars = 'sd_ID__c.Worst')^2 
Best_var = posterior_samples(Pos_lv2_Evnt, pars = 'sd_ID__c.Best')^2

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Worst = posterior_samples(Pos_lv2_Evnt, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_lv2_Evnt, pars = 'sd_ID__c.Worst')*
  posterior_samples(Pos_lv2_Evnt, pars = 'cor_ID__Intercept__c.Worst')

cov_Int_Best = posterior_samples(Pos_lv2_Evnt, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_lv2_Evnt, pars = 'sd_ID__c.Best')*
  posterior_samples(Pos_lv2_Evnt, pars = 'cor_ID__Intercept__c.Best') 

cov_Worst_Best = posterior_samples(Pos_lv2_Evnt, pars = 'sd_ID__c.Worst')*
  posterior_samples(Pos_lv2_Evnt, pars = 'sd_ID__c.Best')*
  posterior_samples(Pos_lv2_Evnt, pars = 'cor_ID__c.Worst__c.Best')

#Getting level 1 error variance
sigma<-posterior_samples(Pos_lv2_Evnt, pars = 'sigma')
post_samples<-data.frame(Intercept,
                         M.Worst,
                         M.Best,
                         Worst, 
                         Best, 
                         Int_var, 
                         Worst_var, 
                         Best_var, 
                         cov_Int_Worst, 
                         cov_Int_Best, 
                         cov_Worst_Best, 
                         sigma)

colnames(post_samples)<-c('Intercept', 
                          'M.Worst',
                          'M.Best',
                          'Worst', 
                          'Best', 
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

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  D<-sample(sampling_list[[i]], size = 1000, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Worst', 'Best')]))       #Make sure the effects line up - in order of within_cov
    Gamma_b<-as.vector(
      as.matrix(
        post_samples[D[d], c('Intercept', 'M.Worst', 'M.Best')]))    #level fixed intercept and any level 2 fixed effects
    tau<-c(post_samples[D[d], c('Int_var', 'cov_Int_Worst', 'cov_Int_Best')], 
           post_samples[D[d], c('cov_Int_Worst', 'Worst_var', 'cov_Worst_Best')], 
           post_samples[D[d], c('cov_Int_Best', 'cov_Worst_Best', 'Best_var')])
    
    tau<-matrix(unlist(tau), 
                byrow=TRUE, 
                ncol=3)
    
    #Matrix rows/columns ordered starting with tau.00 (intercept variance)
    #then add columns/rows for random effects of slopes in order of within_cov
    
    Sigma2<-post_samples[D[d], 'sigma']
    
    r2mlm_Evnt<-r2MLM(data=dat.study1.list[[i]], 
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
    P_Evnt_alone$between_var<-c(P_Evnt_alone$between_var, as.numeric(r2mlm_Evnt$Decompositions['fixed, between', 'total'])+
                                as.numeric(r2mlm_Evnt$Decompositions['mean variation', 'total']))
    P_Evnt_alone$within_var<-c(P_Evnt_alone$within_var, as.numeric(r2mlm_Evnt$Decompositions['sigma2', 'total'])+
                               as.numeric(r2mlm_Evnt$Decompositions['slope variation', 'total'])+
                               as.numeric(r2mlm_Evnt$Decompositions['fixed, within', 'total']))
    P_Evnt_alone$between_All_tot<-c(P_Evnt_alone$between_All_tot, as.numeric(r2mlm_Evnt$R2s['f2', 'total']))
    P_Evnt_alone$between_All_btw<-c(P_Evnt_alone$between_All_btw, as.numeric(r2mlm_Evnt$R2s['f2', 'between']))
    P_Evnt_alone$between_res_btw<-c(P_Evnt_alone$between_res_btw, as.numeric(r2mlm_Evnt$R2s['m', 'between']))
    P_Evnt_alone$within_fix_wthn<-c(P_Evnt_alone$within_fix_wthn, as.numeric(r2mlm_Evnt$Decompositions['fixed, within', 'within']))
    P_Evnt_alone$within_fix_tot<-c(P_Evnt_alone$within_fix_tot, as.numeric(r2mlm_Evnt$Decompositions['fixed, within', 'total']))
    P_Evnt_alone$within_slope_var_wthn<-c(P_Evnt_alone$within_slope_var_wthn, as.numeric(r2mlm_Evnt$Decompositions['slope variation', 'within']))
    P_Evnt_alone$within_res_wthn<-c(P_Evnt_alone$within_res_wthn, as.numeric(r2mlm_Evnt$Decompositions['sigma2', 'within']))
    P_Evnt_alone$within_unmod_tot<-c(P_Evnt_alone$within_unmod_tot, as.numeric(r2mlm_Evnt$Decompositions['sigma2', 'total'])+
                                     as.numeric(r2mlm_Evnt$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}
hist(P_Evnt_alone$between_var)
mean(P_Evnt_alone$between_var)
hist(P_Evnt_alone$within_var)
mean(P_Evnt_alone$within_var)

#Should be TRUE
mean(P_Evnt_alone$between_var)+mean(P_Evnt_alone$within_var)==1

#Total variance explained by level 2 effects in model:
mean(P_Evnt_alone$between_All_tot)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=
#Combined Level 2 Model (these are context effects on intercept)
#Attempting Variance Partionining for DN alone
within_cov<-c(21,23)                    #Columns with group-mean centered predictors
between_cov<-c(6, 22, 24)                       #Columns with between-subject predictors
random_cov<-c(21,23)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Pos_lv2_All, pars = 'b_Intercept') 
DN = posterior_samples(Pos_lv2_All, pars = 'b_c.DN')
M.Worst = posterior_samples(Pos_lv2_All, pars = 'b_mean.Worst')
M.Best = posterior_samples(Pos_lv2_All, pars = 'b_mean.Best')
Worst = posterior_samples(Pos_lv2_All, pars = 'b_c.Worst') 
Best = posterior_samples(Pos_lv2_All, pars = 'b_c.Best') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Pos_lv2_All, pars = 'sd_ID__Intercept')^2
Worst_var = posterior_samples(Pos_lv2_All, pars = 'sd_ID__c.Worst')^2 
Best_var = posterior_samples(Pos_lv2_All, pars = 'sd_ID__c.Best')^2

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Worst = posterior_samples(Pos_lv2_All, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_lv2_All, pars = 'sd_ID__c.Worst')*
  posterior_samples(Pos_lv2_All, pars = 'cor_ID__Intercept__c.Worst')

cov_Int_Best = posterior_samples(Pos_lv2_All, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_lv2_All, pars = 'sd_ID__c.Best')*
  posterior_samples(Pos_lv2_All, pars = 'cor_ID__Intercept__c.Best') 

cov_Worst_Best = posterior_samples(Pos_lv2_All, pars = 'sd_ID__c.Worst')*
  posterior_samples(Pos_lv2_All, pars = 'sd_ID__c.Best')*
  posterior_samples(Pos_lv2_All, pars = 'cor_ID__c.Worst__c.Best')

#Getting level 1 error variance
sigma<-posterior_samples(Pos_lv2_All, pars = 'sigma')
post_samples<-data.frame(Intercept,
                         DN,
                         M.Worst,
                         M.Best,
                         Worst, 
                         Best, 
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
                          'Best', 
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
P_Lv2_comb<-list(between_var=vector(), 
               within_var=vector(), 
               between_All_tot=vector(),
               between_All_btw=vector(),
               between_res_btw=vector(),
               within_fix_wthn=vector(),
               within_fix_tot=vector(),
               within_slope_var_wthn=vector(),
               within_res_wthn=vector(), 
               within_unmod_tot=vector())

for(i in 1:length(dat.study1.list)){
  print(paste('Starting w/ imputed set', i, 'out of', M))
  D<-sample(sampling_list[[i]], size = 1000, replace = FALSE)
  for(d in 1:length(D)){
    #browser()
    Gamma_w<-as.vector(
      as.matrix(
        post_samples[D[d], c('Worst', 'Best')]))       #Make sure the effects line up - in order of within_cov
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
    
    r2mlm_All<-r2MLM(data=dat.study1.list[[i]], 
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
    P_Lv2_comb$between_var<-c(P_Lv2_comb$between_var, as.numeric(r2mlm_All$Decompositions['fixed, between', 'total'])+
                              as.numeric(r2mlm_All$Decompositions['mean variation', 'total']))
    P_Lv2_comb$within_var<-c(P_Lv2_comb$within_var, as.numeric(r2mlm_All$Decompositions['sigma2', 'total'])+
                             as.numeric(r2mlm_All$Decompositions['slope variation', 'total'])+
                             as.numeric(r2mlm_All$Decompositions['fixed, within', 'total']))
    P_Lv2_comb$between_All_tot<-c(P_Lv2_comb$between_All_tot, as.numeric(r2mlm_All$R2s['f2', 'total']))
    P_Lv2_comb$between_All_btw<-c(P_Lv2_comb$between_All_btw, as.numeric(r2mlm_All$R2s['f2', 'between']))
    P_Lv2_comb$between_res_btw<-c(P_Lv2_comb$between_res_btw, as.numeric(r2mlm_All$R2s['m', 'between']))
    P_Lv2_comb$within_fix_wthn<-c(P_Lv2_comb$within_fix_wthn, as.numeric(r2mlm_All$Decompositions['fixed, within', 'within']))
    P_Lv2_comb$within_fix_tot<-c(P_Lv2_comb$within_fix_tot, as.numeric(r2mlm_All$Decompositions['fixed, within', 'total']))
    P_Lv2_comb$within_slope_var_wthn<-c(P_Lv2_comb$within_slope_var_wthn, as.numeric(r2mlm_All$Decompositions['slope variation', 'within']))
    P_Lv2_comb$within_res_wthn<-c(P_Lv2_comb$within_res_wthn, as.numeric(r2mlm_All$Decompositions['sigma2', 'within']))
    P_Lv2_comb$within_unmod_tot<-c(P_Lv2_comb$within_unmod_tot, as.numeric(r2mlm_All$Decompositions['sigma2', 'total'])+
                                   as.numeric(r2mlm_All$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}

hist(P_Lv2_comb$between_var)
mean(P_Lv2_comb$between_var)
hist(P_Lv2_comb$within_var)
mean(P_Lv2_comb$within_var)

#Should be TRUE
mean(P_Lv2_comb$between_var)+mean(P_Lv2_comb$within_var)==1

#Total variance explained by level 2 effects in model:
mean(P_Lv2_comb$between_All_tot)

#=-=-=-=-=-=-=-=-=-=
#Cross-Level Model 
within_cov<-c(21, 23, 25, 26)           #Columns with group-mean centered predictors & cross-level interactions
between_cov<-c(6, 22, 24)               #Columns with between-subject predictors
#Note I have added a series of cross level interactions
#Variables for these cross-level interactions are calculated in imputed data sets
random_cov<-c(21,23)                    #Level 1 predictors with random slopes

#Getting distributions of fixed effects
Intercept = posterior_samples(Pos_cross, pars = 'b_Intercept') 
DN = posterior_samples(Pos_cross, pars = 'b_c.DN')
M.Worst = posterior_samples(Pos_cross, pars = 'b_mean.Worst')
M.Best = posterior_samples(Pos_cross, pars = 'b_mean.Best')
Worst = posterior_samples(Pos_cross, pars = 'b_c.Worst') 
Best = posterior_samples(Pos_cross, pars = 'b_c.Best') 
DNxWorst = posterior_samples(Pos_cross, pars = 'b_c.Worst:c.DN') 
DNxBest = posterior_samples(Pos_cross, pars = 'b_c.Best:c.DN') 

#Getting between-subjects variance terms for random effects
Int_var = posterior_samples(Pos_cross, pars = 'sd_ID__Intercept')^2
Worst_var = posterior_samples(Pos_cross, pars = 'sd_ID__c.Worst')^2 
Best_var = posterior_samples(Pos_cross, pars = 'sd_ID__c.Best')^2

#Getting covariances of random effects (i.e., tau matrix)
cov_Int_Worst = posterior_samples(Pos_cross, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_cross, pars = 'sd_ID__c.Worst')*
  posterior_samples(Pos_cross, pars = 'cor_ID__Intercept__c.Worst')

cov_Int_Best = posterior_samples(Pos_cross, pars = 'sd_ID__Intercept')*
  posterior_samples(Pos_cross, pars = 'sd_ID__c.Best')*
  posterior_samples(Pos_cross, pars = 'cor_ID__Intercept__c.Best') 

cov_Worst_Best = posterior_samples(Pos_cross, pars = 'sd_ID__c.Worst')*
  posterior_samples(Pos_cross, pars = 'sd_ID__c.Best')*
  posterior_samples(Pos_cross, pars = 'cor_ID__c.Worst__c.Best')

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

for(i in 1:length(dat.study1.list)){
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
    
    dat.study1.list[[i]]$DNxWorst<-dat.study1.list[[i]]$c.DN*dat.study1.list[[i]]$c.Worst
    dat.study1.list[[i]]$DNxBest<-dat.study1.list[[i]]$c.DN*dat.study1.list[[i]]$c.Best
    
    r2mlm_All<-r2MLM(data=dat.study1.list[[i]], 
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
    P_Lv2_cross$between_var<-c(P_Lv2_cross$between_var, as.numeric(r2mlm_All$Decompositions['fixed, between', 'total'])+
                               as.numeric(r2mlm_All$Decompositions['mean variation', 'total']))
    P_Lv2_cross$within_var<-c(P_Lv2_cross$within_var, as.numeric(r2mlm_All$Decompositions['sigma2', 'total'])+
                              as.numeric(r2mlm_All$Decompositions['slope variation', 'total'])+
                              as.numeric(r2mlm_All$Decompositions['fixed, within', 'total']))
    P_Lv2_cross$between_All_tot<-c(P_Lv2_cross$between_All_tot, as.numeric(r2mlm_All$R2s['f2', 'total']))
    P_Lv2_cross$between_All_btw<-c(P_Lv2_cross$between_All_btw, as.numeric(r2mlm_All$R2s['f2', 'between']))
    P_Lv2_cross$between_res_btw<-c(P_Lv2_cross$between_res_btw, as.numeric(r2mlm_All$R2s['m', 'between']))
    P_Lv2_cross$within_fix_wthn<-c(P_Lv2_cross$within_fix_wthn, as.numeric(r2mlm_All$Decompositions['fixed, within', 'within']))
    P_Lv2_cross$within_fix_tot<-c(P_Lv2_cross$within_fix_tot, as.numeric(r2mlm_All$Decompositions['fixed, within', 'total']))
    P_Lv2_cross$within_slope_var_wthn<-c(P_Lv2_cross$within_slope_var_wthn, as.numeric(r2mlm_All$Decompositions['slope variation', 'within']))
    P_Lv2_cross$within_res_wthn<-c(P_Lv2_cross$within_res_wthn, as.numeric(r2mlm_All$Decompositions['sigma2', 'within']))
    P_Lv2_cross$within_unmod_tot<-c(P_Lv2_cross$within_unmod_tot, as.numeric(r2mlm_All$Decompositions['sigma2', 'total'])+
                                    as.numeric(r2mlm_All$Decompositions['slope variation', 'total']))
    
    print(paste("Sampling", d, "out of", length(D), '- Dataset:', i))
  }
  print(paste("*** Finshed imputed data set", i, 'out of', M, "***"))
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=
#Using Combined model as reference point for variance 
#Slight differences in final estimates can be achieved with different model combinations 
#This is because there are slight, non-meaningful differences in the total between & within variances from model to model
P_btw_DN_uni_tot<-mean(P_Lv2_comb$between_All_tot)-mean(P_Evnt_alone$between_All_tot)
P_btw_shared_tot<-mean(P_DN_alone$between_All_tot)-P_btw_DN_uni_tot
P_btw_context_uni_tot<-mean(P_Lv2_comb$between_All_tot)-mean(P_DN_alone$between_All_tot)

sum(c(P_btw_DN_uni_tot, 
      P_btw_shared_tot,
      P_btw_context_uni_tot))

#The above total should equal: 
mean(P_Lv2_comb$between_All_tot)

#Now getting the between subjects residual variance from final model
P_btw_res_tot<-mean(P_Lv2_comb$between_var-P_Lv2_comb$between_All_tot)
P_btw_res_tot+sum(c(P_btw_DN_uni_tot, 
                    P_btw_shared_tot,
                    P_btw_context_uni_tot))
#Final getting the within subjects variance: 
P_wthn_best_worst_tot<-mean(P_Lv2_comb$within_fix_tot)
P_wthn_res_tot<-mean(P_Lv2_comb$within_var)-P_wthn_best_worst_tot

#Should be equal to 1
sum(c(P_wthn_best_worst_tot, 
      P_wthn_res_tot, 
      mean(P_Lv2_comb$between_var)))

P_tot_btw<-mean(P_Lv2_comb$between_var)
P_tot_wthn<-mean(P_Lv2_comb$within_var)  

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

png(paste0(study1.graphics, '/S1_PosMood_river.png'), 
    units = 'in', 
    res = 1200, 
    height = 10, 
    width = 10)
riverplot(riv, 
          nodewidth = 3, 
          plot_area = .95)
title(ylab = 'Study 1 - Riverplot of Total Variance Decomposition Estimated from Positive Mood Models')
dev.off()

#--
P_tot_var_react<-mean(P_Lv2_cross$within_fix_tot) - mean(P_Lv2_comb$within_fix_tot) 
P_DN_total <- P_btw_DN_uni_tot + P_btw_shared_tot + P_tot_var_react 
P_DN_unique_DN<-P_btw_DN_uni_tot/P_DN_total
P_DN_indirect_DN<-P_btw_shared_tot/P_DN_total
P_DN_react_DN<-P_tot_var_react/P_DN_total

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

png(paste0(study1.graphics, '/S1_PosMood_DN_River.png'), 
    units = 'in', 
    res = 900, 
    height = 10, 
    width = 10)
riverplot(riv, 
          nodewidth = 3, 
          plot_area = .95)
title(ylab = 'Study 1 - Riverplot of DN Effect Decomposition Derived from Positive Mood Model')
dev.off()
