####################################################################################################
# Study 2 Modeling Script: Generating Imputation Data Set

# Description: 
#   Most of this script involves a data preparation and cleaning. First summary values were 
#   were generated for the relevant EMA variables and these summary variables ultimately formed the 
#   basis for generating posterior distributions for each individual using a series of multilevel
#   imputation models. 

# Modeling Notes: 
#   Efforts to re-create these analyses are best done on a machine with sufficient memory and the 
#   user would likely benefit from running in an independent terminal if s/he wants to be able to 
#   continue separate lines of development/coding work. 
####################################################################################################

#---------------------------------------------------------------------------------------------------
# Package import
library(pan)
library(mitml)
library(mice)
library(tidyverse)
library(corrgram)
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------
# Folder setup
user<-ifelse(Sys.getenv("USERPROFILE")=="", "~", Sys.getenv("USERPROFILE"))
wd<-paste0(user, '/Dropbox/UMD/Shackman Lab/EMA_MS')
data.folder<-paste0(wd, '/Data')
study2.out<-paste0(wd, '/Study 2 output')
study2.graphics<-paste0(study2.out, '/Graphics')
study2.model<-paste0(study2.out, '/Model summaries')
stan.code<-paste0(wd, '/Stan_code')
#---------------------------------------------------------------------------------------------------

dat.study2_lv1 <- read.csv(paste0(data.folder, "/PAX_T1_all_clean.csv"), 
                           header = TRUE, 
                           stringsAsFactors = FALSE)

head(dat.study2_lv1)
table(dat.study2_lv1$SIG_Corrected)
IDs <- unique(dat.study2_lv1$ID)

tmp.DF <- data.frame(ID = rep(IDs, each = 56), 
                     SIG_Corrected = rep(101:108, 7*length(IDs)), 
                     DAY_Corrected = rep(rep(1:7, each = 8))) 

# DAY_Corrected is the day within the EMA window (1-7) 
# SIG_Corrected is the survey signal value within a given day (per software output range from 
# 101-108)
lv1_tmp <- merge(tmp.DF, dat.study2_lv1, by = c("ID", "DAY_Corrected", "SIG_Corrected"), all = TRUE)

#---------------------------------------------------------------------------------------------------
# Creating Composites 
psych::describe(lv1_tmp)  # View basic descriptives for inspection
str(lv1_tmp)
lv1_tmp$ID <- as.character(lv1_tmp$ID)
table(lv1_tmp$ID)

lv1_tmp$NegEvnt <- ifelse(lv1_tmp$EMA_13_NegEvent == 0, "No", "Yes")
lv1_tmp$PosEvnt <- ifelse(lv1_tmp$EMA_12_PosEvent == 0, "No", "Yes")

dat.study2_lv2 <- foreign::read.spss(paste0(data.folder, "/PAX_lv2.sav"), 
                                     to.data.frame = TRUE, 
                                     stringsAsFactors = FALSE)

colnames(dat.study2_lv2)[1] <- "ID"
psych::describe(dat.study2_lv2)
str(dat.study2_lv2)
dat.study2_lv2$ID <- trimws(dat.study2_lv2$ID)

dat_comb_tmp <- merge(lv1_tmp, dat.study2_lv2, by = "ID")  

dat.study2_model <- dat_comb_tmp %>% 
  select(ID, DAY_Corrected, SIG_Corrected, 
         EMA_2a_Enthus, EMA_2b_Joy, EMA_2c_Cheer, 
         EMA_2d_Calm, EMA_2e_Cont, EMA_2f_Relax,
         EMA_2g_Nerv, EMA_2h_Worry, EMA_2i_Afraid, 
         EMA_2j_Annoy, EMA_2k_Angry, EMA_2l_Slug, 
         EMA_2m_Sad, EMA_2n_Tired, EMA_2o_Hopeless,
         PosEvnt, NegEvnt, SexNum, DN_ScreenAndLab1)

# High Positive Vars: 
# Adding 1 to place on variables on a 1-5 scale - relevant for log-transforms later
dat.study2_model$EMA_2a_Enthus <- dat.study2_model$EMA_2a_Enthus + 1
dat.study2_model$EMA_2b_Joy <- dat.study2_model$EMA_2b_Joy + 1
dat.study2_model$EMA_2c_Cheer <- dat.study2_model$EMA_2c_Cheer + 1
dat.study2_model$JOY <- rowMeans(dat.study2_model[c('EMA_2a_Enthus', 
                                                    'EMA_2b_Joy', 
                                                    'EMA_2c_Cheer')], 
                                 na.rm = TRUE)

# Low Positive Vars:
dat.study2_model$EMA_2d_Calm <- dat.study2_model$EMA_2d_Calm + 1
dat.study2_model$EMA_2e_Cont <- dat.study2_model$EMA_2e_Cont + 1
dat.study2_model$EMA_2f_Relax <- dat.study2_model$EMA_2f_Relax + 1
dat.study2_model$CALM <- rowMeans(dat.study2_model[c('EMA_2d_Calm', 
                                                     'EMA_2e_Cont', 
                                                     'EMA_2f_Relax')], 
                                  na.rm = TRUE)

# High Negative Vars:
dat.study2_model$EMA_2g_Nerv <- dat.study2_model$EMA_2g_Nerv + 1
dat.study2_model$EMA_2h_Worry <- dat.study2_model$EMA_2h_Worry + 1
dat.study2_model$EMA_2i_Afraid <- dat.study2_model$EMA_2i_Afraid + 1
dat.study2_model$ANX <- rowMeans(dat.study2_model[c('EMA_2g_Nerv', 
                                                    'EMA_2h_Worry', 
                                                    'EMA_2i_Afraid')], 
                                 na.rm = TRUE)

# Low Negative Vars: 
dat.study2_model$EMA_2m_Sad <- dat.study2_model$EMA_2m_Sad + 1
dat.study2_model$EMA_2o_Hopeless <- dat.study2_model$EMA_2o_Hopeless + 1
dat.study2_model$DEP <- rowMeans(dat.study2_model[c('EMA_2m_Sad', 
                                                    'EMA_2o_Hopeless')], 
                                 na.rm = TRUE)

# Angry/Frustrated Vars:
dat.study2_model$EMA_2j_Annoy <- dat.study2_model$EMA_2j_Annoy + 1
dat.study2_model$EMA_2k_Angry <- dat.study2_model$EMA_2k_Angry + 1
dat.study2_model$ANG <- rowMeans(dat.study2_model[c('EMA_2j_Annoy', 
                                                    'EMA_2k_Angry')], 
                                 na.rm = TRUE)

# Tired Vars:
dat.study2_model$EMA_2l_Slug <- dat.study2_model$EMA_2l_Slug + 1
dat.study2_model$EMA_2n_Tired <- dat.study2_model$EMA_2n_Tired + 1
dat.study2_model$TRD <- rowMeans(dat.study2_model[c('EMA_2l_Slug', 
                                                    'EMA_2n_Tired')], 
                                 na.rm = TRUE)

colnames(dat.study2_model)[22] <- "c.DN"  # Renaming for consistency with Study 1 scripts
dat.study2_model$c.DN <- as.numeric(scale(dat.study2_model$c.DN))

#---------------------------------------------------------------------------------------------------
# Getting aggregate scores for each ID
dat_EMA_means = dat.study2_model %>% 
  group_by(ID) %>% 
  summarize(JOY_m = mean(JOY, na.rm = TRUE), 
            JOY_sd = sd(JOY, na.rm = TRUE), 
            CALM_m = mean(CALM, na.rm = TRUE), 
            CALM_sd = sd(CALM, na.rm = TRUE),
            ANX_m = mean(ANX, na.rm = TRUE), 
            ANX_sd = sd(ANX, na.rm = TRUE),
            DEP_m = mean(DEP, na.rm = TRUE), 
            DEP_sd = sd(DEP, na.rm = TRUE),
            ANG_m = mean(ANG, na.rm = TRUE), 
            ANG_sd = sd(ANG, na.rm = TRUE),
            TRD_m = mean(TRD, na.rm = TRUE), 
            TRD_sd = sd(TRD, na.rm = TRUE))

# Exploratory data analysis
png(paste0(study2.graphics, '/S2_EDA_corrgram_summary_vars.png'), 
    res = 900, 
    units = "in", 
    height = 10, 
    width = 10)

corrgram(dat_EMA_means[,2:ncol(dat_EMA_means)], 
         order = TRUE, 
         lower.panel = panel.shade, 
         upper.panel = panel.conf, 
         col.regions = colorRampPalette(c("navy", "royalblue", "white", 
                                          "salmon", "red")))
title("Plot of Mood Variable Means & SDs for Study 2")
dev.off()

#---------------------------------------------------------------------------------------------------
dat.study2_model <- merge(dat.study2_model, dat_EMA_means, by = "ID")

psych::describe(dat.study2_model)

dat.study2_model$NegEvnt <- as.factor(dat.study2_model$NegEvnt)
dat.study2_model$PosEvnt <- as.factor(dat.study2_model$PosEvnt)

# The modeling syntax usese all summary data from the EMA mood items to generate posterior 
# distributions for individual missing values with random intercepts allowing for individual 
# differences in posterior distributions. Note that this imputation can take a while... 
fml <- NegEvnt + PosEvnt + JOY + CALM + ANX + DEP + TRD + ANG ~ c.DN + JOY_m + JOY_sd + CALM_m + CALM_sd +
  ANX_m + ANX_sd + DEP_m + DEP_sd + ANG_m + ANG_sd + TRD_m + TRD_sd + (1|ID) 

M <- 10 
imp <- jomoImpute(dat.study2_model, 
                  formula = fml, 
                  n.burn = 1000, 
                  n.iter = 1250, 
                  m = M, 
                  seed = 7022019)

dat.imp <- mitmlComplete(imp)
dat.study2_list <- list() 

for(i in 1:M){
  dat.study2_list[[i]] <- dat.imp[[i]]
  dat.study2_list[[i]]$NegEvnt_dich <- ifelse(dat.study2_list[[i]]$NegEvnt == 'No', 0, 1)
  dat.study2_list[[i]]$PosEvnt_dich <- ifelse(dat.study2_list[[i]]$PosEvnt == 'No', 0, 1)
  
  # Assigning only positive values to variables that will receive lognormal priors 
  # Included all negative momentary moood variables and JOY which was positively skewed
  # Adding some randomness in that if the value is below 0, it will be assigned a value from a 
  # uniform distribution betwee .01 and .99
  dat.study2_list[[i]]$ANX[dat.study2_list[[i]]$ANX <= 0] <- runif(sum(dat.study2_list[[i]]$ANX <= 0), .01, .99)
  dat.study2_list[[i]]$DEP[dat.study2_list[[i]]$DEP <= 0] <- runif(sum(dat.study2_list[[i]]$DEP <= 0), .01, .99)
  dat.study2_list[[i]]$ANG[dat.study2_list[[i]]$ANG <= 0] <- runif(sum(dat.study2_list[[i]]$ANG <= 0), .01, .99)  
  dat.study2_list[[i]]$TRD[dat.study2_list[[i]]$TRD <= 0] <- runif(sum(dat.study2_list[[i]]$TRD <= 0), .01, .99)
  dat.study2_list[[i]]$JOY[dat.study2_list[[i]]$JOY <= 0] <- runif(sum(dat.study2_list[[i]]$JOY <= 0), .01, .99)
  
  # Weighted averages from the imputed data values - there were 3 positive mood items per subscale
  # There were 3 anxiety and 2 depression items - necessitating the (slightly) more complicated 
  # weighting presented below. 
  dat.study2_list[[i]]$POS <- (dat.study2_list[[i]]$JOY + dat.study2_list[[i]]$CALM)/2
  dat.study2_list[[i]]$NEG <- (dat.study2_list[[i]]$ANX*3 + dat.study2_list[[i]]$DEP*2)/5 # need a weighted average here. 
  
  for(j in 1:length(IDs)){
    dat.study2_list[[i]]$prop.NegEvnt[dat.study2_list[[i]]$ID == IDs[j]] <- mean(dat.study2_list[[i]]$NegEvnt_dich[dat.study2_list[[i]]$ID == IDs[j]])
    dat.study2_list[[i]]$prop.PosEvnt[dat.study2_list[[i]]$ID == IDs[j]] <- mean(dat.study2_list[[i]]$PosEvnt_dich[dat.study2_list[[i]]$ID == IDs[j]])
  }
  dat.study2_list[[i]]$c.NegEvnt <- dat.study2_list[[i]]$NegEvnt_dich - dat.study2_list[[i]]$prop.NegEvnt
  dat.study2_list[[i]]$c.PosEvnt <- dat.study2_list[[i]]$PosEvnt_dich - dat.study2_list[[i]]$prop.PosEvnt
}

