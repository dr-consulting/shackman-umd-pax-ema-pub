# Script used to create summary data for Study 2 to include in Supplement
# Involved loading a pre-saved but incomplete version of the EMA_Supplemental Data

load(paste0(data.folder, "/Emotion MS environment.RData"))

supp_descrip_study1 <- dat %>% 
  group_by(subid) %>% 
  summarize(DN = mean(ZAP_Both, na.rm=TRUE),
            NegEvnt_mean = mean(WorstEvent_Neg, na.rm = TRUE), 
            PosEvnt_mean = mean(BestEvent_Pos, na.rm = TRUE), 
            NEG_mean = mean(MAFS_NA, na.rm=TRUE),
            POS_mean = mean(MAFS_PA, na.rm=TRUE), 
            NegEvnt_sd = sd(WorstEvent_Neg, na.rm = TRUE), 
            PosEvnt_sd = sd(BestEvent_Pos, na.rm = TRUE), 
            NEG_sd = sd(MAFS_NA, na.rm=TRUE),
            POS_sd = sd(MAFS_PA, na.rm=TRUE))

load(paste0(data.folder, "/Study2_Clean.RData"))
supp_DN_study2 <- dat.lv2 %>% 
  rename(DN = DN_ScreenAndLab1) %>% 
  dplyr::select(ID, DN, Gender)

supp_Mood_items_study2 <- dat.study2 %>% 
  dplyr::select(ID, Enthus, Joy, Cheer, Calm, Content, Relax, Nerv, Worry, Afraid, Annoy, Angry, Slug,
                Tired, Sad, Hopeless)

supp_Mood_comps_study2 <- dat.study2 %>% 
  dplyr::select(ID, JOY, CALM, POS, NEG, ANX, DEP, ANG, TRD)

supp_descrip_study2 <- dat2_descrip

save(list = c("supp_DN_study1", "supp_Mood_study1", "supp_DN_study2", "supp_Mood_comps_study2", 
              "supp_Mood_items_study2", "dat.CFA", "dat.EFA", "supp_descrip_study2", 
              "supp_descrip_study1"), 
     file = "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Supplemental.RData")