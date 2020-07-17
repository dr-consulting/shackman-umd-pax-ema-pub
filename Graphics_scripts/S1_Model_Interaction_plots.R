source("~/dr-consulting_GH/shackman-umd-pax-ema-pub/Graphics_scripts/bayes_variance_riverplots_utils.R")
library(bayesplot)

POSTERIOR_PATH <- "/media/matthew/My Book/EMA_S1_Bayesian_Posteriors"
DATA_FILEPATH <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Data/study1_data.RData"
GRAPHICS_DIR <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/plot_outputs/"
data_loader(POSTERIOR_PATH, "S1_NEG_PosEvnt_Rct", data_filename = DATA_FILEPATH)
data_loader(POSTERIOR_PATH, "S1_POS_PosEvnt_Rct", data_filename = DATA_FILEPATH)

# IDs are made up - will assume new levels
# Holding constant all other properties between the two groups
# Will only var DN and the occurrence of a negative or positiv event...

sd_PosEvnt <- sd(dat.study1_model$c.PosEvnt, na.rm = TRUE)
m_PosEvnt <- mean(dat.study1_model$PosEvnt, na.rm = TRUE)
sd_NegEvnt <- sd(dat.study1_model$c.NegEvnt, na.rm = TRUE)
m_NegEvnt <- mean(dat.study1_model$NegEvnt, na.rm = TRUE)

new_data <- data.frame(ID=c("PAX999", "PAX999", "PAX998", "PAX998"), 
                       c.DN=c(-1, -1, 1, 1), 
                       c.PosEvnt=c(-sd_PosEvnt, sd_PosEvnt, -sd_PosEvnt, sd_PosEvnt), 
                       c.NegEvnt=c(-sd_NegEvnt, sd_NegEvnt, -sd_NegEvnt, sd_NegEvnt),
                       DN_cat=c("-1 SD", "-1 SD", "+1 SD", "+1 SD"), 
                       PosEvnt_cat=c("Low Intensity\nPositive Event", "High Intensity\nPositive Event", 
                                      "Low Intensity\nPositive Event", "High Intensity\nPositive Event"),
                       NegEvnt_cat=c("Low Intensity\nNegative Event", "High Intensity\nNegative Event", 
                                      "Low Intensity\nNegative Event", "High Intensity\nNegative Event"),
                       m.POS = mean(dat.study1_model$m.POS), 
                       m.NEG = mean(dat.study1_model$m.NEG), 
                       sd.POS = mean(dat.study1_model$sd.POS), 
                       sd.NEG = mean(dat.study1_model$sd.NEG), 
                       sd.PosEvnt = sd(dat.study1_model$c.PosEvnt, na.rm = TRUE),
                       m.PosEvnt = mean(dat.study1_model$PosEvnt, na.rm = TRUE),
                       sd.NegEvnt = sd(dat.study1_model$c.NegEvnt, na.rm = TRUE),
                       m.NegEvnt = mean(dat.study1_model$NegEvnt, na.rm = TRUE),
                       stringsAsFactors = FALSE)

#----------------------------------------------------------------------------------------------------------------------
# Negative Mood, Positive Event x DN Model
posterior_df <- posterior_predict(S1_NEG_PosEvnt_Rct, newdata=new_data, allow_new_levels=TRUE)[,,"NEG"]
colnames(posterior_df) <- paste(new_data$DN_cat, new_data$PosEvnt_cat, sep="_")

S1_NEG_DN_x_PosEvnt_plot <- posterior_df %>%
    as_tibble() %>% 
    select_all() %>% 
    pivot_longer(cols = all_of(colnames(posterior_df)), 
                 values_to = "p_NEG", 
                 names_to = c("DN_cat", "PosEvnt_cat"), 
                 names_sep = "_") %>% 
    mutate(
      PosEvnt_cat = forcats::fct_relevel(PosEvnt_cat, c("Low Intensity\nPositive Event", 
                                                        "High Intensity\nPositive Event"))
    ) %>% 
    ggplot(aes(x=DN_cat, y=p_NEG, fill=PosEvnt_cat))+
    geom_violin(alpha=.15)+
    stat_summary(fun=mean, geom="bar", alpha=.85, color="black", position=position_dodge(.9))+
    stat_summary(fun.data=mean_hdci, position=position_dodge(.9), geom="pointrange", show.legend = FALSE)+
    stat_summary(fun.data=mean_hdci, position=position_dodge(.9), geom="linerange", fun.args=list(.width=.8), 
                 size=1.5, show.legend = FALSE)+
    scale_fill_manual(values = c("Low Intensity\nPositive Event"=RColorBrewer::brewer.pal(9, "Blues")[5], 
                                 "High Intensity\nPositive Event"=RColorBrewer::brewer.pal(9, "Reds")[5]))+
    theme_bw()+
    labs(y="Posterior Predicted Negative Mood Composite", 
         x="Dispositional Negativity Score", fill="Positive Event Rating", 
         title = "Study 1 DN x Positive Event Interaction: Negative Mood")+
    scale_x_discrete(limits=c("-1 SD", "+1 SD"))+
    theme(legend.position = "bottom")+
    guides(fill=guide_legend(title.position = "top", title.hjust = .5))+
    coord_cartesian(ylim=c(0, 6))

# Save png version of plot
png(paste(GRAPHICS_DIR, "S1_NEG_DN_x_PosEvnt_ppd_interaction_plot_w_bars.png"), units = "in", width = 6, height = 6, 
    res=300)  
S1_NEG_DN_x_PosEvnt_plot
dev.off()

# Save postscript (.eps) version of plot
postscript(paste(GRAPHICS_DIR, "S1_NEG_DN_x_PosEvnt_ppd_interaction_plot_w_bars.eps"), width = 6, height = 6)  
S1_NEG_DN_x_PosEvnt_plot
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Positive Mood, Positive Event x DN Model
posterior_df <- posterior_predict(S1_POS_PosEvnt_Rct, newdata=new_data, allow_new_levels=TRUE)[,,"POS"]
colnames(posterior_df) <- paste(new_data$DN_cat, new_data$PosEvnt_cat, sep="_")

S1_POS_DN_x_PosEvnt_plot <- posterior_df %>%
  as_tibble() %>% 
  select_all() %>% 
  pivot_longer(cols = all_of(colnames(posterior_df)), 
               values_to = "p_POS", 
               names_to = c("DN_cat", "PosEvnt_cat"), 
               names_sep = "_") %>% 
  mutate(
    PosEvnt_cat = forcats::fct_relevel(PosEvnt_cat, c("Low Intensity\nPositive Event", 
                                                      "High Intensity\nPositive Event"))
  ) %>% 
  ggplot(aes(x=DN_cat, y=p_POS, fill=PosEvnt_cat))+
  geom_violin(alpha=.15)+
  stat_summary(fun=mean, geom="bar", alpha=.85, color="black", position=position_dodge(.9))+
  stat_summary(fun.data=mean_hdci, position=position_dodge(.9), geom="pointrange", show.legend = FALSE)+
  stat_summary(fun.data=mean_hdci, position=position_dodge(.9), geom="linerange", fun.args=list(.width=.8), 
               size=1.5, show.legend = FALSE)+
  scale_fill_manual(values = c("Low Intensity\nPositive Event"=RColorBrewer::brewer.pal(9, "Blues")[5], 
                               "High Intensity\nPositive Event"=RColorBrewer::brewer.pal(9, "Reds")[5]))+
  theme_bw()+
  labs(y="Posterior Predicted Positive Mood Composite", 
       x="Dispositional Negativity Score", fill="Positive Event Rating", 
       title = "Study 1 DN x Positive Event Interaction: Positive Mood")+
  scale_x_discrete(limits=c("-1 SD", "+1 SD"))+
  theme(legend.position = "bottom")+
  guides(fill=guide_legend(title.position = "top", title.hjust = .5))+
  coord_cartesian(ylim=c(0, 6))

# Save png version of plot
png(paste(GRAPHICS_DIR, "S1_POS_DN_x_PosEvnt_ppd_interaction_plot_w_bars.png"), units = "in", width = 6, height = 6, 
    res=300)  
S1_POS_DN_x_PosEvnt_plot
dev.off()

# Save postscript (.eps) version of plot
postscript(paste(GRAPHICS_DIR, "S1_POS_DN_x_PosEvnt_ppd_interaction_plot_w_bars.eps"), width = 6, height = 6)  
S1_POS_DN_x_PosEvnt_plot
dev.off()
