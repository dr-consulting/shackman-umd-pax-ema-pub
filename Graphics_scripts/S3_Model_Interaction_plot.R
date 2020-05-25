library(bayesplot)
library(tidyverse)
library(tidybayes)
library(brms)

DATA_FILEPATH <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Data/study3_data.RData"
GRAPHICS_DIR <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/plot_outputs"
load(DATA_FILEPATH)
# IDs are made up - will assume new levels
# Holding constant all other properties between the two groups
# Will only var DN and the occurrence of a negative or positiv event...

new_data <- data.frame(ID=c("PAX999", "PAX999", "PAX999", "PAX999", "PAX998", "PAX998", "PAX998", "PAX998"), 
                       DN=c(-1, -1, -1, -1, 1, 1, 1, 1), 
                       certainty=rep(c("Certain", "Uncertain", "Certain", "Uncertain")),
                       valence=rep(c("Safety", "Safety", "Threat", "Threat")),
                       DN_cat=c("-1 SD", "-1 SD", "-1 SD", "-1 SD", "+1 SD", "+1 SD", "+1 SD", "+1 SD"), 
                       cond_cat=rep(c("Certain Safe", "Uncertain Safe", "Certain Threat", "Uncertain Threat")),
                       stringsAsFactors = FALSE)

#----------------------------------------------------------------------------------------------------------------------
# Negative Mood, Positive Event x DN Model
posterior_df <- posterior_predict(fit_log, newdata=new_data, allow_new_levels=TRUE)
colnames(posterior_df) <- paste(new_data$DN_cat, new_data$cond_cat, sep="_")

png(paste0(GRAPHICS_DIR, "/S3_Anx_Rating_DN_x_cond_ppd_interaction_plot_w_bars.png"), 
    units = "in", width = 6, height = 6, res=300)
posterior_df %>%
  as_tibble() %>% 
  select_all() %>% 
  pivot_longer(cols = all_of(colnames(posterior_df)), 
               values_to = "p_rating", 
               names_to = c("DN_cat", "cond_cat"), 
               names_sep = "_") %>% 
  mutate(cond_cat=fct_relevel(cond_cat, "Certain Safe", "Uncertain Safe", "Certain Threat", 
                              "Uncertain Threat")) %>% 
  ggplot(aes(x=DN_cat, y=p_rating, fill=cond_cat))+
  geom_violin(alpha=.15)+
  stat_summary(fun=mean, geom="bar", alpha=.85, color="black", position=position_dodge(.9))+
  stat_summary(fun.data=mean_hdci, position=position_dodge(.9), geom="pointrange", show.legend = FALSE)+
  stat_summary(fun.data=mean_hdci, position=position_dodge(.9), geom="linerange", fun.args=list(.width=.8), 
               size=1.5, show.legend = FALSE)+
  scale_fill_manual(values = c("Certain Safe"=RColorBrewer::brewer.pal(9, "Blues")[3],
                               "Uncertain Safe"=RColorBrewer::brewer.pal(9, "Blues")[7],
                               "Certain Threat"=RColorBrewer::brewer.pal(9, "Reds")[3],
                               "Uncertain Threat"=RColorBrewer::brewer.pal(9, "Reds")[7]), 
                    limits=c("Certain Safe", "Uncertain Safe", "Certain Threat", "Uncertain Threat"))+
  theme_bw()+
  labs(y="Posterior Predicted Anxiety Rating", 
       x="Dispositional Negativity", fill="Certainty/Threat Condition", 
       title = "Study 3 DN x Certainty/Threat Condition: Anxiety Rating")+
  scale_x_discrete(limits=c("-1 SD", "+1 SD"))+
  theme(legend.position = "bottom")+
  guides(fill=guide_legend(title.position = "top", title.hjust = .5))+
  coord_cartesian(ylim=c(0, 6))
dev.off()
