library(bayesplot)
library(tidyverse)
library(tidybayes)
library(brms)

GRAPHICS_DIR <- "~/github/ATNL/shackman-umd-pax-ema-pub/plot_outputs"
# IDs are made up - will assume new levels
# Holding constant all other properties between the two groups
# Will only var DN and the occurrence of a negative or positive event...

new_data <- data.frame(ID=c("PAX999", "PAX999", "PAX999", "PAX999", "PAX998", "PAX998", "PAX998", "PAX998"), 
                       DN=c(-1, -1, -1, -1, 1, 1, 1, 1), 
                       certainty=rep(c("Certain", "Uncertain", "Certain", "Uncertain")),
                       valence=rep(c("Safety", "Safety", "Threat", "Threat")),
                       DN_cat=c("-1 SD", "-1 SD", "-1 SD", "-1 SD", "+1 SD", "+1 SD", "+1 SD", "+1 SD"), 
                       cond_cat=rep(c("Certain Safe", "Uncertain Safe", "Certain Threat", "Uncertain Threat")),
                       stringsAsFactors = FALSE)

#----------------------------------------------------------------------------------------------------------------------
# ANOVA model for S3 
# posterior_df <- dat.study3[[c('ID', ]]
posterior_df <- posterior_epred(fit_gamma, newdata = new_data, allow_new_levels=TRUE, re_formula = NA) %>% 
  as.data.frame()

colnames(posterior_df) <- paste(new_data$DN_cat, new_data$cond_cat, sep="_")

S3_ANOVA_plot <- posterior_df %>%
  as_tibble() %>% 
  select_all() %>% 
  pivot_longer(cols = all_of(colnames(posterior_df)), 
               values_to = "p_rating", 
               names_to = c("DN_cat", "cond_cat"), 
               names_sep = "_") %>% 
  mutate(cond_cat=fct_relevel(cond_cat, "Certain Safe", "Uncertain Safe", "Certain Threat", 
                              "Uncertain Threat")) %>% 
  ggplot(aes(x=DN_cat, y=p_rating, fill=cond_cat))+
  stat_summary(fun=mean, geom="bar", color="black", position=position_dodge(.9), show.legend = FALSE) +
  stat_summary(fun.data=mean_hdi, position=position_dodge(.9), geom="pointrange", show.legend = FALSE, size=.5) +
  stat_summary(fun.data=mean_hdi, position=position_dodge(.9), geom="linerange", fun.args=list(.width=.8), 
               size=1.5, show.legend = FALSE)+
  scale_fill_manual(values = c("Certain Safe"=RColorBrewer::brewer.pal(9, "Blues")[5],
                               "Uncertain Safe"=RColorBrewer::brewer.pal(9, "Purples")[5], 
                               "Certain Threat"=RColorBrewer::brewer.pal(9, "Oranges")[5],
                               "Uncertain Threat"=RColorBrewer::brewer.pal(9, "Reds")[5]))+
  theme_bw()+
  labs(y="Posterior Predicted Anxiety Rating", 
       x="Dispositional Negativity Score", fill="Certainty/Threat", 
       title = "Study 3 DN x Certainty/Threat Condition: Anxiety Rating")+
  scale_x_discrete(limits=c("-1 SD", "+1 SD"))+
  scale_y_continuous(limits=c(0, 4)) +
  geom_hline(yintercept = 1, lty='dashed', color='black') +
  coord_flip() +
  annotate('text', x=1.60, y=0, label='Certain Safe', hjust=-.05) +
  annotate('text', x=1.825, y=0, label='Uncertain Safe', hjust=-.05) +
  annotate('text', x=2.05, y=0, label='Certain Threat', hjust=-.05) +
  annotate('text', x=2.275, y=0, label='Uncertain Threat', hjust=-.05) +
  annotate('label', x=2.225, y=1.5, label='Minimum scale value', angle=90, fill='grey')

S3_ANOVA_plot

bayesplot::mcmc_areas(posterior_df)

# png version of the plot... 
png(paste0(GRAPHICS_DIR, "/S3_Anx_Rating_DN_x_cond_ppd_interaction_plot_w_bars.png"), 
    units = "in", width = 8, height = 6, res=300)
S3_ANOVA_plot
dev.off()

# eps version of the plot - no transparency
postscript(paste0(GRAPHICS_DIR, "/S3_Anx_Rating_DN_x_cond_ppd_interaction_plot_w_bars.eps"), width = 8, height = 6)
S3_ANOVA_plot
dev.off()
