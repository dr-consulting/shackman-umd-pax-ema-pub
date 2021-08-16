source("~/github/ATNL/shackman-umd-pax-ema-pub/Graphics_scripts/bayes_variance_riverplots_utils.R")
library(bayesplot)

POSTERIOR_PATH <- "/media/dr-owner/HDD1/EMA_S1_Bayesian_Posteriors"
DATA_FILEPATH <- "~/github/ATNL/shackman-umd-pax-ema-pub/Data/study1_data.RData"
GRAPHICS_DIR <- "~/github/ATNL/shackman-umd-pax-ema-pub/plot_outputs"
data_loader(POSTERIOR_PATH, "/gamma/S1_NEG_PosEvnt_Rct", data_filename = DATA_FILEPATH)
data_loader(POSTERIOR_PATH, "/gaussian/S1_POS_PosEvnt_Rct", data_filename = DATA_FILEPATH)

# Bringing in some utilities to get imputed mean values of certain effects
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

prop_NegEvnt <- get_imputed_mean(dat.study1_list, "prop_NegEvnt", FALSE)
NegEvnt_yes <- 1 - prop_NegEvnt
NegEvnt_no <- 0 - prop_NegEvnt

prop_PosEvnt <- get_imputed_mean(dat.study1_list, "prop_PosEvnt", FALSE)
PosEvnt_yes <- 1 - prop_PosEvnt
PosEvnt_no <- 0 - prop_PosEvnt

new_data <- data.frame(ID=c("PAX999", "PAX999", "PAX998", "PAX998"), 
                       c.DN=c(-1, -1, 1, 1), 
                       c.PosEvnt=c(PosEvnt_no, PosEvnt_yes, PosEvnt_no, PosEvnt_yes), 
                       c.NegEvnt=c(NegEvnt_no, NegEvnt_yes, NegEvnt_no, NegEvnt_yes),
                       DN_cat=c("-1 SD", "-1 SD", "+1 SD", "+1 SD"), 
                       PosEvnt_cat=c("Low Intensity\nPositive Event", "High Intensity\nPositive Event", 
                                      "Low Intensity\nPositive Event", "High Intensity\nPositive Event"),
                       NegEvnt_cat=c("Low Intensity\nNegative Event", "High Intensity\nNegative Event", 
                                      "Low Intensity\nNegative Event", "High Intensity\nNegative Event"),
                       prop_PosEvnt = prop_PosEvnt, 
                       prop_NegEvnt = prop_NegEvnt,
                       stringsAsFactors = FALSE)

#----------------------------------------------------------------------------------------------------------------------
# Negative Mood, Positive Event x DN Model
posterior_df <- posterior_epred(S1_NEG_PosEvnt_Rct, newdata = new_data, allow_new_levels=TRUE, re_formula = NA) %>% 
  as.data.frame()
posterior_df <- posterior_df - 1 
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
    stat_summary(fun=mean, geom="bar", color="black", position=position_dodge(.9), show.legend = FALSE)+
    stat_summary(fun.data=mean_hdci, position=position_dodge(.9), geom="pointrange", show.legend = FALSE, size=.75)+
    stat_summary(fun.data=mean_hdci, position=position_dodge(.9), geom="linerange", fun.args=list(.width=.8), 
                 size=1.5, show.legend = FALSE)+
    scale_fill_manual(values = c("Low Intensity\nPositive Event"="#426ebd", 
                                 "High Intensity\nPositive Event"="#ad580b"))+
    theme_bw()+
    labs(y="Posterior Predicted Negative Mood Composite", 
         x="Dispositional Negativity Score", fill="Positive Event Rating", 
         title = "Study 1 DN x Positive Event Interaction: Negative Mood")+
    scale_x_discrete(limits=c("-1 SD", "+1 SD"))+
    scale_y_continuous(limits=c(0, 4)) +
    coord_flip() +
    annotate('text', x=2.25, y=0, label='High Intensity\nPositive Event', hjust=-.05) +
    annotate('text', x=1.80, y=0, label='Low Intensity\nPositive Event', hjust=-.05) + 
    annotate('label', x=2, y=2, label='Larger difference in negative mood', fill='grey') +
    geom_segment(x = .5, y = 2, xend = .5, yend = 3.95,
                 arrow = arrow(length = unit(0.25, "cm"), type='closed'), 
                 color = 'grey') +
    geom_hline(yintercept = 4, lty='dashed', color='grey') +
    annotate('text', x=.5, y=1.75, label='Max rating', color='grey')
  

S1_NEG_DN_x_PosEvnt_plot

# Save png version of plot
png(paste0(GRAPHICS_DIR, "/S1_NEG_DN_x_PosEvnt_ppd_interaction_plot_w_bars.png"), units = "in", width = 8, height = 6, 
    res=300)  
S1_NEG_DN_x_PosEvnt_plot
dev.off()

# Save postscript (.eps) version of plot
postscript(paste0(GRAPHICS_DIR, "/S1_NEG_DN_x_PosEvnt_ppd_interaction_plot_w_bars.eps"), width = 8, height = 6)  
S1_NEG_DN_x_PosEvnt_plot
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# Positive Mood, Positive Event x DN Model
posterior_df <- posterior_epred(S1_POS_PosEvnt_Rct, newdata = new_data, allow_new_levels=TRUE, re_formula = NA) %>% 
  as.data.frame()
posterior_df <- posterior_df - 1
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
  stat_summary(fun=mean, geom="bar", color="black", position=position_dodge(.9), show.legend = FALSE)+
  stat_summary(fun.data=mean_hdci, position=position_dodge(.9), geom="pointrange", show.legend = FALSE, size=.75)+
  stat_summary(fun.data=mean_hdci, position=position_dodge(.9), geom="linerange", fun.args=list(.width=.8), 
               size=1.5, show.legend = FALSE)+
  scale_fill_manual(values = c("Low Intensity\nPositive Event"="#426ebd", 
                               "High Intensity\nPositive Event"="#ad580b"))+
  theme_bw()+
  labs(y="Posterior Predicted Positive Mood Composite", 
       x="Dispositional Negativity Score", fill="Positive Event Rating", 
       title = "Study 1 DN x Positive Event Interaction: Positive Mood")+
  scale_x_discrete(limits=c("-1 SD", "+1 SD"))+
  scale_y_continuous(limits=c(0, 4)) +
  coord_flip() +
  annotate('text', x=2.25, y=0, label='High Intensity\nPositive Event', hjust=-.05) +
  annotate('text', x=1.80, y=0, label='Low Intensity\nPositive Event', hjust=-.05) + 
  annotate('label', x=2, y=3, label='Larger difference in positive mood', fill='grey') +
  geom_segment(x = .5, y = 3, xend = .5, yend = 3.95,
               arrow = arrow(length = unit(0.25, "cm"), type='closed'), 
               color = 'grey') +
  geom_hline(yintercept = 4, lty='dashed', color='grey') +
  annotate('text', x=.5, y=2.75, label='Max rating', color='grey')

S1_POS_DN_x_PosEvnt_plot

# Save png version of plot
png(paste0(GRAPHICS_DIR, "/S1_POS_DN_x_PosEvnt_ppd_interaction_plot_w_bars.png"), units = "in", width = 8, height = 6, 
    res=300)  
S1_POS_DN_x_PosEvnt_plot
dev.off()

# Save postscript (.eps) version of plot
postscript(paste0(GRAPHICS_DIR, "/S1_POS_DN_x_PosEvnt_ppd_interaction_plot_w_bars.eps"), width = 8, height = 6)  
S1_POS_DN_x_PosEvnt_plot
dev.off()
