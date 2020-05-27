# Study   NA/PA   NE/PE   Tonic%  Exposure%   Reactivity%   Model_DN  Model_Agg_Exp   Model_Moment_Exp  Model_Reactivity
library(tidyverse)
library(cowplot)
# Created this matrix manually to try to summarize all of model effects

summary_matrix <-c("Study 1", "NA", "NE", 1.24, 1.63, 0.18, "pos", "pos", "pos", "nsg",
                   "Study 2", "NA", "NE", 0.73, 0.03, 0.05, "pos", "pos", "pos", "pos",
                   "Study 1", "NA", "PE", 2.72, 0.06, 0.33, "pos", "nsg", "neg", "neg",
                   "Study 2", "NA", "PE", 0.76, 0.00, 0.21, "pos", "nsg", "neg", "neg",
                   "Study 1", "PA", "NE", 5.34, 1.23, 0.00, "neg", "nsg", "neg", "nsg",
                   "Study 2", "PA", "NE", 6.35, 0.00, 0.00, "neg", "nsg", "neg", "nsg",
                   "Study 1", "PA", "PE", 6.25, 0.24, 0.38, "neg", "pos", "pos", "pos", 
                   "Study 2", "PA", "PE", 4.97, 0.83, 0.30, "neg", "pos", "pos", "nsg")

summary_matrix <- matrix(summary_matrix, byrow = TRUE, nrow = 8)

colnames(summary_matrix) <- c("Study", "Mood", "Event", "Tonic_var", "Exposure_var", "Reactivity_var", 
                              "DN_credible", "Agg_exp_credible", "Moment_exp_credible", "Reactivity_credible")

summary_df <- data.frame(summary_matrix, stringsAsFactors = FALSE)

# Horizontal bar plot of Study 1 and Study 2 Effects for Negative Mood/Negative Events Decomposition
plot_NA_NE <- 
summary_df %>% 
  mutate(
    Tonic_var = as.numeric(Tonic_var), 
    Exposure_var = as.numeric(Exposure_var), 
    Reactivity_var = as.numeric(Reactivity_var)
  ) %>% 
  mutate(
    Tonic_var = Tonic_var/rowSums(select(.,Tonic_var, Exposure_var, Reactivity_var)), 
    Exposure_var = Exposure_var/rowSums(select(.,Tonic_var, Exposure_var, Reactivity_var)),
    Reactivity_var = Reactivity_var/rowSums(select(.,Tonic_var, Exposure_var, Reactivity_var))
  ) %>% 
  select(Study, Mood, Event, Tonic_var, Exposure_var, Reactivity_var) %>% 
  rename(
    Tonic = Tonic_var, 
    Exposure = Exposure_var, 
    Reactivity = Reactivity_var
  ) %>% 
  pivot_longer(cols = Tonic:Reactivity, names_to = "var_source", values_to = "prop_var") %>% 
  mutate(
    var_source = fct_relevel(var_source, "Tonic", "Exposure", "Reactivity")
  ) %>% 
  filter(Mood == "NA" & Event == "NE") %>% 
  ggplot(aes(x = Study, y = prop_var))+
  geom_bar(aes(fill = var_source), stat = "identity")+
  scale_x_discrete(limits = c("Study 2", "Study 1"))+
  scale_fill_manual(name = NULL, 
                    values=c("Tonic"=RColorBrewer::brewer.pal(9, "Blues")[7], 
                             "Exposure"=RColorBrewer::brewer.pal(9, "Purples")[7], 
                             "Reactivity"=RColorBrewer::brewer.pal(9, "Reds")[7]), 
                    breaks = c("Tonic", "Exposure", "Reactivity"))+
  labs(y=NULL, x=NULL)+
  theme_void()+
  theme(plot.title = element_text(hjust = .5), axis.text.y = element_text())+
  coord_flip()

# Horizontal bar plot of Study 1 and Study 2 Effects for Negative Mood/Positive Events Decomposition
plot_NA_PE <- 
  summary_df %>% 
  mutate(
    Tonic_var = as.numeric(Tonic_var), 
    Exposure_var = as.numeric(Exposure_var), 
    Reactivity_var = as.numeric(Reactivity_var)
  ) %>% 
  mutate(
    Tonic_var = Tonic_var/rowSums(select(.,Tonic_var, Exposure_var, Reactivity_var)), 
    Exposure_var = Exposure_var/rowSums(select(.,Tonic_var, Exposure_var, Reactivity_var)),
    Reactivity_var = Reactivity_var/rowSums(select(.,Tonic_var, Exposure_var, Reactivity_var))
  ) %>% 
  select(Study, Mood, Event, Tonic_var, Exposure_var, Reactivity_var) %>% 
  rename(
    Tonic = Tonic_var, 
    Exposure = Exposure_var, 
    Reactivity = Reactivity_var
  ) %>% 
  pivot_longer(cols = Tonic:Reactivity, names_to = "var_source", values_to = "prop_var") %>% 
  mutate(
    var_source = fct_relevel(var_source, "Tonic", "Exposure", "Reactivity")
  ) %>% 
  filter(Mood == "NA" & Event == "PE") %>% 
  ggplot(aes(x = Study, y = prop_var))+
  geom_bar(aes(fill = var_source), stat = "identity")+
  scale_x_discrete(limits = c("Study 2", "Study 1"))+
  scale_fill_manual(name = NULL, 
                    values=c("Tonic"=RColorBrewer::brewer.pal(9, "Blues")[7], 
                             "Exposure"=RColorBrewer::brewer.pal(9, "Purples")[7], 
                             "Reactivity"=RColorBrewer::brewer.pal(9, "Reds")[7]), 
                    breaks = c("Tonic", "Exposure", "Reactivity"))+
  labs(y=NULL, x=NULL)+
  theme_void()+
  theme(plot.title = element_text(hjust = .5), axis.text.y = element_text())+
  coord_flip()

# Horizontal bar plot of Study 1 and Study 2 Effects for Negative Mood/Negative Events Decomposition
plot_PA_NE <- 
  summary_df %>% 
  mutate(
    Tonic_var = as.numeric(Tonic_var), 
    Exposure_var = as.numeric(Exposure_var), 
    Reactivity_var = as.numeric(Reactivity_var)
  ) %>% 
  mutate(
    Tonic_var = Tonic_var/rowSums(select(.,Tonic_var, Exposure_var, Reactivity_var)), 
    Exposure_var = Exposure_var/rowSums(select(.,Tonic_var, Exposure_var, Reactivity_var)),
    Reactivity_var = Reactivity_var/rowSums(select(.,Tonic_var, Exposure_var, Reactivity_var))
  ) %>% 
  select(Study, Mood, Event, Tonic_var, Exposure_var, Reactivity_var) %>% 
  rename(
    Tonic = Tonic_var, 
    Exposure = Exposure_var, 
    Reactivity = Reactivity_var
  ) %>% 
  pivot_longer(cols = Tonic:Reactivity, names_to = "var_source", values_to = "prop_var") %>% 
  mutate(
    var_source = fct_relevel(var_source, "Tonic", "Exposure", "Reactivity")
  ) %>% 
  filter(Mood == "PA" & Event == "NE") %>% 
  ggplot(aes(x = Study, y = prop_var))+
  geom_bar(aes(fill = var_source), stat = "identity")+
  scale_x_discrete(limits = c("Study 2", "Study 1"))+
  scale_fill_manual(name = NULL, 
                    values=c("Tonic"=RColorBrewer::brewer.pal(9, "Blues")[7], 
                             "Exposure"=RColorBrewer::brewer.pal(9, "Purples")[7], 
                             "Reactivity"=RColorBrewer::brewer.pal(9, "Reds")[7]), 
                    breaks = c("Tonic", "Exposure", "Reactivity"))+
  labs(y=NULL, x=NULL)+
  theme_void()+
  theme(plot.title = element_text(hjust = .5), axis.text.y = element_text())+
  coord_flip()

# Horizontal bar plot of Study 1 and Study 2 Effects for Negative Mood/Positive Events Decomposition
plot_PA_PE <- 
  summary_df %>% 
  mutate(
    Tonic_var = as.numeric(Tonic_var), 
    Exposure_var = as.numeric(Exposure_var), 
    Reactivity_var = as.numeric(Reactivity_var)
  ) %>% 
  mutate(
    Tonic_var = Tonic_var/rowSums(select(.,Tonic_var, Exposure_var, Reactivity_var)), 
    Exposure_var = Exposure_var/rowSums(select(.,Tonic_var, Exposure_var, Reactivity_var)),
    Reactivity_var = Reactivity_var/rowSums(select(.,Tonic_var, Exposure_var, Reactivity_var))
  ) %>% 
  select(Study, Mood, Event, Tonic_var, Exposure_var, Reactivity_var) %>% 
  rename(
    Tonic = Tonic_var, 
    Exposure = Exposure_var, 
    Reactivity = Reactivity_var
  ) %>% 
  pivot_longer(cols = Tonic:Reactivity, names_to = "var_source", values_to = "prop_var") %>% 
  mutate(
    var_source = fct_relevel(var_source, "Tonic", "Exposure", "Reactivity")
  ) %>% 
  filter(Mood == "PA" & Event == "PE") %>% 
  ggplot(aes(x = Study, y = prop_var))+
  geom_bar(aes(fill = var_source), stat = "identity")+
  scale_x_discrete(limits = c("Study 2", "Study 1"))+
  scale_fill_manual(name = NULL, 
                    values=c("Tonic"=RColorBrewer::brewer.pal(9, "Blues")[7], 
                             "Exposure"=RColorBrewer::brewer.pal(9, "Purples")[7], 
                             "Reactivity"=RColorBrewer::brewer.pal(9, "Reds")[7]), 
                    breaks = c("Tonic", "Exposure", "Reactivity"))+
  labs(y=NULL, x=NULL)+
  theme_void()+
  theme(plot.title = element_text(hjust = .5), axis.text.y = element_text())+
  coord_flip()


var_legend <- cowplot::get_legend(
  plot_NA_NE + 
    guides(fill = guide_legend(title = "Variance Pathway", 
                               title.position = "bottom", title.hjust = .5, nrow=1, reverse = TRUE)) +
    theme(legend.position = "bottom")
)

cowplot::plot_grid(plot_NA_NE + theme(legend.position = "none"), 
                   plot_NA_PE + theme(legend.position = "none"), 
                   plot_PA_NE + theme(legend.position = "none"), 
                   plot_PA_PE + theme(legend.position = "none"),
                   var_legend,
                   rel_heights = c(3, 3, 3, 3, 1),
                   ncol = 1)


#######################################################################################################################
# Now attempting to plot significance profiles for model effects
eff_plot_NA_NE <- 
summary_df %>% 
  select(Study, Mood, Event, DN_credible, Agg_exp_credible, Moment_exp_credible, Reactivity_credible) %>% 
  pivot_longer(cols = DN_credible:Reactivity_credible, names_to = "model_eff", values_to = "credible") %>% 
  mutate(
    model_eff = fct_relevel(model_eff, "DN_credible", "Agg_exp_credible", "Moment_exp_credible", 
                            "Reactivity_credible")
  ) %>% 
  filter(Mood == "NA" & Event == "NE") %>% 
  ggplot(aes(x=model_eff, y=Study, fill=credible)) + 
  geom_tile(width=.90, height=.90)+
  theme_void()+
  theme(axis.text = element_text(), 
        plot.title = element_text(hjust = .5))+
  scale_y_discrete(limits=c("Study 2", "Study 1"))+
  scale_x_discrete(labels=c("DN_credible"="DN", "Agg_exp_credible"="Aggregated\nExposure", 
                            "Moment_exp_credible"="Momentary\nExposure", "Reactivity_credible"="Reactivity"))+
  scale_fill_manual(values=c("nsg"="grey50", "pos"=RColorBrewer::brewer.pal(9, "Reds")[7], 
                             "neg"=RColorBrewer::brewer.pal(9, "Blues")[7]))

eff_plot_NA_PE <- 
summary_df %>% 
  select(Study, Mood, Event, DN_credible, Agg_exp_credible, Moment_exp_credible, Reactivity_credible) %>% 
  pivot_longer(cols = DN_credible:Reactivity_credible, names_to = "model_eff", values_to = "credible") %>% 
  mutate(
    model_eff = fct_relevel(model_eff, "DN_credible", "Agg_exp_credible", "Moment_exp_credible", 
                            "Reactivity_credible")
  ) %>% 
  filter(Mood == "NA" & Event == "PE") %>% 
  ggplot(aes(x=model_eff, y=Study, fill=credible)) + 
  geom_tile(width=.90, height=.90)+
  theme_void()+
  theme(axis.text = element_text(), 
        plot.title = element_text(hjust = .5))+
  scale_y_discrete(limits=c("Study 2", "Study 1"))+
  scale_x_discrete(labels=c("DN_credible"="DN", "Agg_exp_credible"="Aggregated\nExposure", 
                            "Moment_exp_credible"="Momentary\nExposure", "Reactivity_credible"="Reactivity"))+
  scale_fill_manual(values=c("nsg"="grey50", "pos"=RColorBrewer::brewer.pal(9, "Reds")[7], 
                             "neg"=RColorBrewer::brewer.pal(9, "Blues")[7]))

eff_plot_PA_NE <- 
summary_df %>% 
  select(Study, Mood, Event, DN_credible, Agg_exp_credible, Moment_exp_credible, Reactivity_credible) %>% 
  pivot_longer(cols = DN_credible:Reactivity_credible, names_to = "model_eff", values_to = "credible") %>% 
  mutate(
    model_eff = fct_relevel(model_eff, "DN_credible", "Agg_exp_credible", "Moment_exp_credible", 
                            "Reactivity_credible")
  ) %>% 
  filter(Mood == "PA" & Event == "NE") %>% 
  ggplot(aes(x=model_eff, y=Study, fill=credible)) + 
  geom_tile(width=.90, height=.90)+
  theme_void()+
  theme(axis.text = element_text(), 
        plot.title = element_text(hjust = .5))+
  scale_y_discrete(limits=c("Study 2", "Study 1"))+
  scale_x_discrete(labels=c("DN_credible"="DN", "Agg_exp_credible"="Aggregated\nExposure", 
                            "Moment_exp_credible"="Momentary\nExposure", "Reactivity_credible"="Reactivity"))+
  scale_fill_manual(values=c("nsg"="grey50", "pos"=RColorBrewer::brewer.pal(9, "Reds")[7], 
                             "neg"=RColorBrewer::brewer.pal(9, "Blues")[7]))

eff_plot_PA_PE <- 
summary_df %>% 
  select(Study, Mood, Event, DN_credible, Agg_exp_credible, Moment_exp_credible, Reactivity_credible) %>% 
  pivot_longer(cols = DN_credible:Reactivity_credible, names_to = "model_eff", values_to = "credible") %>% 
  mutate(
    model_eff = fct_relevel(model_eff, "DN_credible", "Agg_exp_credible", "Moment_exp_credible", 
                            "Reactivity_credible")
  ) %>% 
  filter(Mood == "PA" & Event == "PE") %>% 
  ggplot(aes(x=model_eff, y=Study, fill=credible)) + 
  geom_tile(width=.90, height=.90)+
  theme_void()+
  theme(axis.text = element_text(), 
        plot.title = element_text(hjust = .5))+
  scale_y_discrete(limits=c("Study 2", "Study 1"))+
  scale_x_discrete(labels=c("DN_credible"="DN", "Agg_exp_credible"="Aggregated\nExposure", 
                            "Moment_exp_credible"="Momentary\nExposure", "Reactivity_credible"="Reactivity"))+
  scale_fill_manual(values=c("nsg"="grey50", "pos"=RColorBrewer::brewer.pal(9, "Reds")[7], 
                             "neg"=RColorBrewer::brewer.pal(9, "Blues")[7]), 
                    labels=c("neg"="Negative", "nsg"="Not Credible", "pos"="Positive"))

eff_legend <- cowplot::get_legend(
  eff_plot_PA_PE + 
    guides(fill = guide_legend(title = "Model Effect Direction", 
                               title.position = "bottom", title.hjust = .5, nrow=1)) +
    theme(legend.position = "bottom")
)

var_col_title <- ggdraw()+
  draw_label(
    "Variance Decomposition", 
    hjust = .5
  ) +   
  theme(
        plot.margin = margin(0, 0, 0, 7)
  )


eff_col_title <- ggdraw()+
  draw_label(
    "Model Effects: Directionality/Credibility", 
    hjust = .5
  )+   
  theme(
    plot.margin = margin(0, 0, 0, 7)
  )


main_title_row <- plot_grid(var_col_title, eff_col_title, 
                            align = "vh", 
                            axis="tblr", 
                            ncol=2)

row_1_title <- ggdraw() +
  draw_label(
    "Negative Mood and Negative Events Models", 
    hjust = .5
  )
  
plot_row_1 <- cowplot::plot_grid(plot_NA_NE + theme(legend.position = "none"), 
                                 eff_plot_NA_NE + theme(legend.position = "none", 
                                                        axis.text.x = element_blank(), 
                                                        axis.text.y = element_blank()),
                                 align = "vh", 
                                 axis="tblr", 
                                 ncol=2)
row_2_title <- ggdraw()+
  draw_label(
    "Negative Mood and Positive Events Models", 
    hjust = .5
  )

plot_row_2 <- cowplot::plot_grid(plot_NA_PE + theme(legend.position = "none"), 
                                 eff_plot_NA_PE + theme(legend.position = "none", 
                                                        axis.text.x = element_blank(), 
                                                        axis.text.y = element_blank()),
                                 align = "vh", 
                                 axis="tblr", 
                                 ncol=2)

row_3_title <- ggdraw()+
  draw_label(
    "Positive Mood and Negative Events Models", 
    hjust = .5
  )

plot_row_3 <- cowplot::plot_grid(plot_PA_NE + theme(legend.position = "none"), 
                                 eff_plot_PA_NE + theme(legend.position = "none", 
                                                        axis.text.x = element_blank(), 
                                                        axis.text.y = element_blank()),
                                 align = "vh", 
                                 axis="tblr", 
                                 ncol=2)

row_4_title <- ggdraw()+
  draw_label(
    "Positive Mood and Positive Events Models", 
    hjust = .5
  )

plot_row_4 <- cowplot::plot_grid(plot_PA_PE + theme(legend.position = "none"), 
                                 eff_plot_PA_PE + theme(legend.position = "none", 
                                                        axis.text.y = element_blank()),
                                 align = "vh", 
                                 axis="tblr", 
                                 ncol=2)

plot_row_5 <- cowplot::plot_grid(var_legend, 
                                 eff_legend,
                                 align = "v", 
                                 axis="tb", 
                                 ncol=2)

png("~/dr-consulting_GH/shackman-umd-pax-ema-pub/plot_outputs/First_draft_summary_effects.png", 
    height = 14, width = 11, units="in", res=600)
cowplot::plot_grid(main_title_row, 
                   row_1_title,
                   plot_row_1,
                   row_2_title,
                   plot_row_2,
                   row_3_title,
                   plot_row_3,
                   row_4_title,
                   plot_row_4,
                   plot_row_5,
                   align = "vh",
                   axis = "tblr",
                   rel_heights = c(1, .5, 3, .5, 3, .5, 3, .5, 3, 1),
                   ncol = 1)
dev.off()