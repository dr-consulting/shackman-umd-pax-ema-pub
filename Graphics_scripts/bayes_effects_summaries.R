############################ NOT USED IN THE MANUSCRIPT AND NOT FULLY VETTED ##########################################
# Ultimately decided on a simpler summary table approach. 


source("~/dr-consulting_GH/shackman-umd-pax-ema-pub/Graphics_scripts/bayes_variance_riverplots_utils.R")
library(tidyverse)
library(tidybayes)
library(cowplot)

GRAPHICS_DIR <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/plot_outputs/"

green_check_mark <- ggplot(data = data.frame(x=.5, y=.5, label = emojifont::emoji("heavy_check_mark")), 
                           aes(x=x, y=y, label=label)) +
  geom_text(family="OpenSansEmoji", color = RColorBrewer::brewer.pal(9, "Greens")[7], size=20)+
  theme_void()

red_cross_mark <- ggplot(data = data.frame(x=.5, y=.5, label = emojifont::emoji("x")), 
                         aes(x=x, y=y, label=label)) +
  geom_text(family="OpenSansEmoji", color = RColorBrewer::brewer.pal(9, "Reds")[7], size=20)+
  theme_void()

# Study 1
POSTERIOR_PATH <- "/media/matthew/My Book/EMA_S1_Bayesian_Posteriors"
DATA_FILEPATH <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Data/study1_data.RData"
SUMMARY_DIRPATH <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Study_1_model_summaries/"
data_loader(POSTERIOR_PATH, "S1_NEG_NegEvnt_Rct", DATA_FILEPATH)

DN <- posterior_samples(S1_NEG_NegEvnt_Rct, pars="b_NEG_c.DN")
Agg_NegEvnt <- posterior_samples(S1_NEG_NegEvnt_Rct, pars="b_NEG_m.NegEvnt")
NegEvnt <- posterior_samples(S1_NEG_NegEvnt_Rct, pars="bsp_NEG_mic.NegEvnt")[1]
NegEvnt_react <- posterior_samples(S1_NEG_NegEvnt_Rct, pars="bsp_NEG_mic.NegEvnt:c.DN")

S1_NA_NE <- data.frame(DN, 
                       Agg_NegEvnt, 
                       NegEvnt, 
                       NegEvnt_react)

names(S1_NA_NE) <- c("DN", "Mean NE Rating", "NE Rating", "DN x NE Reactivity")

# Setting fill color
fill_clr <- "grey50"
if(0 < quantile(S1_NA_NE[,"DN"], .025)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Reds")[7]
}
if(0 > quantile(S1_NA_NE[,"DN"], .975)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Blues")[7]
}


S1_DN <- ggplot(data = S1_NA_NE["DN"], aes(x=DN, y=0, fill="DN"))+
  stat_halfeyeh(show.legend=FALSE, alpha=.7, adjust=1.75)+
  geom_vline(xintercept = 0, color="black", lty="dashed")+
  scale_fill_manual(values = c(fill_clr))+
  scale_x_continuous(breaks = c(0, .075, .15))+
  theme_bw()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank())

S1_DN

# Setting fill color
fill_clr <- "grey50"
if(0 < quantile(S1_NA_NE[,"Mean NE Rating"], .025)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Reds")[7]
}
if(0 > quantile(S1_NA_NE[,"Mean NE Rating"], .975)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Blues")[7]
}

S1_Agg_NegEvnt <- ggplot(data = S1_NA_NE["Mean NE Rating"], aes(x=`Mean NE Rating`, y=0, fill="Mean NE Rating"))+
  stat_halfeyeh(show.legend=FALSE, alpha=.7, adjust=1.75)+
  geom_vline(xintercept = 0, color="black", lty="dashed")+
  scale_fill_manual(values = c(fill_clr))+
  scale_x_continuous(breaks = c(0, .15, .3))+
  theme_bw()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank())

S1_Agg_NegEvnt

# Setting fill color
fill_clr <- "grey50"
if(0 < quantile(S1_NA_NE[,"NE Rating"], .025)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Reds")[7]
}
if(0 > quantile(S1_NA_NE[,"NE Rating"], .975)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Blues")[7]
}

S1_NegEvnt <- ggplot(data = S1_NA_NE["NE Rating"], aes(x=`NE Rating`, y=0, fill="NE Rating"))+
  stat_halfeyeh(show.legend=FALSE, alpha=.7, adjust=1.75)+
  geom_vline(xintercept = 0, color="black", lty="dashed")+
  scale_fill_manual(values = c(fill_clr))+
  scale_x_continuous(breaks = c(0, .05, .1))+
  theme_bw()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank())

S1_NegEvnt

# Setting fill color
fill_clr <- "grey50"
if(0 < quantile(S1_NA_NE[,"DN x NE Reactivity"], .025)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Reds")[7]
}
if(0 > quantile(S1_NA_NE[,"DN x NE Reactivity"], .975)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Blues")[7]
}

S1_NegEvnt_react  <- ggplot(data = S1_NA_NE["DN x NE Reactivity"], aes(x=`DN x NE Reactivity`, y=0, 
                                                                       fill="DN x NE Reactivity"))+
  stat_halfeyeh(show.legend=FALSE, alpha=.7, adjust=1.75)+
  geom_vline(xintercept = 0, color="black", lty="dashed")+
  scale_fill_manual(values = c(fill_clr))+
  scale_x_continuous(limits = c(-0.02, 0.09), breaks = c(0, .04, .08))+
  theme_bw()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank())

S1_NegEvnt_react 


# Study 2 
POSTERIOR_PATH <- "/media/matthew/My Book/EMA_S2_Bayesian_Posteriors"
DATA_FILEPATH <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Data/study2_data.RData"
SUMMARY_DIR <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Study_2_model_summaries/"
data_loader(POSTERIOR_PATH, "S2_NEG_NegEvnt_x_DN_prop.NegEvnt", DATA_FILEPATH)

DN <- posterior_samples(S2_NEG_NegEvnt_x_DN_prop.NegEvnt, pars="b_c.DN")
Agg_NegEvnt <- posterior_samples(S2_NEG_NegEvnt_x_DN_prop.NegEvnt, pars="b_prop.NegEvnt")
NegEvnt <- posterior_samples(S2_NEG_NegEvnt_x_DN_prop.NegEvnt, pars="b_c.NegEvnt")[1]
NegEvnt_react <- posterior_samples(S2_NEG_NegEvnt_x_DN_prop.NegEvnt, pars="b_c.NegEvnt:c.DN")

S2_NA_NE <- data.frame(DN, 
                       Agg_NegEvnt, 
                       NegEvnt, 
                       NegEvnt_react)

names(S2_NA_NE) <- c("DN", "NE Relative Frequency", "NE", "DN x NE Reactivity")

# Setting fill color
fill_clr <- "grey50"
if(0 < quantile(S2_NA_NE[,"DN"], .025)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Reds")[7]
}
if(0 > quantile(S2_NA_NE[,"DN"], .975)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Blues")[7]
}

S2_DN <- ggplot(data = S2_NA_NE["DN"], aes(x=DN, y=0, fill="DN"))+
  stat_halfeyeh(show.legend=FALSE, alpha=.7, adjust=1.75)+
  geom_vline(xintercept = 0, color="black", lty="dashed")+
  scale_fill_manual(values = c(fill_clr))+
  scale_x_continuous(breaks = c(0, .06, .12))+
  theme_bw()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank())

S2_DN

# Setting fill color
fill_clr <- "grey50"
if(0 < quantile(S2_NA_NE[,"NE Relative Frequency"], .025)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Reds")[7]
}
if(0 > quantile(S2_NA_NE[,"NE Relative Frequency"], .975)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Blues")[7]
}

S2_Agg_NegEvnt <- ggplot(data = S2_NA_NE["NE Relative Frequency"], aes(x=`NE Relative Frequency`, y=0, 
                                                                       fill="NE Relative Frequency"))+
  stat_halfeyeh(show.legend=FALSE, alpha=.7, adjust=1.75)+
  geom_vline(xintercept = 0, color="black", lty="dashed")+
  scale_fill_manual(values = c(fill_clr))+
  scale_x_continuous(limits = c(NA, 1.30), breaks=c(0, .5, 1))+
  theme_bw()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank())

S2_Agg_NegEvnt

# Setting fill color
fill_clr <- "grey50"
if(0 < quantile(S2_NA_NE[,"NE"], .025)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Reds")[7]
}
if(0 > quantile(S2_NA_NE[,"NE"], .975)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Blues")[7]
}

S2_NegEvnt <- ggplot(data = S2_NA_NE["NE"], aes(x=`NE`, y=0, fill="NE"))+
  stat_halfeyeh(show.legend=FALSE, alpha=.7, adjust=1.75)+
  geom_vline(xintercept = 0, color="black", lty="dashed")+
  scale_fill_manual(values = c(fill_clr))+
  scale_x_continuous(breaks = c(0, .10, .20))+
  theme_bw()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank())

S2_NegEvnt

# Setting fill color
fill_clr <- "grey50"
if(0 < quantile(S2_NA_NE[,"DN x NE Reactivity"], .025)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Reds")[7]
}
if(0 > quantile(S2_NA_NE[,"DN x NE Reactivity"], .975)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Blues")[7]
}

S2_NegEvnt_react  <- ggplot(data = S2_NA_NE["DN x NE Reactivity"], aes(x=`DN x NE Reactivity`, y=0, 
                                                                       fill="DN x NE Reactivity"))+
  stat_halfeyeh(show.legend=FALSE, alpha=.7, adjust=1.75)+
  geom_vline(xintercept = 0, color="black", lty="dashed")+
  scale_fill_manual(values = c(fill_clr))+
  scale_x_continuous(limits = c(-0.02, 0.09), breaks = c(0, .04, .08))+
  theme_bw()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank())

S2_NegEvnt_react 

row_1_title <- ggdraw() +
  draw_label(
    "Study 1: Posterior Distribution of Negative Affect/Negative Events Model Parameters", 
    hjust = .5
  )

row_3_title <- ggdraw() +
  draw_label(
    "Study 2: Posterior Distribution of Negative Affect/Negative Events Model Parameters", 
    hjust = .5
    )

S1_grid <- plot_grid(S1_DN, S1_Agg_NegEvnt, S1_NegEvnt, S1_NegEvnt_react, 
                     ncol=4)
S2_grid <- plot_grid(S2_DN, S2_Agg_NegEvnt, S2_NegEvnt, S2_NegEvnt_react, 
                     ncol=4)

compare_row <- plot_grid(green_check_mark, green_check_mark, green_check_mark, red_cross_mark, 
                         ncol = 4)

NA_NE_plot <- plot_grid(row_1_title, 
                        S1_grid,
                        row_3_title,
                        S2_grid,
                        nrow = 4, align = "hv", axis = "tblr", 
                        rel_heights = c(.2, 1, .2, 1))


NA_NE_plot<- plot_grid(NA_NE_plot, 
                       compare_row, 
                       nrow = 2, 
                       rel_heights = c(2.4,1))

png(paste0(GRAPHICS_DIR, "S1_S2_NEG_NegEvnt_Model_Effects_Summary.png"), 
    units="in", height = 5, width=15, res=600)
NA_NE_plot
dev.off()

#######################################################################################################################
# Study 1
POSTERIOR_PATH <- "/media/matthew/My Book/EMA_S1_Bayesian_Posteriors"
DATA_FILEPATH <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Data/study1_data.RData"
SUMMARY_DIRPATH <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Study_1_model_summaries/"
data_loader(POSTERIOR_PATH, "S1_NEG_PosEvnt_Rct", DATA_FILEPATH)

DN <- posterior_samples(S1_NEG_PosEvnt_Rct, pars="b_NEG_c.DN")
Agg_PosEvnt <- posterior_samples(S1_NEG_PosEvnt_Rct, pars="b_NEG_m.PosEvnt")
PosEvnt <- posterior_samples(S1_NEG_PosEvnt_Rct, pars="bsp_NEG_mic.PosEvnt")[1]
PosEvnt_react <- posterior_samples(S1_NEG_PosEvnt_Rct, pars="bsp_NEG_mic.PosEvnt:c.DN")

S1_NA_PE <- data.frame(DN, 
                       Agg_PosEvnt, 
                       PosEvnt, 
                       PosEvnt_react)

names(S1_NA_PE) <- c("DN", "Mean PE Rating", "PE Rating", "DN x PE Reactivity")

# Setting fill color
fill_clr <- "grey50"
if(0 < quantile(S1_NA_PE[,"DN"], .025)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Reds")[7]
}
if(0 > quantile(S1_NA_PE[,"DN"], .975)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Blues")[7]
}

S1_DN <- ggplot(data = S1_NA_PE["DN"], aes(x=DN, y=0, fill="DN"))+
  stat_halfeyeh(show.legend=FALSE, alpha=.7, adjust=1.75)+
  geom_vline(xintercept = 0, color="black", lty="dashed")+
  scale_fill_manual(values = c(fill_clr))+
  scale_x_continuous(breaks = c(0, .1, .2))+
  theme_bw()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank())

S1_DN

# Setting fill color
fill_clr <- "grey50"
if(0 < quantile(S1_NA_PE[,"Mean PE Rating"], .025)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Reds")[7]
}
if(0 > quantile(S1_NA_PE[,"Mean PE Rating"], .975)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Blues")[7]
}

S1_Agg_PosEvnt <- ggplot(data = S1_NA_PE["Mean PE Rating"], aes(x=`Mean PE Rating`, y=0, fill="Mean PE Rating"))+
  stat_halfeyeh(show.legend=FALSE, alpha=.7, adjust=1.75)+
  geom_vline(xintercept = 0, color="black", lty="dashed")+
  scale_fill_manual(values = c(fill_clr))+
  scale_x_continuous(limits = c(-.25, .31), breaks = c(-.15, 0, .15, .30))+
  theme_bw()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank())

S1_Agg_PosEvnt

# Setting fill color
fill_clr <- "grey50"
if(0 < quantile(S1_NA_PE[,"PE Rating"], .025)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Reds")[7]
}
if(0 > quantile(S1_NA_PE[,"PE Rating"], .975)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Blues")[7]
}

S1_PosEvnt <- ggplot(data = S1_NA_PE["PE Rating"], aes(x=`PE Rating`, y=0, fill="PE Rating"))+
  stat_halfeyeh(show.legend=FALSE, alpha=.7, adjust=1.75)+
  geom_vline(xintercept = 0, color="black", lty="dashed")+
  scale_fill_manual(values = c(fill_clr))+
  scale_x_continuous(breaks = c(-.08, -.04, 0))+
  theme_bw()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank())

S1_PosEvnt

# Setting fill color
fill_clr <- "grey50"
if(0 < quantile(S1_NA_PE[,"DN x PE Reactivity"], .025)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Reds")[7]
}
if(0 > quantile(S1_NA_PE[,"DN x PE Reactivity"], .975)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Blues")[7]
}


S1_PosEvnt_react  <- ggplot(data = S1_NA_PE["DN x PE Reactivity"], aes(x=`DN x PE Reactivity`, y=0, 
                                                                       fill="DN x PE Reactivity"))+
  stat_halfeyeh(show.legend=FALSE, alpha=.7, adjust=1.75)+
  geom_vline(xintercept = 0, color="black", lty="dashed")+
  scale_fill_manual(values = c(fill_clr))+
  scale_x_continuous(limits = c(-.051, .021), breaks=c(-.04, -.02, 0, .02))+
  theme_bw()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank())

S1_PosEvnt_react 


# Study 2 
POSTERIOR_PATH <- "/media/matthew/My Book/EMA_S2_Bayesian_Posteriors"
DATA_FILEPATH <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Data/study2_data.RData"
SUMMARY_DIR <- "~/dr-consulting_GH/shackman-umd-pax-ema-pub/Study_2_model_summaries/"
data_loader(POSTERIOR_PATH, "S2_NEG_PosEvnt_x_DN_prop.PosEvnt", DATA_FILEPATH)

DN <- posterior_samples(S2_NEG_PosEvnt_x_DN_prop.PosEvnt, pars="b_c.DN")
Agg_PosEvnt <- posterior_samples(S2_NEG_PosEvnt_x_DN_prop.PosEvnt, pars="b_prop.PosEvnt")
PosEvnt <- posterior_samples(S2_NEG_PosEvnt_x_DN_prop.PosEvnt, pars="b_c.PosEvnt")[1]
PosEvnt_react <- posterior_samples(S2_NEG_PosEvnt_x_DN_prop.PosEvnt, pars="b_c.PosEvnt:c.DN")

S2_NA_PE <- data.frame(DN, 
                       Agg_PosEvnt, 
                       PosEvnt, 
                       PosEvnt_react)

names(S2_NA_PE) <- c("DN", "PE Relative Frequency", "PE", "DN x PE Reactivity")

# Setting fill color
fill_clr <- "grey50"
if(0 < quantile(S2_NA_PE[,"DN"], .025)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Reds")[7]
}
if(0 > quantile(S2_NA_PE[,"DN"], .975)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Blues")[7]
}


S2_DN <- ggplot(data = S2_NA_PE["DN"], aes(x=DN, y=0, fill="DN"))+
  stat_halfeyeh(show.legend=FALSE, alpha=.7, adjust=1.75)+
  geom_vline(xintercept = 0, color="black", lty="dashed")+
  scale_fill_manual(values = c(fill_clr))+
  scale_x_continuous(breaks = c(0, .05, .10))+
  theme_bw()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank())

S2_DN

# Setting fill color
fill_clr <- "grey50"
if(0 < quantile(S2_NA_PE[,"PE Relative Frequency"], .025)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Reds")[7]
}
if(0 > quantile(S2_NA_PE[,"PE Relative Frequency"], .975)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Blues")[7]
}


S2_Agg_PosEvnt <- ggplot(data = S2_NA_PE["PE Relative Frequency"], aes(x=`PE Relative Frequency`, y=0, 
                                                                       fill="PE Relative Frequency"))+
  stat_halfeyeh(show.legend=FALSE, alpha=.7, adjust=1.75)+
  geom_vline(xintercept = 0, color="black", lty="dashed")+
  scale_fill_manual(values = c(fill_clr))+
  scale_x_continuous(limits = c(-.25, .31), breaks=c(-.15, 0, .15, .3))+
  theme_bw()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank())

S2_Agg_PosEvnt

# Setting fill color
fill_clr <- "grey50"
if(0 < quantile(S2_NA_PE[,"PE"], .025)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Reds")[7]
}
if(0 > quantile(S2_NA_PE[,"PE"], .975)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Blues")[7]
}

S2_PosEvnt <- ggplot(data = S2_NA_PE["PE"], aes(x=`PE`, y=0, fill="PE"))+
  stat_halfeyeh(show.legend=FALSE, alpha=.7, adjust=1.75)+
  geom_vline(xintercept = 0, color="black", lty="dashed")+
  scale_fill_manual(values = c(fill_clr))+
  scale_x_continuous(breaks = c(-0.08, -0.04, 0))+
  theme_bw()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank())

S2_PosEvnt

# Setting fill color
fill_clr <- "grey50"
if(0 < quantile(S2_NA_PE[,"DN x PE Reactivity"], .025)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Reds")[7]
}
if(0 > quantile(S2_NA_PE[,"DN x PE Reactivity"], .975)){ 
  fill_clr <- RColorBrewer::brewer.pal(9, "Blues")[7]
}

S2_PosEvnt_react  <- ggplot(data = S2_NA_PE["DN x PE Reactivity"], aes(x=`DN x PE Reactivity`, y=0, fill="DN x PE Reactivity"))+
  stat_halfeyeh(show.legend=FALSE, alpha=.7, adjust=1.75)+
  geom_vline(xintercept = 0, color="black", lty="dashed")+
  scale_fill_manual(values = c(fill_clr))+
  scale_x_continuous(limits = c(-.051, .021), breaks=c(-.04, -.02, 0, .02))+
  theme_bw()+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank())

S2_PosEvnt_react 

row_1_title <- ggdraw() +
  draw_label(
    "Study 1: Posterior Distribution of Negative Affect/Positive Events Model Parameters", 
    hjust = .5
  )

row_3_title <- ggdraw() +
  draw_label(
    "Study 2: Posterior Distribution of Negative Affect/Positive Events Model Parameters", 
    hjust = .5
  )

S1_grid <- plot_grid(S1_DN, S1_Agg_PosEvnt, S1_PosEvnt, S1_PosEvnt_react, 
                     ncol=4)
S2_grid <- plot_grid(S2_DN, S2_Agg_PosEvnt, S2_PosEvnt, S2_PosEvnt_react, 
                     ncol=4)
compare_row <- plot_grid(green_check_mark, green_check_mark, green_check_mark, green_check_mark, 
                         ncol = 4)

NA_PE_plot <- plot_grid(row_1_title, 
                        S1_grid,
                        row_3_title,
                        S2_grid,
                        nrow = 4, align = "hv", axis = "tblr", 
                        rel_heights = c(.2, 1, .2, 1))

NA_PE_plot<- plot_grid(NA_PE_plot, 
                       compare_row, 
                       nrow = 2, 
                       rel_heights = c(2.4,1))

png(paste0(GRAPHICS_DIR, "S1_S2_NEG_PosEvnt_Model_Effects_Summary.png"), 
    units="in", height = 5, width=15, res=600)
NA_PE_plot
dev.off()

png(paste0(GRAPHICS_DIR, "S1_S2_NEG_Model_Effects_Summary.png"), 
    units="in", height = 5, width=15, res=600)
plot_grid(NA_NE_plot, 
          NA_PE_plot, 
          ncol=2)
dev.off()