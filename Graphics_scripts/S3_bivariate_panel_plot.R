library(tidyverse)
library(ggExtra)
library(glue)
library(cowplot)

ROOT_DIR <- "~/github/ATNL/shackman-umd-pax-ema-pub"
DATA_PATH <- glue("{ROOT_DIR}/Data/study3_data.RData")
PLOT_DIR <- glue("{ROOT_DIR}/plot_outputs")
load(DATA_PATH)

source("~/R_utils/ggplot_utils/summary_functions.R")

df_threat <- dat.study3 %>% 
    group_by(ID, valence) %>% 
    summarize(
        anx_rating = mean(rating, na.rm=TRUE),
        DN = mean(DN, na.rm=TRUE)
    ) %>%
    ungroup() %>% 
    select(DN, valence, anx_rating) %>%
    filter(valence == "Threat")

df_safe <- dat.study3 %>% 
    group_by(ID, valence) %>% 
    summarize(
        anx_rating = mean(rating, na.rm=TRUE),
        DN = mean(DN, na.rm=TRUE)
    ) %>%
    ungroup() %>% 
    select(DN, valence, anx_rating) %>%
    filter(valence == "Safety")

df_comb <- dat.study3 %>% 
    group_by(ID) %>% 
    summarize(
        anx_rating = mean(rating, na.rm=TRUE),
        DN = mean(DN, na.rm=TRUE)
    ) %>%
    ungroup() %>% 
    select(DN, anx_rating)
    

threat_fit <- regression(df_threat, "DN", "anx_rating")
safe_fit <- regression(df_safe, "DN", "anx_rating")
comb_fit <- regression(df_comb, "DN", "anx_rating")

g1 <- ggplot(data = df_threat, aes(x = DN, y = anx_rating)) +
    geom_point(color = RColorBrewer::brewer.pal(9, "Reds")[5]) +
    stat_smooth(method = "lm", se = FALSE, color = RColorBrewer::brewer.pal(9, "Blues")[5], lwd = 1.5) +
    labs(y = "Anxiety Rating During 'Threat' Trials", 
         x = "Standardized DN Scores", 
         title = "") +
    scale_y_continuous(limits = c(.5, 4.5)) +
    scale_x_continuous(limits = c(-2.25, 2.25)) +
    geom_label(data = threat_fit, 
               inherit.aes = FALSE, 
               aes(x = label_pos(-2.25, 2.25, .1), 
                   y = label_pos(.5, 4.5, .9), 
                   label = paste("R^2 ==", as.character(R2))), 
               parse = TRUE) +
    theme_bw() +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(color = RColorBrewer::brewer.pal(9, "Reds")[5]))

g2 <- ggplot(data = df_safe, aes(x = DN, y = anx_rating)) +
    geom_point(color = RColorBrewer::brewer.pal(9, "Blues")[5]) +
    stat_smooth(method = "lm", se = FALSE, color = RColorBrewer::brewer.pal(9, "Reds")[5], lwd = 1.5) +
    labs(y = "Anxiety Rating During 'Safe' Trials", 
         x = "Standardized DN Scores", 
         title = "") +
    scale_y_continuous(limits = c(.5, 4.5)) +
    scale_x_continuous(limits = c(-2.25, 2.25)) +
    geom_label(data = safe_fit, 
               inherit.aes = FALSE, 
               aes(x = label_pos(-2.25, 2.25, .1), 
                   y = label_pos(.5, 4.5, .9), 
                   label = paste("R^2 ==", as.character(R2))), 
               parse = TRUE) +
    theme_bw() +
    theme(axis.title.y = element_text(color = RColorBrewer::brewer.pal(9, "Blues")[5]))

g3 <- ggplot(data = df_comb, aes(x = DN, y = anx_rating)) +
    geom_point(color = RColorBrewer::brewer.pal(9, "Purples")[5]) +
    stat_smooth(method = "lm", se = FALSE, color = RColorBrewer::brewer.pal(9, "Reds")[5], lwd = 1.5) +
    labs(y = "Mean Anxiety Ratings", 
         x = "Standardized DN Scores", 
         title = "") +
    scale_y_continuous(limits = c(.5, 4.5)) +
    scale_x_continuous(limits = c(-2.25, 2.25)) +
    geom_label(data = comb_fit, 
               inherit.aes = FALSE, 
               aes(x = label_pos(-2.25, 2.25, .1), 
                   y = label_pos(.5, 4.5, .9), 
                   label = paste("R^2 ==", as.character(R2))), 
               parse = TRUE) +
    theme_bw() +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),)

g <- plot_grid(
    ggMarginal(g3, margins = "both", type = "histogram", fill = RColorBrewer::brewer.pal(9, "Purples")[5]),
    ggMarginal(g1, margins = "y", type = "histogram", fill = RColorBrewer::brewer.pal(9, "Reds")[5]), 
    ggMarginal(g2, margins = "y", type = "histogram", fill = RColorBrewer::brewer.pal(9, "Blues")[5]), 
    ncol = 1, align = "hv", axis = "tbl", labels = c("a.", "b.", "c."), rel_heights = c(1/(5/6), 1, 1)
    )

png(filename = glue("{PLOT_DIR}/S3_tri_panel_bivariate_col.png"), units = "in", height = 13.5, width = 4.5, res = 900)
g
dev.off()

postscript(glue("{PLOT_DIR}/S3_tri_panel_bivariate_col.eps"), height = 13.5, width = 4.5)
g
dev.off()