## Code description
#######################################################################################################################################
## Script for Patching the Figure 2
## 
#######################################################################################################################################

## Load libraries
#######################################################################################################################################
library(readr)
library(vegan)
library(dplyr)
library(stringr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggtext)
library(patchwork)
library(ggtree)
library(lmerTest)
library(rstatix)
library(ggpubr)
#######################################################################################################################################

## Set the working directory
#######################################################################################################################################
setwd("C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\data")
#######################################################################################################################################

## Load the  data
#######################################################################################################################################
df_for_figure2c <- as.data.frame(read_tsv('df_for_figure2c.tsv'))
df_for_figure2d <- as.data.frame(read_tsv('df_for_figure2d.tsv'))
df_for_figure2e <- as.data.frame(read_tsv('df_for_figure2e.tsv'))
df_for_figure2e$Timepoint <- factor(df_for_figure2e$Timepoint, levels = c("M0", "M1", "M2", "M3", "M4", "M6", "M12"))
df_for_figure2e <- df_for_figure2e[df_for_figure2e$Dataset == "RPKM_count", ]

df_for_figure2d <- df_for_figure2d %>%
  mutate(number_NCs = ifelse(cohort %in% c("maqsood", "shah"), 8, NA),
         number_NCs = ifelse(cohort %in% c("garmaeva"), 1, number_NCs),
         number_NCs = ifelse(cohort %in% c("liang"), 20, number_NCs))

ggplot(df_for_figure2d, aes(x=number_NCs, y=same_cohort_NC_presence)) +
  geom_point() +  # Adjusted point size and added transparency
  geom_smooth() # Use linewidth instead of size
  
summary(lmer(same_cohort_NC_presence ~ number_NCs + (1|cohort/nc_subject_group), REML = F, data = df_for_figure2d))
summary(lmer(same_cohort_NC_presence ~ number_NCs + Type + Timepoint + (1|nc_subject_group), REML = F, data = df_for_figure2d))
df_for_figure2d$Type <- factor(df_for_figure2d$Type, levels = c("same", "different"))
df_for_figure2d$Timepoint_numeric
#######################################################################################################################################

## Set the working directory
#######################################################################################################################################
setwd("C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\plots\\neg_ctrl_sharing")
#######################################################################################################################################

## Patching the Figure 2
#######################################################################################################################################
figure_2A <- readRDS(file="burkholderia.rds")
figure_2A$labels$tag <- "a"

figure_2B <- readRDS(file="phix.rds")
figure_2B$labels$tag <- "b"

figure_2C <- readRDS(file="micro_bacteroides.rds")
figure_2C$labels$tag <- "c"

figure_2D <- ggplot(df_for_figure2c, aes(x = type_cohort, y = Distance)) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outlier.alpha = 0) +
  # scale_y_log10() +
  facet_grid(cohort_nc ~ category, scales = "free", labeller = labeller(
    category = c(
      "NCs" = "NCs and NCs",
      "Samples" = "NCs and Samples"
    ),
    cohort_nc = c(
      "garmaeva" = "Samples Garmaeva *et al.* <br>",
      "liang" = "Samples Liang *et al.* <br>",
      "maqsood" = "Samples Maqsood *et al.* <br>",
      "shah" = "Samples Shah *et al.* <br>"
    )
  )) +
  labs(y = "log10(1 - (Bray-Curtis dissimilarity))", x = "Type of cohort \n (samples from the same or different cohort compared)", tag="d") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    plot.title = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 8), 
    axis.text.y = element_text(size = 6),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = ggtext::element_markdown(face = "bold")
  )

stat.test2d <- df_for_figure2c[df_for_figure2c$cohort_nc!="garmaeva",] %>%
  group_by(cohort_nc, category) %>%
  t_test(Distance ~ type_cohort) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test2d <- stat.test2d %>% 
  add_xy_position(x = "type_cohort", dodge = 0.8)

stat.test2d[7,] <- stat.test2d[3,]
stat.test2d[8,] <- stat.test2d[4,]
stat.test2d[7,"cohort_nc"] <- "garmaeva"
stat.test2d[8,"cohort_nc"] <- "garmaeva"
stat.test2d[7,"y.position"] <- log(max(df_for_figure2c[df_for_figure2c$cohort_nc=="garmaeva" & df_for_figure2c$category == "NCs",]$Distance + 1))
stat.test2d[8,"y.position"] <- log(max(df_for_figure2c[df_for_figure2c$cohort_nc=="garmaeva" & df_for_figure2c$category == "Samples",]$Distance + 1))

stat.test2d$xmin <- 1
stat.test2d$xmax <- 2

stat.test2d$p.signif <- c("****", "****", "****", "****", "****", "****", "NA", "****")

figure_2D <- figure_2D + 
  stat_pvalue_manual(stat.test2d, tip.length = 0.02, size=2.5, label = "p.signif")
figure_2D

correlations <- df_for_figure2d %>%
  group_by(cohort) %>%
  summarize(cor = cor(different_cohort_NC_presence, same_cohort_NC_presence, method = "spearman"))

figure_2E <- ggplot(df_for_figure2d, aes(x=different_cohort_NC_presence, y=same_cohort_NC_presence)) +
  geom_point(size = 0.5, color="#2E236C", alpha=0.75) +  # Adjusted point size and added transparency
  geom_smooth(method=lm, color="#2E236C", fill="#C8ACD6", se=TRUE, linewidth=0.5) +  # Use linewidth instead of size
  facet_wrap(~ cohort, nrow = 2, ncol = 2, scales = "free", labeller = labeller(
    cohort = c(
      "garmaeva" = "Samples Garmaeva *et al.* <br>",
      "liang" = "Samples Liang *et al.* <br>",
      "maqsood" = "Samples Maqsood *et al.* <br>",
      "shah" = "Samples Shah *et al.* <br>"
    )
  )) +
  labs(x = "% shared vOTUs with NCs \n from different studies", y = "% shared vOTUs with NCs \n from same study", tag="e") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = ggtext::element_markdown(face = "bold")
  ) +
  # Add Spearman correlation values as text on the facets
  geom_text(data = correlations, aes(x = Inf, y = Inf, label = paste("ρ = ", round(cor, 2))),
            hjust = 1.1, vjust = 1.5, size = 3)

figure_2E

figure_2F <- ggplot(df_for_figure2e, aes(x = Timepoint, y = value)) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outlier.alpha = 0) +
  scale_y_log10() +
  facet_grid(. ~ cohort, scales = "free", labeller = labeller(
    cohort = c(
      "garmaeva" = "Samples Garmaeva *et al.* <br>",
      "liang" = "Samples Liang *et al.* <br>"
    )
  )) +
  labs(y = "log10(% shared vOTUs)", tag="f") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = ggtext::element_markdown(face = "bold")
  )

stat.test2f <- df_for_figure2e %>%
  group_by(cohort) %>%
  t_test(value ~ Timepoint) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test2f <- stat.test2f %>% add_xy_position(x = "value")
stat.test2f$y.position <- log10(stat.test2f$y.position)

stat.test2f$xmin <- 1
stat.test2f$xmax[stat.test2f$cohort == "garmaeva"] <- 5
stat.test2f$xmax[stat.test2f$cohort == "liang"] <- 3
stat.test2f <- stat.test2f[c(4, 12), ]

stat.test2f$p.signif <- c("ns", "***")
figure_2F <- figure_2F + stat_pvalue_manual(stat.test2f, tip.length = 0, size=2.5, label = "p.signif")
figure_2F


# Combine the plots using patchwork
combined_plot2 <- (figure_2A + figure_2B + figure_2C + plot_layout(nrow=1, guides = "collect")) / (figure_2D + figure_2E + figure_2F) +
  plot_layout(heights = c(2.2, 1))# Save the combined plot as a PDF

ggsave("combined_figure2.png", combined_plot2, width = 45/2.54, height = 30/2.54)
#######################################################################################################################################

