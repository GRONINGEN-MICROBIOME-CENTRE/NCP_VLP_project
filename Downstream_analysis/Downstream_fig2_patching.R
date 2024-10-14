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
setwd("~/Desktop/Projects_2024/AMG_paper/")
#######################################################################################################################################

## Load the  data
#######################################################################################################################################
df_for_figure2c <- as.data.frame(read_tsv('df_for_figure2c.tsv'))
df_for_figure2d <- as.data.frame(read_tsv('df_for_figure2d.tsv'))
df_for_figure2e <- as.data.frame(read_tsv('df_for_figure2e.tsv'))
df_for_figure2e$Timepoint <- factor(df_for_figure2e$Timepoint, levels = c("M0", "M1", "M2", "M3", "M4", "M6", "M12"))
df_for_figure2e <- df_for_figure2e[df_for_figure2e$Dataset == "RPKM_count", ]
df_for_figure2d$Type <- factor(df_for_figure2d$Type, levels = c("same", "different"))
#######################################################################################################################################

## Set the working directory
#######################################################################################################################################
#setwd("C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\plots\\neg_ctrl_sharing")
#######################################################################################################################################

## Patching the Figure 2
#######################################################################################################################################
# FIGURE 2a
figure_2A <- readRDS(file="burkholderia.rds")
figure_2A$labels$tag <- "a"
figure_2A

# FIGURE 2b
figure_2B <- readRDS(file="phix.rds")
figure_2B$labels$tag <- "b"

# FIGURE 2c
figure_2C <- ggplot(df_for_figure2c, aes(x = type_cohort, y = Distance )) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outlier.alpha = 0) +
  scale_y_log10() +
  facet_grid( ~ category, scales = "free", labeller = labeller(
    category = c(
      "NCs" = "NCs and NCs",
      "Samples" = "NCs and Samples"
    ))) +
  labs(y = "log10(1 - (Bray-Curtis dissimilarity))", x = "Type of cohort \n (samples from the same or different cohort compared)", tag="c") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=7),
    #plot.title = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 7), 
    axis.text.y = element_text(size = 7),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = ggtext::element_markdown(face = "bold")
  )

stat.test2c <- df_for_figure2c[df_for_figure2c$cohort_nc!="garmaeva",] %>%
  group_by(category) %>%
  t_test(Distance ~ type_cohort) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test2c <- stat.test2c %>% 
  add_xy_position(x = "type_cohort", dodge = 0.8)

stat.test2c$y.position <- log10(min(1-df_for_figure2c$Distance))

stat.test2d$xmin <- 1
stat.test2d$xmax <- 2

stat.test2d$p.adj.signif <- "***"

figure_2C <- figure_2C + 
  stat_pvalue_manual(stat.test2c, tip.length = 0.02, size=2.5, label = "p.adj.signif")
figure_2C

# FIGURE 2d
correlations <- df_for_figure2d %>%
  group_by(cohort) %>%
  summarize(cor = cor(different_cohort_NC_presence, same_cohort_NC_presence, method = "spearman"))

figure_2D <- ggplot(df_for_figure2d, aes(x=different_cohort_NC_presence, y=same_cohort_NC_presence)) +
  geom_point(size = 0.7, color="#2E236C", alpha=0.75) +  # Adjusted point size and added transparency
  geom_smooth(method="lm", color="#2E236C", fill="#C8ACD6", se=TRUE, linewidth=0.5) +  # Use linewidth instead of size
  facet_wrap(~ cohort, nrow = 2, ncol = 2, scales = "free", labeller = labeller(
    cohort = c(
      "garmaeva" = "Samples Garmaeva *et al.* <br>",
      "liang" = "Samples Liang *et al.* <br>",
      "maqsood" = "Samples Maqsood *et al.* <br>",
      "shah" = "Samples Shah *et al.* <br>"
    )
  )) +
  labs(x = "% shared vOTUs with NCs from different studies", y = "% shared vOTUs with NCs from same study", tag="d") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=7),
    #plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = ggtext::element_markdown(face = "bold")
  ) +
  # Add Spearman correlation values as text on the facets
  geom_text(data = correlations, aes(x = Inf, y = Inf, label = paste("rho = ", round(cor, 2))),
            hjust = 1.1, vjust = 1.5, size = 3)

figure_2D
# FIGURE 2e

figure_2E <- ggplot(df_for_figure2e, aes(x = Timepoint, y = value)) +
  geom_jitter(width = 0.1, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outlier.alpha = 0, width=0.5) +
  scale_y_log10() +
  facet_grid(. ~ cohort, scales = "free_x",
             space='free_x', labeller = labeller(
    cohort = c(
      "garmaeva" = "Samples Garmaeva *et al.* <br>",
      "liang" = "Samples Liang *et al.* <br>")
  )) +
  labs(y = "log10(% shared vOTUs)", tag="e") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=6),
    #plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = ggtext::element_markdown(face = "bold")
  )

stat.test2e <- df_for_figure2e %>%
  group_by(cohort) %>%
  t_test(value ~ Timepoint) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test2e <- stat.test2e %>% add_xy_position(x = "value")
stat.test2e$y.position <- log10(stat.test2e$y.position)

stat.test2e$xmin <- 1
stat.test2e$xmax[stat.test2e$cohort == "garmaeva"] <- 5
stat.test2e$xmax[stat.test2e$cohort == "liang"] <- 3
stat.test2e <- stat.test2e[c(4, 12), ]

stat.test2e$p.signif <- c("ns", "***")
figure_2E <- figure_2E + stat_pvalue_manual(stat.test2e, tip.length = 0.009, size=2.5, label = "p.signif")
figure_2E


# FIGURE 2f
figure_2F <- readRDS(file="micro_bacteroides.rds")
figure_2F$labels$tag <- "f"




# Combine the plots using patchwork
figure_2BCD <- (figure_2B | figure_2C | figure_2D) + plot_layout(widths = c(2,3,5))
figure_2EF <- (figure_2E | figure_2F) + plot_layout(widths = c(3.5,6.5))
combined_plot2 <- figure_2A / figure_2BCD / figure_2EF + plot_layout(heights  = c(3, 4, 3))

ggsave("combined_figure2.pdf", combined_plot2, width = 21/2.54, height = 24.7/2.54)
#######################################################################################################################################

