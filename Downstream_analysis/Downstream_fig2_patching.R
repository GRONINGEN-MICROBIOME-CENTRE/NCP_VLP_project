## Code description
#######################################################################################################################################
## Script for Patching the Figure 2
## 
#######################################################################################################################################

## Load libraries
#######################################################################################################################################
library(ggplot2)
library(ComplexUpset)
library(readr)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(patchwork)

# library(vegan)
# library(dplyr)
# library(stringr)
# library(reshape2)
# library(tidyr)
# library(ggplot2)
# library(ggrepel)
# library(ggtext)
# library(patchwork)
# library(ggtree)
# library(lmerTest)
# 
# 

#######################################################################################################################################

## Set the working directory
#######################################################################################################################################
setwd("~/Desktop/Projects_2024/AMG_paper/")
#######################################################################################################################################

## Load the  data
#######################################################################################################################################
nc_sharing <- read.table('data_upset_plot_NCs_sharing.txt', sep='\t', header=T, check.names = F)
studies <- colnames(nc_sharing)

df_for_figure2c <- as.data.frame(read_tsv('df_for_figure2c.tsv'))
#######################################################################################################################################

## Set the working directory
#######################################################################################################################################
#setwd("C:\\Users\\Natal\\Documents\\UMCG\\amg_paper\\plots\\neg_ctrl_sharing")
#######################################################################################################################################

## Patching the Figure 2
#######################################################################################################################################
# FIGURE 2a, UpSet plot:
tada <- colnames(nc_sharing)
F2A <- upset(nc_sharing, studies, name = NULL, 
      set_sizes=(
  upset_set_size(mapping=aes(fill='bars_color'))
  + ylab('N vOTUs in NCs') +  
    scale_fill_manual(values=c('bars_color'='#8174A0'), guide='none' ) )
    ) +
  labs(tag="a") +
  theme(
    plot.tag = element_text(face = "bold"), 
    plot.tag.position = c(-0.37, 2.5),
    axis.text.y = ggtext::element_markdown()
  )

#
figure_2B <- readRDS(file="burkholderia.rds")
figure_2B$labels$tag <- "b"
figure_2B$labels$title <- "Burkholderia phage L41225"


# FIGURE 2b
figure_2C <- readRDS(file="phix.rds")
figure_2C$labels$tag <- "c"

# FIGURE 2c
figure_2D <- ggplot(df_for_figure2c, aes(x = type_cohort, y = Distance )) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outlier.alpha = 0) +
  scale_y_log10() +
  facet_grid( ~ category, scales = "free", labeller = labeller(
    category = c(
      "NCs" = "NCs vs NCs",
      "Samples" = "NCs vs Samples"
    ))) +
  labs(y = "log<sub>10</sub>Similarity index", x = "study", tag="d") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=8),
    #plot.title = element_text(size = 10),
    axis.title.x = element_text(size=8),
    axis.title.y = ggtext::element_markdown(size = 9),
    axis.text.x = element_text(size = 8), 
    axis.text.y = element_text(size = 8),
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

stat.test2c$y.position <- 1

stat.test2c$p.adj.signif <- "***"

figure_2D <- figure_2D + 
  stat_pvalue_manual(stat.test2c, tip.length = 0.02, size=2.5, label = "p.adj.signif")
figure_2D

# Combine the plots using patchwork
f2AB <- (F2A / figure_2B) + plot_layout(widths = c(4,6))
f2CD <- (figure_2C | figure_2D) + plot_layout(widths = c(5, 5))

combined_plot2 <- f2AB /  f2CD + plot_layout(heights  = c(3, 3, 4))


ggsave("combined_figure2_UPD2.pdf", combined_plot2, width = 21/2.54, height = 27/2.54)
#######################################################################################################################################

