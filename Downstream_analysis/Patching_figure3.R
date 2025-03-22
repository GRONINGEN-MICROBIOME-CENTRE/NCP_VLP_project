## Code description
#######################################################################################################################################
## Script for Patching the Figure 3
## 
#######################################################################################################################################

## Load libraries
#######################################################################################################################################
library(readr)
library(tidyverse)
library(ggpubr)
library(rstatix)
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
df_for_figure3abe <- as.data.frame(read_tsv('df_figure3abe.tsv'))

df_for_figure3c <- as.data.frame(read_tsv('df_for_figure2e.tsv'))
df_for_figure3c$Timepoint <- factor(df_for_figure3c$Timepoint, levels = c("M0", "M1", "M2", "M3", "M4", "M6", "M12"))
df_for_figure3c <- df_for_figure3c[df_for_figure3c$Dataset == "RPKM_count", ]

infant_metadata <- read.table('data_stacked_barplot_infants_strain_contamination.txt', sep='\t', header=T)
infant_metadata <- infant_metadata[ !is.na(infant_metadata$richness) & infant_metadata$richness >0,]
#######################################################################################################################################

## Patching the Figure 3
#######################################################################################################################################
# FIGURE 3a
# FIGURE 3b
stat.test3a <- df_for_figure3abe %>%
  group_by(cohort) %>%
  t_test(same_cohort_NC_presence ~ Type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test3a <- stat.test3a %>% add_xy_position(x = "same_cohort_NC_presence")

stat.test3a$xmin <- 1
stat.test3a$xmax <- 2

stat.test3a$p.signif <- c("***", "***")

figure_3A <- ggplot(df_for_figure3abe, aes(x = Type, y = same_cohort_NC_presence)) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outlier.alpha = 0) +
  # scale_y_log10() +
  facet_grid(. ~ cohort, scales = "free", labeller = labeller(
    cohort = c(
      "garmaeva" = "Samples<br>Garmaeva *et al.*",
      "maqsood" = "Samples<br>Maqsood *et al.*"
    )
  )) +
  labs(y = "% vOTUs shared with NCs", tag="a") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=8),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=13)
  )+ 
  stat_pvalue_manual(stat.test3a, tip.length = 0.02, size=2.5, label = "p.signif")
figure_3A


stat.test3b <- df_for_figure3abe %>%
  group_by(cohort) %>%
  t_test(same_cohort_NC ~ Type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test3b <- stat.test3b %>% add_xy_position(x = "same_cohort_NC")
stat.test3b$xmin <- 1
stat.test3b$xmax <- 2
stat.test3b$p.signif <- c("***", "***")

figure_3B <- ggplot(df_for_figure3abe, aes(x = Type, y = same_cohort_NC)) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outlier.alpha = 0) +
  facet_grid(. ~ cohort, scales = "free", labeller = labeller(
    cohort = c(
      "garmaeva" = "Samples<br>Garmaeva *et al.*",
      "maqsood" = "Samples<br>Maqsood *et al.*"
    )
  )) +
  labs(y = "N vOTUs shared with NCs", tag="b") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=8),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=13)
  ) + 
  stat_pvalue_manual(stat.test3b, tip.length = 0.02, size=2.5, label = "p.signif")
figure_3B


#  FIGURE 3c
stat.test3c <- df_for_figure3c %>%
  group_by(cohort) %>%
  t_test(value ~ Timepoint) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test3c <- stat.test3c %>% add_xy_position(x = "value")

stat.test3c$xmin <- 1
stat.test3c$xmax[stat.test3c$cohort == "garmaeva"] <- 5
stat.test3c$xmax[stat.test3c$cohort == "liang"] <- 3
stat.test3c <- stat.test3c[c(4, 12), ]
stat.test3c$y.position <- c(30, 105)

stat.test3c$p.signif <- c("ns", "***")

figure_3C <- ggplot(df_for_figure3c, aes(x = Timepoint, y = value)) +
  geom_jitter(width = 0.1, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outlier.alpha = 0, width=0.5) +
  # scale_y_log10() +
  facet_grid(. ~ cohort, scales = "free_x",
             space='free_x', labeller = labeller(
               cohort = c(
                 "garmaeva" = "Samples<br>Garmaeva *et al.*",
                 "liang" = "Samples<br>Liang *et al.*")
             )) +
  labs(y = "% vOTUs shared with NCs", tag="c") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=8),
    #plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = ggtext::element_markdown(face = "bold")
  )+ 
  stat_pvalue_manual(stat.test3c, tip.length = 0.009, size=2.5, label = "p.signif")

#  FIGURE 3e
sum(infant_metadata$shared_st_iNC_perc >= 10, na.rm = T)

dim(infant_metadata) # => 12%

### <debugged till here>

infant_metadata$study_br <- dplyr::recode(infant_metadata$Study,
                                          'garmaeva' = 'Garmaeva *et al*., 2024',
                                          'liang' = 'Liang *et al*., 2020',
                                          'maqsood' = 'Maqsood *et al*., 2019',
                                          'shah' = 'Shah *et al*., 2023')

# richness

infant_metadata$Contamination_Category <- cut(infant_metadata$shared_st_iNC_perc,
                                            breaks = c(-Inf, 0, 10, 20, Inf),
                                            labels = c('0', '=<10%', '10-20%', '>20%'),
                                            right = TRUE, include.lowest = TRUE)
category_counts <- infant_metadata %>%
  group_by(study_br, Contamination_Category) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count) * 100)

# abundance
infant_metadata$Contamination_Category_ab <- cut(infant_metadata$perc_st_iNC_ab,
                                               breaks = c(-Inf, 0, 10, 20, Inf),
                                               labels = c('0', '=<10%', '10-20%', '>20%'),
                                               right = TRUE, include.lowest = TRUE)

category_counts_ab <- infant_metadata %>%
  group_by(study_br, Contamination_Category_ab) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count) * 100)




#median(infant_metadata$shared_st_iNC_perc) #1.02
#median(infant_metadata[infant_metadata$shared_iNC_perc>0,]$perc_st_iNC_ab)

#summary(infant_metadata$perc_st_iNC_ab)
#sum(infant_metadata$perc_st_iNC_ab >=10)/1108 # 11%

# Plotting the stacked bar plot

Figure_3d <- ggplot(category_counts, aes(x = study_br, y = Proportion, fill = Contamination_Category)) +
  geom_bar(stat = "identity") +
  labs(title = "Contaminant-attributed\nrichness across studies",
       x = "Study", 
       y = "Percentage of infant samples",
       fill="Richness contamination",
       tag = "d") +
  scale_fill_manual(values = c('skyblue', 'lightgreen', 'orange', 'lightcoral')) +
  theme_minimal() +
  theme(axis.text.x = ggtext::element_markdown(size=8, angle=45, hjust=1),
        axis.text.y=element_text(size=8),
        axis.title = element_text(size=9),
        legend.position = "bottom",
        legend.text=element_text(size=8),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9, hjust=0.5, vjust=-1),
        plot.tag = element_text(face="bold")) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE, title.position = 'top', title.hjust = 0.5)) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  geom_hline(aes(yintercept=10), linetype="dashed", color = "darkblue")


Figure_3e <- ggplot(category_counts_ab, aes(x = study_br, y = Proportion, fill = Contamination_Category_ab)) +
  geom_bar(stat = "identity") +
  labs(title = "Contaminant-attributed\nabundance across studies",
       x = "Study", 
       y = "Percentage of infant samples",
       fill="Abundance of contaminants",
       tag="e") +
  scale_fill_manual(values = c('skyblue', 'lightgreen', 'orange', 'lightcoral')) +
  theme_minimal() +
  theme(axis.text.x = ggtext::element_markdown(size=8, angle=45, hjust=1),
        axis.text.y=element_text(size=8),
        axis.title = element_text(size=9),
        legend.position = "bottom",
        legend.text=element_text(size=8),
        legend.title = element_text(size=9),
        plot.title = element_text(size=9, hjust=0.5, vjust=-1),
        plot.tag = element_text(face="bold")) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE, title.position = 'top', title.hjust = 0.5)) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  geom_hline(aes(yintercept=10), linetype="dashed", color = "darkblue")


# FIGURE 3e
stat.test3f <- df_for_figure3abe %>%
  group_by(cohort) %>%
  t_test(same_cohort_NC_presence ~ Type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test3f <- stat.test3f %>% add_xy_position(x = "N_strain_shared")
stat.test3f$y.position <- c(47, 12)
stat.test3f$xmin <- 1
stat.test3f$xmax <- 2

stat.test3f$p.signif <- c("****", "****")

figure_3f <- ggplot(df_for_figure3abe, aes(x = Type, y = N_strain_shared)) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outlier.alpha = 0) +
  facet_grid(. ~ cohort, scales = "free", labeller = labeller(
    cohort = c(
      "garmaeva" = "Samples<br>Garmaeva *et al.*",
      "maqsood" = "Samples<br>Maqsood *et al.*"
    )
  )) +
  labs(y = "N strains shared with NCs", tag="f") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=8),
    axis.title = element_text(size = 9),
    axis.title.x = element_text(size=9, vjust=1),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=13)
  )+ 
  stat_pvalue_manual(stat.test3f, tip.length = 0.02, size=2.5, label = "p.signif")



# FIGURE 3f
figure_3g <- readRDS(file="micro_bacteroides.rds")
figure_3g$labels$tag <- "g"
figure_3g$labels$title <- "Bacteroides phage L6428"


# merging panels:
f3ABC <- (figure_3A | figure_3B | figure_3C) + plot_layout(widths = c(3,3,4))
f3DEF <- (Figure_3d | Figure_3e | figure_3f) + plot_layout(widths = c(3.5,3.5,3))

combined_plot3 <- f3ABC /  f3DEF / figure_3g + plot_layout(heights  = c(3.25, 3.25, 3.5))
ggsave("combined_figure3_new.pdf", combined_plot3, width = 21/2.54, height = 27/2.54)

####################################################################################################################################################

