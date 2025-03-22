library(tidyverse)
library(psych)

estimation <- read.table("../data/updated_NC_counts.txt", header=TRUE, sep='\t', quote = "")
metadata_work <- read.table("../data/Sample_metadata_upd.tsv", header=TRUE, sep='\t', quote = "")
df_for_figure9as <- read.table("../data/df_for_figure9as.tsv", header=TRUE, sep='\t', quote = "")

estimation$Sample_name <-  row.names(estimation)
estimation$Sample_name <- gsub('"', '', estimation$Sample_name)
estimation$Sample_name <- gsub('_dedup', '', estimation$Sample_name)
row.names(estimation) <-  NULL
working_table <- merge(metadata_work, estimation, by="Sample_name", all.x=T)

working_table$X.filter_50percCoverage. <- working_table$X.filter_50percCoverage./working_table$clean_reads_comb
working_table$X.filter_75percCoverage. <- working_table$X.filter_75percCoverage./working_table$clean_reads_comb

table_9b <- working_table[working_table$status == "SAMPLES" & !is.na("X.filter_75percCoverage.") & !is.na("presence_strain_shared_perc"), ]

spearman_result_garmaeva_9b <- psych::corr.test(table_9b[table_9b$Study == "garmaeva", c("presence_strain_shared_perc", "X.filter_75percCoverage.")], method = "spearman")
print(spearman_result_garmaeva_9b$r)
print(spearman_result_garmaeva_9b$p)

spearman_result_maqsood_9b <- psych::corr.test(table_9b[table_9b$Study == "maqsood", c("presence_strain_shared_perc", "X.filter_75percCoverage.")], method = "spearman")
print(spearman_result_maqsood_9b$r)
print(spearman_result_maqsood_9b$p)

spearman_result_liang_9b <- psych::corr.test(table_9b[table_9b$Study == "liang", c("presence_strain_shared_perc", "X.filter_75percCoverage.")], method = "spearman")
print(spearman_result_liang_9b$r)
print(spearman_result_liang_9b$p)

spearman_result_shah_9b <- psych::corr.test(table_9b[table_9b$Study == "shah", c("presence_strain_shared_perc", "X.filter_75percCoverage.")], method = "spearman")
print(spearman_result_shah_9b$r)
print(spearman_result_shah_9b$p)


correlations_presence_strain_vs_diffvotu <- working_table[working_table$status == "SAMPLES" & !is.na(working_table$X.filter_75percCoverage.) & !is.na(working_table$presence_strain_shared_perc), ] %>%
  group_by(Study) %>%
  summarize(cor = cor(presence_strain_shared_perc, X.filter_75percCoverage., method = "spearman"))


figure_9A_supp <- ggplot(df_for_figure9as, aes(x = same_diff, y = value)) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outliers = FALSE) +
  facet_grid(cohort ~ pres_abun, scales = "free", labeller = labeller(
    pres_abun = c(
      "presence" = "% vOTUs shared (presence)",
      "abundance" = "Abundance of shared vOTUs "
    ),
    cohort = c(
      "garmaeva" = "Samples Garmaeva *et al.* <br>",
      "liang" = "Samples Liang *et al.* <br>",
      "maqsood" = "Samples Maqsood *et al.* <br>",
      "shah" = "Samples Shah *et al.* <br>"
    )
  )) +
  labs(x="Study", y="% vOTUs shared with NCs", tag="a") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=7),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 8), 
    axis.text.y = element_text(size = 7),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=10)
  )

stat.test9as <- df_for_figure9as %>%
  group_by(cohort, pres_abun) %>%
  t_test(value ~ same_diff) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test9as <- stat.test9as %>% 
  add_xy_position(x = "same_diff", dodge = 0.8)

stat.test9as$xmin <- 1
stat.test9as$xmax <- 2

stat.test9as$p.signif <- c("****", "****", "ns", "ns", "ns", "ns", "****", "****")

figure_9A_supp <- figure_9A_supp + 
  stat_pvalue_manual(stat.test9as, tip.length = 0.02, size=2.5, label = "p.signif")


figure_9B_supp <- ggplot(working_table[working_table$status == "SAMPLES", ], aes(x=presence_strain_shared_perc, y=X.filter_75percCoverage.)) +
  geom_point(size = 1.3, color="#2E236C", alpha=0.75) +
  geom_smooth(method="lm", color="#2E236C", fill="#C8ACD6", se=TRUE, linewidth=0.5) +
  facet_wrap(~ Study, nrow = 2, ncol = 2, scales = "free", labeller = labeller(
    Study = c(
      "garmaeva" = "Samples Garmaeva *et al.* <br>",
      "liang" = "Samples Liang *et al.* <br>",
      "maqsood" = "Samples Maqsood *et al.* <br>",
      "shah" = "Samples Shah *et al.* <br>"
    )
  )) +
  labs(x = "% shared strains with NCs from own studies", y = "Proportion of reads mapped to the vOTU \nrepresentatives from external NCs", tag="b") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=7),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 8), 
    axis.text.y = element_text(size = 7),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=10)
  )+
  geom_text(data = correlations_presence_strain_vs_diffvotu, aes(x = Inf, y = Inf, label = paste("rho = ", round(cor, 2))),
            hjust = 1.1, vjust = 1.5, size = 3)

figure_9_supp <- figure_9A_supp | figure_9B_supp + 
  plot_layout( 
    heights = c(1, 0.5))

ggsave("../plots/neg_ctrl_sharing/Supplementary_figure9_revision.pdf", figure_9_supp, width = 21/2.54, height = 21/2.54)

