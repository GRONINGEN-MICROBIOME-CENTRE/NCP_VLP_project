## Volcano plot ##
library(tidyverse)
library(ggplot2)
library(EnhancedVolcano)

result_diff_ab_l <- read.delim('../data/result_diff_ab_l_upd.tsv')
result_diff_ab_l$new_label <- gsub(".*(N\\d+_L\\d+).*", "\\1", result_diff_ab_l$vOTU)
result_diff_ab_m <- read.delim('../data/result_diff_ab_m_upd.tsv')
result_diff_ab_m$new_label <- gsub(".*(N\\d+_L\\d+).*", "\\1", result_diff_ab_m$vOTU)
result_diff_ab_s <- read.delim('../data/result_diff_ab_s_upd.tsv')
result_diff_ab_s$new_label <- gsub(".*(N\\d+_L\\d+).*", "\\1", result_diff_ab_s$vOTU)



volcano_maqsood <- EnhancedVolcano(result_diff_ab_m,
                                   title = 'Study Maqsood et al',
                                   pCutoff = 0.05,
                                   pointSize = 2.0,
                                   labSize = 2,
                                   boxedLabels = TRUE,
                                   parseLabels = TRUE,
                                   lab = result_diff_ab_m$new_label,
                                   x = 'log2fold_change',
                                   y = 'pval',
                                   subtitle = "Total 6 vOTUs",
                                   caption = "",
                                   axisLabSize = 10,
                                   titleLabSize = 12,
                                   subtitleLabSize = 10,
                                   legendLabSize = 10,
                                   legendIconSize = 5
                                   ) +
  labs(tag="a") +
  theme(plot.tag = element_text(face="bold", size=12)) +
  guides(color = "none")


volcano_liang <- EnhancedVolcano(result_diff_ab_l,
                                 title = 'Study Liang et al',
                                 pCutoff = 0.05,
                                 pointSize = 2.0,
                                 labSize = 4,
                                 lab = result_diff_ab_l$new_label,
                                 x = 'log2fold_change',
                                 y = 'pval',
                                 subtitle = "Total 58 vOTUs",
                                 caption = "",
                                 axisLabSize = 10,
                                 titleLabSize = 12,
                                 subtitleLabSize = 10,
                                 legendLabSize = 10,
                                 legendIconSize = 5) +
  labs(tag="b") +
  theme(plot.tag = element_text(face="bold", size=12))


volcano_shah <- EnhancedVolcano(result_diff_ab_s,
                                title = 'Study Shah et al',
                                pCutoff = 0.05,
                                pointSize = 2.0,
                                labSize = 4,
                                lab = result_diff_ab_s$new_label,
                                x = 'log2fold_change',
                                y = 'pval',
                                subtitle = "Total 2,914 vOTUs",
                                caption = "",
                                axisLabSize = 10,
                                titleLabSize = 12,
                                subtitleLabSize = 10,
                                legendLabSize = 10,
                                legendIconSize = 5) +
  labs(tag="c") +
  theme(plot.tag = element_text(face="bold", size=12))


figure_4_supp <- layout <- (volcano_maqsood + volcano_liang) / (volcano_shah + plot_spacer()) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

figure_4_supp


ggsave("../plots/neg_ctrl_sharing/Supplementary_figure4_diff_votus_revision.pdf", figure_4_supp, width = 20/2.54, height = 27/2.54)

# result_diff_ab_m$pval <- round(result_diff_ab_m$pval, 3)
# result_diff_ab_m$log2fold_change <- round(result_diff_ab_m$log2fold_change, 3)
# result_diff_ab_m <- result_diff_ab_m[c("vOTU", "new_label", "log2fold_change", "pval")]
# 
# result_diff_ab_l$pval <- round(result_diff_ab_l$pval, 3)
# result_diff_ab_l$log2fold_change <- round(result_diff_ab_l$log2fold_change, 3)
# result_diff_ab_l <- result_diff_ab_l[c("vOTU", "new_label", "log2fold_change", "pval")]
# 
# result_diff_ab_s$pval <- formatC(result_diff_ab_s$pval, format = "e", digits = 3)
# result_diff_ab_s$log2fold_change <- round(result_diff_ab_s$log2fold_change, 3)
# result_diff_ab_s <- result_diff_ab_s[c("vOTU", "new_label", "log2fold_change", "pval")]
# 
# write.table(result_diff_ab_m, "../result_diff_ab_m_upd_label.tsv", sep='\t', row.names=F, col.names=T, quote=F)
# write.table(result_diff_ab_l, "../result_diff_ab_l_upd_label.tsv", sep='\t', row.names=F, col.names=T, quote=F)
# write.table(result_diff_ab_s, "../result_diff_ab_s_upd_label.tsv", sep='\t', row.names=F, col.names=T, quote=F)
