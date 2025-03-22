library(tidyverse)
library(reshape2)

RPKM_host <- read.delim('../data/RPKM_host.tsv')
metadata_work <- read.table("../data/Sample_metadata_upd.tsv", header=TRUE, sep='\t', quote = "")

set.seed(123)
ncs <- metadata_work$Sample_name[metadata_work$Type == "Neg_ctrl" & 
                                   metadata_work$Sample_name %in% colnames(RPKM_host)]
samples_l <- sample(metadata_work$Sample_name[metadata_work$Type == "Infant" & metadata_work$Study == "liang" & metadata_work$timepoint_custom_category == "Infant (age < 5 months)" &
                                              metadata_work$Sample_name %in% colnames(RPKM_host)], 10)

samples_g <- sample(metadata_work$Sample_name[metadata_work$Type == "Infant" & metadata_work$Study == "garmaeva" &
                                                metadata_work$Sample_name %in% colnames(RPKM_host)], 1)

samples_m <- sample(metadata_work$Sample_name[metadata_work$Type == "Infant" & metadata_work$Study == "maqsood" &
                                                metadata_work$Sample_name %in% colnames(RPKM_host)], 7)

samples_s <- sample(metadata_work$Sample_name[metadata_work$Type == "Infant" & metadata_work$Study == "shah" &
                                                metadata_work$Sample_name %in% colnames(RPKM_host)], 8)

RPKM_host_filt <- RPKM_host[, colnames(RPKM_host) %in% c(ncs, samples_l, samples_g, samples_m, samples_s)]
RPKM_host_filt <- RPKM_host_filt[rowSums(RPKM_host_filt) > 0, colSums(RPKM_host_filt) > 0]

RPKM_host_filt_normalized <- as.data.frame(apply(RPKM_host_filt, 2, function(x) x / sum(x, na.rm = TRUE)))

for_filt <- as.data.frame(rowSums((RPKM_host_filt > 0) == 1))
df_row_max <- data.frame(MaxValue = apply(RPKM_host_filt_normalized, 1, max, na.rm = TRUE), 
                         row.names = rownames(RPKM_host_filt_normalized))
for_filt <- merge(for_filt, df_row_max, by="row.names")
colnames(for_filt) <- c("bug", "prevalence", "max")
bugs <- for_filt$bug[for_filt$prevalence > 2 & for_filt$max > 0.05]

RPKM_host_filt_normalized_plot <- RPKM_host_filt_normalized
RPKM_host_filt_normalized_plot[nrow(RPKM_host_filt_normalized_plot) +1, ] <- 
  colSums(RPKM_host_filt_normalized_plot[!(row.names(RPKM_host_filt_normalized_plot) %in% bugs), ])
rownames(RPKM_host_filt_normalized_plot)[nrow(RPKM_host_filt_normalized_plot)] <- "Other"
RPKM_host_filt_normalized_plot <- RPKM_host_filt_normalized_plot[row.names(RPKM_host_filt_normalized_plot)
                                                                 %in% c(bugs, "Other"), ]
RPKM_host_filt_normalized_plot$bug <- row.names(RPKM_host_filt_normalized_plot)
row.names(RPKM_host_filt_normalized_plot) <- NULL
RPKM_host_filt_normalized_plot_melt <- melt(RPKM_host_filt_normalized_plot, id.vars = "bug")
colnames(RPKM_host_filt_normalized_plot_melt) <- c("bug", "Sample_name", "ra")
RPKM_host_filt_normalized_plot_melt <- merge(RPKM_host_filt_normalized_plot_melt, 
                                             metadata_work[c("Sample_name", "status", "Study")],
                                             by = "Sample_name", all.x = T)

RPKM_host_filt_normalized_t <- as.data.frame(t(RPKM_host_filt_normalized))
df_pair <- unique(RPKM_host_filt_normalized_plot_melt[c("Sample_name", "status", "Study")])

RPKM_host_filt_normalized_t$Sample_name <- row.names(RPKM_host_filt_normalized_t)
df_pair <- merge(df_pair, RPKM_host_filt_normalized_t[c("Sample_name", "Bacteroides", "Escherichia", "Bifidobacterium", "Burkholderia")], 
                 by="Sample_name", all.x = T)

df_pair <- df_pair %>%
  arrange(status, -Bacteroides, -Escherichia, -Bifidobacterium, -Burkholderia) %>%
  group_by(status) %>% 
  mutate(pair = rep(1:26, length.out = n())) %>%  # Adjust the number of repetitions based on group size
  ungroup()

RPKM_host_filt_normalized_plot_melt <- merge(RPKM_host_filt_normalized_plot_melt,
                                             df_pair[c("Sample_name", "pair")], 
                                             by = "Sample_name", all.x = T)


RPKM_host_filt_normalized_plot_melt$bug <- as.factor(RPKM_host_filt_normalized_plot_melt$bug)
RPKM_host_filt_normalized_plot_melt$bug <- factor(RPKM_host_filt_normalized_plot_melt$bug, 
                                                  levels = c("Bacteroides", "Escherichia", "Bifidobacterium", "Burkholderia", "Alistipes", 
                                                             "Blautia_A", "Clostridium", "Enterococcus", "Faecalibacterium", "Flavonifractor",
                                                             "Haemophilus_D", "Hungatella", "Klebsiella", "Lawsonibacter", "Phocaeicola",
                                                             "Prevotella", "Pseudomonas_E", "Streptococcus", "Veillonella",
                                                             "Other", "Unclassified"))

host_aggregates_ra <- ggplot(RPKM_host_filt_normalized_plot_melt, aes(x = pair, y = ra, fill = bug)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("tomato2", "steelblue", "yellow2", "green3", "purple3", "turquoise3",
                               "gold2", "aquamarine2", "orchid",  "red", "blue",
                               "violetred1", "darkgoldenrod1", "firebrick", "palegreen3", "seagreen",
                               "orange2", "pink", "navajowhite", "slategray3", "grey35")) +
  facet_grid(status ~ ., scales="free") +
  theme(
    strip.text = ggtext::element_markdown(size=9),
    strip.background = element_rect(fill = "transparent"),
    axis.text.x=element_blank(), 
    axis.ticks.x=element_blank(),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"),
    legend.position = "bottom"
  ) +
  labs(x = "", y = "Host aggregated vOTU abundance", fill = "Host bacterial genus") +
  guides(fill = guide_legend(title.position = "top"))
ggsave("../plots/neg_ctrl_sharing/host_aggregates_ra.png", host_aggregates_ra, width = 21/2.54, height = 21/2.54)

df_pair <- merge(df_pair, metadata_work[c("Sample_name", "Type", "Timepoint")], by="Sample_name")
