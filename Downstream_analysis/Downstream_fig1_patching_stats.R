## Code description
#######################################################################################################################################
## Script for the negativeome paper
#  If you would like to reproduce the results outlined in the paper, please note that column names and filenames have been modified in
#  the deposited versions of the files in the database to improve readability. Additionally, some columns were removed from the initial
#  tables, as indicated at the end of the script, since they were not used in the analysis. The code for renaming columns and saving
#  "Sample_metadata.tsv" file is located at the end of the script. The code for "vOTU_repesentatives_metadata.tsv" creation can be
#  found here: https://github.com/GRONINGEN-MICROBIOME-CENTRE/NCP_VLP_project/blob/master/Upstream_analysis/NCP_all/dereplication_stat.R
#######################################################################################################################################

## Load libraries
#######################################################################################################################################
.libPaths("/home1/p309176/R_libraries")
library(readr)
library(vegan)
library(dplyr)
library(stringr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggtext)
library(VennDiagram)
library(lmerTest)
library(patchwork)
library(MuMIn)
library(psych)
library(ggpubr)
library(rstatix)
library(MetBrewer)
library(ggExtra)
#######################################################################################################################################

## Set the working directory
#######################################################################################################################################
setwd("/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/plots")
#######################################################################################################################################

## Load the  data
#######################################################################################################################################
extended_tof <- read.delim('../../VIR_DB/virus_contigs/MERGED_Extended_TOF_NCP')
row.names(extended_tof) <- extended_tof$New_CID
dim(extended_tof)
extended_tof <- extended_tof %>%
  mutate(virus_group = ifelse(grepl("Duplodnaviria", taxonomy), "dsDNA", "Unclassified"),
         virus_group = ifelse(grepl("Inoviridae", taxonomy), "ssDNA", virus_group),
         virus_group = ifelse(grepl("Microviridae", taxonomy), "ssDNA", virus_group),
         virus_group = ifelse(grepl("Cressdnaviricota", taxonomy), "ssDNA", virus_group),
         virus_group = ifelse(grepl("Cossaviricota", taxonomy), "dsDNA", virus_group),
         virus_group = ifelse(grepl("Riboviria", taxonomy), "RNA", virus_group),
         virus_group = ifelse(grepl("Nodaviridae", taxonomy), "RNA", virus_group),
         virus_group = ifelse(grepl("Tolivirales", taxonomy), "RNA", virus_group),
         virus_group = ifelse(grepl("Retroviridae", taxonomy), "RNA", virus_group),
         virus_group = ifelse(grepl("Cystoviridae", taxonomy), "RNA", virus_group),
         virus_group = ifelse(grepl("Picornaviridae", taxonomy), "RNA", virus_group),
         virus_group = ifelse(grepl("Astroviridae", taxonomy), "RNA", virus_group),
         virus_group = ifelse(grepl("Caliciviridae", taxonomy), "RNA", virus_group),
         virus_group = ifelse(grepl("Varidnaviria", taxonomy), "dsDNA", virus_group),
         virus_group = ifelse(grepl("Anelloviridae", taxonomy), "ssDNA", virus_group),
         virus_group = ifelse(grepl("Portogloboviridae", taxonomy), "dsDNA", virus_group),
         virus_group = ifelse(grepl("Bicaudaviridae", taxonomy), "dsDNA", virus_group),
         host_group = ifelse(grepl("Caudoviricetes", taxonomy), "Prokaryotes", "Unclassified"),
         host_group = ifelse(grepl("Herviviricetes", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Inoviridae", taxonomy), "Prokaryotes", host_group),
         host_group = ifelse(grepl("Microviridae", taxonomy), "Prokaryotes", host_group),
         host_group = ifelse(grepl("Genomoviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Circoviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Nanoviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Geminiviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Smacoviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Cossaviricota", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Nodaviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Tombusviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Retroviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Cystoviridae", taxonomy), "Prokaryotes", host_group),
         host_group = ifelse(grepl("Picornaviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Astroviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Caliciviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Phycodnaviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Mimiviridae", taxonomy), "Eukaryotes(giant_virus)", host_group),
         host_group = ifelse(grepl("Adenoviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Corticoviridae", taxonomy), "Prokaryotes", host_group),
         host_group = ifelse(grepl("Iridoviridae", taxonomy), "Prokaryotes", host_group),
         host_group = ifelse(grepl("Lavidaviridae", taxonomy), "viral phages", host_group),
         host_group = ifelse(grepl("Adintoviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Tectiviridae", taxonomy), "Prokaryotes", host_group),
         host_group = ifelse(grepl("Sphaerolipoviridae", taxonomy), "Prokaryotes", host_group),
         host_group = ifelse(grepl("Autolykiviridae", taxonomy), "Prokaryotes", host_group),
         host_group = ifelse(grepl("Marseilleviridae", taxonomy), "Eukaryotes(giant_virus)", host_group),
         host_group = ifelse(grepl("Poxviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Anelloviridae", taxonomy), "Eukaryotes", host_group),
         host_group = ifelse(grepl("Portogloboviridae", taxonomy), "Prokaryotes", host_group),
         host_group = ifelse(grepl("Bicaudaviridae", taxonomy), "Prokaryotes", host_group),
         checkv_quality_and_plasmid = ifelse(checkv_quality %in% c("Complete", "High-quality", "Medium-quality"), "CHM", NA),
         checkv_quality_and_plasmid = ifelse(checkv_quality %in% c("Not-determined", "Low-quality"), "LU", checkv_quality_and_plasmid),
         checkv_quality_and_plasmid = ifelse(plasmid == "Yes", "plasmid", checkv_quality_and_plasmid)
  )

RPKM_contigs_keep <- extended_tof$New_CID[extended_tof$POST_CHV_length >= 1000 & (extended_tof$viral_genes >= extended_tof$host_genes) & extended_tof$plasmid=="No"]

RPKM_initial <- read.delim("../RPKM_counts_VLP_NCP.txt")
dim(RPKM_initial)

RPKM <- RPKM_initial[row.names(RPKM_initial) %in% RPKM_contigs_keep, ]
dim(RPKM)
RPKM <- RPKM[rowSums(RPKM) > 0, colSums(RPKM) > 0]
dim(RPKM)

RPKM_count <- RPKM
RPKM_count[RPKM_count > 0] <- 1

inverse_rank_transform <- function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
RPKM_INT <- as.data.frame(apply(RPKM, 1, inverse_rank_transform))  # do inverse rank transformation (PER ROW)

meta_all_with_qc_curated <- as.data.frame(read_tsv('../../metadata_with_qc_NCPv2.tsv'))
dim(meta_all_with_qc_curated)
meta_all_with_qc_curated <- meta_all_with_qc_curated %>%
  mutate(Subject_ID = ifelse(grepl("kid", Sample_name), Sample_name, Subject_ID),
         nc_subject_group = ifelse(grepl("bctrl1", Sample_name), "NC_maqsood_buffer", "SAMPLE"),
         nc_subject_group = ifelse(grepl("bctrl4", Sample_name), "NC_maqsood_orsay", nc_subject_group),
         nc_subject_group = ifelse(grepl("ctl", Sample_name), "NC_shah_buffer", nc_subject_group),
         nc_subject_group = ifelse(grepl("LN_7C08_VL_405", Sample_name), "NC_garmaeva_buffer", nc_subject_group),
         nc_subject_group = ifelse(Sample_ID %in% c("LT1_D", "LT4_D", "LT5_D", "LT2_D", "LT3_D"), "NC_liang_tube", nc_subject_group),
         nc_subject_group = ifelse(Sample_ID %in% c("LT1_R", "LT2_R", "LT3_R", "LT4_R", "LT5_R"), "NC_liang_tube", nc_subject_group),
         nc_subject_group = ifelse(Sample_ID %in% c("LB1_D", "LB2_D", "LB3_D", "LB4_D", "LB5_D", "LNCB3_D", "LNCB2_D", "LNCB1_D"), "NC_liang_buffer", nc_subject_group),
         nc_subject_group = ifelse(Sample_ID %in% c("LB1_R", "LB2_R", "LB3_R", "LB4_R", "LB5_R", "LNCB1_R", "LNCB3_R", "LNCB2_R"), "NC_liang_buffer", nc_subject_group),
         nc_subject_group = ifelse(Sample_ID %in% c("LD1_D", "LD2_D", "LDH_D", "LDN1_D", "LDN2_D"), "NC_liang_diaper", nc_subject_group),
         nc_subject_group = ifelse(Sample_ID %in% c("LD1_R", "LD2_R", "LDH_R", "LDN1_R", "LDN2_R"), "NC_liang_diaper", nc_subject_group),
         nc_subject_group = ifelse(Sample_ID %in% c("LMDNC_D", "LMDNC_R"), "NC_liang_MDNC", nc_subject_group),
         rna_dna = ifelse(grepl("_D", Sample_ID), "DNA", NA),
         rna_dna = ifelse(grepl("_R", Sample_ID), "RNA", rna_dna),
         ncvssample = ifelse(Type == "Neg_ctrl", "NCs", "SAMPLES"),  # Changed for SAMPLES and NCs here; original - SAMPLE & NC
         nc_subject_group = ifelse(ncvssample == "SAMPLES", Subject_ID, nc_subject_group),
         timepoint_type = ifelse(Type == "Neg_ctrl", "NC", NA),
         timepoint_type = ifelse(Type == "Infant" & Timepoint %in% c("M0", "M1", "M2", "M3", "M4"), "Infant (age < 5 months)", timepoint_type),
         timepoint_type = ifelse(Type == "Infant" & Timepoint %in% c("M6", "M12", "Y2-5"), "Infant (age > 5 months)", timepoint_type),
         timepoint_type = ifelse(Type == "Mother", "Mother", timepoint_type),
         Timepoint_numeric = as.integer(gsub("M", "", Timepoint)),
         Timepoint_numeric = ifelse(grepl("Y2-5", Timepoint), 42, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother", 384, Timepoint_numeric)
  )

meta_all_with_qc_curated$nc_subject_group <- as.factor(meta_all_with_qc_curated$nc_subject_group)
meta_all_with_qc_curated$cohort <- as.factor(meta_all_with_qc_curated$cohort)
meta_all_with_qc_curated$ncvssample <- as.factor(meta_all_with_qc_curated$ncvssample)
meta_all_with_qc_curated$ncvssample <- factor(meta_all_with_qc_curated$ncvssample, levels = c("SAMPLES", "NCs"))

meta_working <- meta_all_with_qc_curated[meta_all_with_qc_curated$Sample_name %in% colnames(RPKM), ]
meta_working$ncvssample <- factor(meta_working$ncvssample, levels = c("SAMPLES", "NCs"))

host_prediction <- read.csv('../../VIR_DB/host_prediction_w_neg_der95_NCP/results/MERGED_Host_prediction_to_genus_m90_v2.csv')
dim(host_prediction)

filtered_host_prediction <- host_prediction %>%
  group_by(Virus) %>%
  slice_max(order_by = Confidence.score, with_ties = FALSE) %>%
  ungroup()
filtered_host_prediction <- as.data.frame(filtered_host_prediction)
dim(filtered_host_prediction)

strains_df_ini <- read.delim('/scratch/hb-llnext/VLP_public_data/nc_project/instrain/instrain_compare_ALL/instrain_compare_ALL_comparisonsTable.tsv')
#######################################################################################################################################

## Adding alpha diversity and richness to meta_working
#######################################################################################################################################
row.names(meta_working) <- meta_working$Sample_name
diversity_df <- as.data.frame(diversity(as.data.frame(t(RPKM)), index='shannon'))
meta_working <- merge(meta_working, diversity_df, by="row.names")
colnames(meta_working)[colnames(meta_working) == "diversity(as.data.frame(t(RPKM)), index = \"shannon\")"] <- "diversity"
meta_working$Row.names <- NULL
richness_table <- as.data.frame(colSums(RPKM_count))
row.names(meta_working) <- meta_working$Sample_name
meta_working <- merge(meta_working, richness_table, by="row.names")
colnames(meta_working)[colnames(meta_working) == "colSums(RPKM_count)"] <- "richness"
meta_working$Row.names <- NULL

meta_all_with_qc_curated <- merge(meta_all_with_qc_curated, meta_working[c("Sample_name", "richness", "diversity")], all.x = T)
meta_all_with_qc_curated$richness[is.na(meta_all_with_qc_curated$richness)] <- 0
meta_all_with_qc_curated$diversity[is.na(meta_all_with_qc_curated$diversity)] <- 0
#######################################################################################################################################

## Analysis for the vOTUs composition in NC vs samples: richness and CHM ratio detected
#######################################################################################################################################
RPKM_count_ch1 <- merge(RPKM_count, extended_tof[colnames(extended_tof) %in% c("checkv_quality_and_plasmid")], by="row.names")
RPKM_count_ch1$Row.names <- NULL

summarized_df <- RPKM_count_ch1 %>%
  group_by(checkv_quality_and_plasmid) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

summarized_df <- as.data.frame(t(summarized_df))
colnames(summarized_df) <- c("CHM", "LU")
summarized_df <- summarized_df[row.names(summarized_df) != "checkv_quality_and_plasmid", ]

summarized_df <- summarized_df %>%
  mutate(across(where(is.character), ~ as.numeric(.)))

summarized_df$CHM_LU_ratio <- summarized_df$CHM / (summarized_df$LU + summarized_df$CHM)

row.names(meta_working) <- meta_working$Sample_name
meta_working <- merge(meta_working, summarized_df, by="row.names")
meta_working$Row.names <- NULL
#######################################################################################################################################

## Analysis for the vOTUs composition in NC vs samples: richness and CHM ratio discovered
#######################################################################################################################################
get_study_and_sample_of_origin <- function(cid) {
  parts <- strsplit(cid, "_")[[1]]
  num_parts <- length(parts)
  sample_of_origin <- paste(parts[1:(num_parts - 6)], collapse = "_")
  split_parts <- strsplit(sample_of_origin, "_", fixed = TRUE)[[1]]
  study_of_origin <- split_parts[1]
  sample_of_origin <- paste(split_parts[-1], collapse = "_")
  
  return(c(study_of_origin, sample_of_origin))
}

split_results <- t(sapply(extended_tof$New_CID, get_study_and_sample_of_origin))

extended_tof$study_of_origin <- split_results[, 1]
extended_tof$sample_of_origin <- split_results[, 2]

summarized_discovered_vOTUs <- extended_tof %>%
  filter(POST_CHV_length >= 1000 & plasmid == "No") %>%
  group_by(sample_of_origin, study_of_origin, checkv_quality_and_plasmid) %>%
  summarise(count = n())

summarized_discovered_vOTUs <- summarized_discovered_vOTUs %>%
  pivot_wider(names_from = checkv_quality_and_plasmid, values_from = count, values_fill = list(count = 0))

summarized_discovered_vOTUs <- as.data.frame(summarized_discovered_vOTUs)
colnames(summarized_discovered_vOTUs) <- c("Sample_name", "cohort", "CHM_discovered", "LU_discovered")

summarized_discovered_vOTUs$CHM_LU_richness_discovered_ratio <- summarized_discovered_vOTUs$CHM_discovered / (summarized_discovered_vOTUs$LU_discovered + summarized_discovered_vOTUs$CHM_discovered)
row.names(summarized_discovered_vOTUs) <- summarized_discovered_vOTUs$Sample_name

row.names(meta_working) <- meta_working$Sample_name
meta_working <- merge(meta_working, summarized_discovered_vOTUs[, !colnames(summarized_discovered_vOTUs) %in% c("cohort", "Sample_name")], all.x = T, by="row.names")
meta_working$Row.names <- NULL
#######################################################################################################################################

## Analysis for the virus groups richness + STATS
#######################################################################################################################################
CHM_viruses <- row.names(extended_tof[extended_tof$checkv_quality_and_plasmid == "CHM", ])
RPKM_count_ch2 <- RPKM_count[row.names(RPKM_count) %in% CHM_viruses, ]
RPKM_count_ch2 <- merge(RPKM_count_ch2, extended_tof[colnames(extended_tof) %in% c("virus_group")], by="row.names")
RPKM_count_ch2$Row.names <- NULL

summarized_df2 <- RPKM_count_ch2 %>%
  group_by(virus_group) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

summarized_df2 <- as.data.frame(t(summarized_df2))
colnames(summarized_df2) <- c("RNA", "Unclassified", "dsDNA", "ssDNA")
summarized_df2 <- summarized_df2[row.names(summarized_df2) != "virus_group", ]

summarized_df2 <- summarized_df2 %>%
  mutate(across(where(is.character), ~ as.numeric(.)))

summarized_df2 <- summarized_df2[rowSums(summarized_df2) > 0, ]

fraction_per_sample <- function(row) {
  row / sum(row)
}

summarized_df2_frac <- as.data.frame(t(apply(summarized_df2, 1, fraction_per_sample)))

summarized_df2_frac$Sample_name <- row.names(summarized_df2_frac)
summarized_df2_frac <- merge(summarized_df2_frac, meta_working[, colnames(meta_working) %in% c("Sample_name", "ncvssample", "cohort", "nc_subject_group")], all.x = T, by="Sample_name")

summarized_df2_frac_stats <- melt(summarized_df2_frac, id.vars = c("Sample_name", "ncvssample", "cohort", "nc_subject_group"))

summarized_df2_frac_stats$value_log_trnsfrmd <- summarized_df2_frac_stats$value + (min(summarized_df2_frac_stats$value[summarized_df2_frac_stats$value > 0])/2)
summarized_df2_frac_stats$value_log_trnsfrmd <- log(summarized_df2_frac_stats$value_log_trnsfrmd)
summarized_df2_frac_stats$ncvssample <- factor(summarized_df2_frac_stats$ncvssample, levels = c("SAMPLES", "NCs"))

summary(lmer(value_log_trnsfrmd ~ ncvssample + (1|cohort/nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$cohort %in% c("liang", "shah") & summarized_df2_frac_stats$variable == "RNA", ]))
summary(lmer(value_log_trnsfrmd ~ ncvssample + (1|cohort/nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$variable == "dsDNA", ]))
summary(lmer(value_log_trnsfrmd ~ ncvssample + (1|cohort/nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$variable == "ssDNA", ]))

summary(lmer(value_log_trnsfrmd ~ ncvssample + (1|cohort/nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$cohort %in% c("maqsood", "shah", "garmaeva") & summarized_df2_frac_stats$variable == "ssDNA", ]))


summary(lmer(value_log_trnsfrmd ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$cohort == "liang" & summarized_df2_frac_stats$variable == "RNA", ]))
summary(lmer(value_log_trnsfrmd ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$cohort == "liang" & summarized_df2_frac_stats$variable == "dsDNA", ]))
summary(lmer(value_log_trnsfrmd ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$cohort == "liang" & summarized_df2_frac_stats$variable == "ssDNA", ]))

summary(lmer(value_log_trnsfrmd ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$cohort == "maqsood" & summarized_df2_frac_stats$variable == "dsDNA", ]))
summary(lmer(value_log_trnsfrmd ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$cohort == "maqsood" & summarized_df2_frac_stats$variable == "ssDNA", ]))

summary(lmer(value_log_trnsfrmd ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$cohort == "shah" & summarized_df2_frac_stats$variable == "RNA", ]))
summary(lmer(value_log_trnsfrmd ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$cohort == "shah" & summarized_df2_frac_stats$variable == "dsDNA", ]))
summary(lmer(value_log_trnsfrmd ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df2_frac_stats[summarized_df2_frac_stats$cohort == "shah" & summarized_df2_frac_stats$variable == "ssDNA", ]))

p.adjust(c(0.39, 0.0696, 0.0046, 0.4735, 0.569, 2.85e-06, 0.41734, 0.8717), method="BH")

summarized_df2_frac_stats$variable <- as.factor(summarized_df2_frac_stats$variable)
summarized_df2_frac_stats$variable <- factor(summarized_df2_frac_stats$variable, levels = c("RNA", "ssDNA", "dsDNA", "Unclassified"))
#######################################################################################################################################

## Analysis for the virus host richness -> prokaryotes vs eukaryotes + stats
#######################################################################################################################################
RPKM_count_ch3 <- RPKM_count[row.names(RPKM_count) %in% CHM_viruses, ]
RPKM_count_ch3 <- merge(RPKM_count_ch3, extended_tof[colnames(extended_tof) %in% c("host_group")], by="row.names")
RPKM_count_ch3$Row.names <- NULL

summarized_df3 <- RPKM_count_ch3 %>%
  group_by(host_group) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

summarized_df3 <- as.data.frame(t(summarized_df3))
colnames(summarized_df3) <- c("Eukaryotes", "Prokaryotes", "Unclassified")
summarized_df3 <- summarized_df3[row.names(summarized_df3) != "host_group", ]

summarized_df3 <- summarized_df3 %>%
  mutate(across(where(is.character), ~ as.numeric(.)))

summarized_df3 <- summarized_df3[rowSums(summarized_df3) > 0, ]

summarized_df3_frac <- as.data.frame(t(apply(summarized_df3, 1, fraction_per_sample)))

summarized_df3_frac$Sample_name <- row.names(summarized_df3_frac)
summarized_df3_frac <- merge(summarized_df3_frac, meta_working[, colnames(meta_working) %in% c("Sample_name", "ncvssample", "cohort", "nc_subject_group")], all.x = T, by="Sample_name")
summarized_df3_frac$Sample_name <- NULL

summarized_df3_frac_stats <- melt(summarized_df3_frac, id.vars = c("ncvssample", "cohort", "nc_subject_group"))

summarized_df3_frac_stats$ncvssample <- factor(summarized_df3_frac_stats$ncvssample, levels = c("SAMPLES", "NCs"))

summarized_df3_frac_stats$value_log_trnsfrmd <- summarized_df3_frac_stats$value + (min(summarized_df3_frac_stats$value[summarized_df3_frac_stats$value > 0])/2)
summarized_df3_frac_stats$value_log_trnsfrmd <- log(summarized_df3_frac_stats$value_log_trnsfrmd)

summary(lmer(value_log_trnsfrmd ~ ncvssample + (1|cohort/nc_subject_group), REML = F, data = summarized_df3_frac_stats[summarized_df3_frac_stats$variable == "Eukaryotes", ]))
summary(lmer(value_log_trnsfrmd ~ ncvssample + (1|cohort/nc_subject_group), REML = F, data = summarized_df3_frac_stats[summarized_df3_frac_stats$variable == "Prokaryotes", ]))
p.adjust(c(0.061, 0.177), method="BH")

summary(lmer(value_log_trnsfrmd ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df3_frac_stats[summarized_df3_frac_stats$cohort == "liang" & summarized_df3_frac_stats$variable == "Eukaryotes", ]))
summary(lmer(value_log_trnsfrmd ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df3_frac_stats[summarized_df3_frac_stats$cohort == "liang" & summarized_df3_frac_stats$variable == "Prokaryotes", ]))

summary(lmer(value_log_trnsfrmd ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df3_frac_stats[summarized_df3_frac_stats$cohort == "maqsood" & summarized_df3_frac_stats$variable == "Eukaryotes", ]))
summary(lmer(value_log_trnsfrmd ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df3_frac_stats[summarized_df3_frac_stats$cohort == "maqsood" & summarized_df3_frac_stats$variable == "Prokaryotes", ]))

summary(lmer(value_log_trnsfrmd ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df3_frac_stats[summarized_df3_frac_stats$cohort == "shah" & summarized_df3_frac_stats$variable == "Eukaryotes", ]))
summary(lmer(value_log_trnsfrmd ~ ncvssample + (1|nc_subject_group), REML = F, data = summarized_df3_frac_stats[summarized_df3_frac_stats$cohort == "shah" & summarized_df3_frac_stats$variable == "Prokaryotes", ]))
p.adjust(c(0.093, 0.413, 0.229289, 0.443, 0.507, 3.04e-06), method="BH")

# Answering if prokaryotes dominating eukaryotes

summarized_df3_frac_stats2 <- summarized_df3_frac

summarized_df3_frac_stats2$Prokaryotes_log_trnsfrmd <- summarized_df3_frac_stats2$Prokaryotes + (min(summarized_df3_frac_stats2$Prokaryotes[summarized_df3_frac_stats2$Prokaryotes > 0])/2)
summarized_df3_frac_stats2$Eukaryotes_log_trnsfrmd <- summarized_df3_frac_stats2$Eukaryotes + (min(summarized_df3_frac_stats2$Eukaryotes[summarized_df3_frac_stats2$Eukaryotes > 0])/2)

summarized_df3_frac_stats2$Prokaryotes_log_trnsfrmd <- log(summarized_df3_frac_stats2$Prokaryotes_log_trnsfrmd)
summarized_df3_frac_stats2$Eukaryotes_log_trnsfrmd <- log(summarized_df3_frac_stats2$Eukaryotes_log_trnsfrmd)

prok_vs_euk_combined <- summary(lmer(Prokaryotes_log_trnsfrmd ~ Eukaryotes_log_trnsfrmd + (1|cohort/nc_subject_group), REML = F, data = summarized_df3_frac_stats2))
prok_vs_euk_combined$coefficients
prok_vs_euk_liang <- summary(lmer(Prokaryotes_log_trnsfrmd ~ Eukaryotes_log_trnsfrmd + (1|nc_subject_group), REML = F, data = summarized_df3_frac_stats2[summarized_df3_frac_stats2$cohort == "liang", ]))
prok_vs_euk_liang$coefficients
prok_vs_euk_maqsood <- summary(lmer(Prokaryotes_log_trnsfrmd ~ Eukaryotes_log_trnsfrmd + (1|nc_subject_group), REML = F, data = summarized_df3_frac_stats2[summarized_df3_frac_stats2$cohort == "maqsood", ]))
prok_vs_euk_maqsood$coefficients
prok_vs_euk_shah <- summary(lmer(Prokaryotes_log_trnsfrmd ~ Eukaryotes_log_trnsfrmd + (1|nc_subject_group), REML = F, data = summarized_df3_frac_stats2[summarized_df3_frac_stats2$cohort == "shah", ]))
prok_vs_euk_shah$coefficients
p.adjust(c(5.363049e-31, 1.41e-06, 7.386918e-112))

## Redo for long format
summarized_df3_frac_stats$variable <- as.factor(summarized_df3_frac_stats$variable)
summarized_df3_frac_stats$variable <- factor(summarized_df3_frac_stats$variable, levels = c("Prokaryotes", "Eukaryotes", "Unclassified"))
#######################################################################################################################################

## Analysis for the phages host richness
#######################################################################################################################################

BU_viruses <- row.names(extended_tof[extended_tof$host_group %in% c("Prokaryotes", "Unclassified"), ])
RPKM_ch4 <- RPKM[row.names(RPKM) %in% CHM_viruses & row.names(RPKM) %in% BU_viruses, ]

row.names(filtered_host_prediction) <- filtered_host_prediction$Virus
filtered_host_prediction$Genus <- sub(".*;g__", "", filtered_host_prediction$Host.genus)
filtered_host_prediction$Genus[filtered_host_prediction$Genus == ""] <- "Unclassified"

RPKM_ch4 <- merge(RPKM_ch4, filtered_host_prediction[colnames(filtered_host_prediction) %in% c("Genus")], all.x = T, by="row.names")
RPKM_ch4$Genus[is.na(RPKM_ch4$Genus)] <- "Unclassified"
RPKM_ch4$Row.names <- NULL

summarized_df4 <- RPKM_ch4 %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

summarized_df4 <- as.data.frame(t(summarized_df4))
colnames(summarized_df4) <- c(t(summarized_df4[1,]))
summarized_df4 <- summarized_df4[row.names(summarized_df4) != "Genus", ]
summarized_df4 <- summarized_df4 %>%
  mutate(across(where(is.character), ~ as.numeric(.)))
RPKM_host_t <- summarized_df4
RPKM_host_t <- RPKM_host_t[rowSums(RPKM_host_t) > 0, colSums(RPKM_host_t) > 0]

# RPKM_host <- as.data.frame(t(RPKM_host_t))
# write.table(RPKM_host, "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/RPKM_host.tsv", sep='\t', row.names=T, col.names=T, quote=F)
#######################################################################################################################################

## Additional analysis: NMDS
#######################################################################################################################################

meta_nmds_all <- meta_working
row.names(meta_nmds_all) <- meta_nmds_all$Sample_name
meta_nmds_all <- meta_nmds_all[ colnames(RPKM), ]  # SORTING FOR ENV FIT
Isfahan2 <- met.brewer('Isfahan2')  # color palette

# This commented block generates the data for NMDS, and needs to be run only once
# ord_all_vir <- metaMDS(t(RPKM), distance = "bray", k=2)
# data.scores.all.vir = as.data.frame(scores(ord_all_vir, "sites"))
# write.table(data.scores.all.vir, "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/datascores_all_vir.tsv", sep='\t', row.names=T, col.names=T, quote=F)
# 
# ord_all_host <- metaMDS(RPKM_host_t, distance = "bray", k=2)  # check previous block to retrieve RPKM
# data.scores.all.hosts = as.data.frame(scores(ord_all_host, "sites"))
# write.table(data.scores.all.hosts, "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/datascores_all_hosts.tsv", sep='\t', row.names=T, col.names=T, quote=F)

data.scores.all.vir <- read.table('../datascores_all_vir.tsv', sep='\t', header=T)
data.scores.all.hosts <- read.table('../datascores_all_hosts.tsv', sep='\t', header=T)

# MERGING NMDS DATA & METADATA 
data.scores.all.vir <- merge(data.scores.all.vir, meta_nmds_all, by='row.names', all.x=T)
row.names(data.scores.all.vir) <- data.scores.all.vir$Row.names
data.scores.all.vir$Row.names <- NULL

data.scores.all.hosts <- merge(data.scores.all.hosts, meta_nmds_all, by='row.names', all.x=T)
row.names(data.scores.all.hosts) <- data.scores.all.hosts$Row.names
data.scores.all.hosts$Row.names <- NULL

data.scores.all.vir$timepoint_type <- factor(data.scores.all.vir$timepoint_type, levels=c("Infant (age < 5 months)", "Infant (age > 5 months)", "Mother", "NC"), ordered = T)
data.scores.all.hosts$timepoint_type <- factor(data.scores.all.hosts$timepoint_type, levels=c("Infant (age < 5 months)", "Infant (age > 5 months)", "Mother", "NC"), ordered = T)

data.scores.all.vir <- data.scores.all.vir[data.scores.all.vir$NMDS1 >= 0.0041 & data.scores.all.vir$NMDS1 <= 0.0059
                                                & data.scores.all.vir$NMDS2 >= 0.004 & data.scores.all.vir$NMDS2 <= 0.006, ]

data.scores.all.hosts <- data.scores.all.hosts[data.scores.all.hosts$NMDS1 >= -0.003 & data.scores.all.hosts$NMDS1 <= 0.009
                                               & data.scores.all.hosts$NMDS2 >= -0.0075 & data.scores.all.hosts$NMDS2 <= 0.009, ]

# Extracting centroids
en.vir = envfit(data.scores.all.vir[c("NMDS1", "NMDS2")], data.scores.all.vir[c("timepoint_type")], permutations = 999, na.rm = TRUE)
en.vir$factors
en.vir$vectors

en.hosts = envfit(data.scores.all.hosts[c("NMDS1", "NMDS2")], data.scores.all.hosts[c("timepoint_type")], permutations = 999, na.rm = TRUE)
en.hosts$factors
en.hosts$vectors

centroids.vir <- as.data.frame(scores(en.vir, "factors"))
centroids.vir$timepoint_type <- c(gsub('timepoint_type', '', row.names(centroids.vir)))

centroids.host <- as.data.frame(scores(en.hosts, "factors"))
centroids.host$timepoint_type <- c(gsub('timepoint_type', '', row.names(centroids.host)))

# Plotting
NMDS_viruses_plot <- ggplot(data = data.scores.all.vir, aes(x = NMDS1, y = NMDS2, color=timepoint_type)) + 
  geom_point(size = 1.5, alpha=0.6) + 
  geom_point(data=centroids.vir, aes(fill=timepoint_type),shape=23, size=4, color='black', ) + 
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = timepoint_type, color=timepoint_type), linetype = 2) +
  xlim(0.0041,0.0059) +
  ylim(0.004, 0.006) +
  theme_bw()+
  labs(color = "Timepoints", title = "Bray-Curtis dissimilarity based on vOTUs") +
  scale_color_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  theme(
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.position = "bottom") +
  guides(fill = "none")
NMDS_viruses_plot <- ggMarginal(NMDS_viruses_plot, type="boxplot", groupFill=T)
ggsave("NMDS_viruses.png", NMDS_viruses_plot, width = 14/2.54, height = 14/2.54)

NMDS_hosts_plot <- ggplot(data = data.scores.all.hosts, aes(x = NMDS1, y = NMDS2, color=timepoint_type)) + 
  geom_point(size = 1.5, alpha=0.6) + 
  geom_point(data=centroids.host, aes(fill=timepoint_type),shape=23, size=4, color='black', ) + 
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = timepoint_type, color=timepoint_type), linetype = 2) +
  xlim(-0.003,0.009) +
  ylim(-0.0075,0.009) +
  theme_bw()+
  labs(color = "Timepoints", title = "Bray-Curtis dissimilarity based on host-based vOTU aggregates") +
  scale_color_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  theme(
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    strip.background = element_rect(fill = "transparent"),
    legend.position = "bottom") +
  guides(fill = "none") 
NMDS_hosts_plot <- ggMarginal(NMDS_hosts_plot, type="boxplot", groupFill=T)
ggsave("NMDS_hosts.png", NMDS_hosts_plot, width = 14/2.54, height = 14/2.54)

no_points_on_vir_NMDs <- data.scores.all.vir[data.scores.all.vir$NMDS1 < 0.0041 | data.scores.all.vir$NMDS1 > 0.0059 | 
                                               data.scores.all.vir$NMDS2 < 0.004 | data.scores.all.vir$NMDS2 > 0.006, ]

no_points_on_host_NMDs <- data.scores.all.hosts[data.scores.all.hosts$NMDS1 < -0.003 | data.scores.all.hosts$NMDS1 > 0.009 | 
                                                  data.scores.all.hosts$NMDS2 < -0.0075 | data.scores.all.hosts$NMDS2 > 0.009, ]

write.table(no_points_on_vir_NMDs[c("Sample_name", "cohort", "Type", "Timepoint")], "../no_points_on_vir_NMDs.txt", sep='\t', row.names=F, col.names=T, quote=F)
write.table(no_points_on_host_NMDs[c("Sample_name", "cohort", "Type", "Timepoint")], "../no_points_on_host_NMDs.txt", sep='\t', row.names=F, col.names=T, quote=F)

#######################################################################################################################################

## Patching the Figure 1
#######################################################################################################################################
write.table(meta_working, "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/df_for_figures1ab.tsv", sep='\t', row.names=F, col.names=T, quote=F)
write.table(summarized_df2_frac_melt, "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/df_for_figure1c.tsv", sep='\t', row.names=F, col.names=T, quote=F)
write.table(summarized_df3_frac_melt, "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/df_for_figure1d.tsv", sep='\t', row.names=F, col.names=T, quote=F)
write.table(summarized_df4_frac_melt, "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/df_for_figures1e.tsv", sep='\t', row.names=F, col.names=T, quote=F)

meta_working$ncvssample <- factor(meta_working$ncvssample, levels = c("NCs", "SAMPLES"))
meta_all_with_qc_curated$ncvssample <- factor(meta_all_with_qc_curated$ncvssample, levels = c("NCs", "SAMPLES"))
meta_all_with_qc_curated$clean_reads_comb_for_plot <- meta_all_with_qc_curated$clean_reads_comb + (min(meta_all_with_qc_curated$clean_reads_comb[meta_all_with_qc_curated$clean_reads_comb > 0])/2)
meta_all_with_qc_curated$richness_for_plot <- meta_all_with_qc_curated$richness + (min(meta_all_with_qc_curated$richness[meta_all_with_qc_curated$richness > 0])/2)

figure_1A <- ggplot(meta_all_with_qc_curated[meta_all_with_qc_curated$clean_reads_comb > 0, ], aes(x=ncvssample, y=clean_reads_comb_for_plot)) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=1.5, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(alpha=0, outliers = FALSE) +
  facet_grid(. ~ cohort, labeller = labeller(
    cohort = c(
      "garmaeva" = "Garmaeva *et al.* <br>",
      "liang" = "Liang *et al.* <br>",
      "maqsood" = "Maqsood *et al.* <br>",
      "shah" = "Shah *et al.* <br>"))) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  scale_y_log10() +
  labs(y = expression(log[10] * "Number of clean reads"), tag="a", fill = "Timepoints") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=7),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.position = "bottom",
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=13)) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 3)
    )
  )

stat.test1 <- meta_working[meta_working$cohort!="garmaeva",] %>%
  group_by(cohort) %>%
  t_test(clean_reads_comb ~ ncvssample) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test1 <- stat.test1 %>% add_xy_position(x = "clean_reads_comb")
stat.test1[4,] <- stat.test1[3,]
stat.test1[4,"cohort"] <- "garmaeva"
stat.test1$y.position <- log10(stat.test1$y.position)
stat.test1[4,"y.position"] <- log10(max(meta_working[meta_working$cohort=="garmaeva",]$clean_reads_comb + 10000) )

stat.test1$xmin <- 1
stat.test1$xmax <- 2

stat.test1$p <- c(0.684, 0.0383, 0.725, NA)
stat.test1$p.adj <- c(0.725, 0.1149, 0.725, NA)
stat.test1$p.signif <- c("ns", "ns", "ns", "NA")

figure_1A <- figure_1A + stat_pvalue_manual(stat.test1, tip.length = 0.02, size=2.5, label = "p.signif")

figure_1B <- ggplot(meta_all_with_qc_curated[meta_all_with_qc_curated$clean_reads_comb > 0, ], aes(x=ncvssample, y=richness_for_plot)) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=1.5, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(alpha=0, outliers = FALSE) +
  facet_grid(. ~ cohort, labeller = labeller(
    cohort = c(
      "garmaeva" = "Garmaeva *et al.* <br>",
      "liang" = "Liang *et al.* <br>",
      "maqsood" = "Maqsood *et al.* <br>",
      "shah" = "Shah *et al.* <br>"))) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  scale_y_log10() +
  labs(y = expression(log[10] * "Richness"), tag="b", fill = "Timepoints") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=7),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=13),
    legend.position = "bottom") +
  guides(
    fill = guide_legend(
      override.aes = list(size = 3)
    )
  )

stat.test2 <- meta_working[meta_working$cohort!="garmaeva",] %>%
  group_by(cohort) %>%
  t_test(richness ~ ncvssample) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test2 <- stat.test2 %>% add_xy_position(x = "richness")
stat.test2[4,] <- stat.test2[3,]
stat.test2[4,"cohort"] <- "garmaeva"
stat.test2$y.position <- log10(stat.test2$y.position)
stat.test2[4,"y.position"] <- log10(max(meta_working[meta_working$cohort=="garmaeva",]$richness + 10000) )

stat.test2$xmin <- 1
stat.test2$xmax <- 2

stat.test2$p <- c(0.00244, 0.54, 0.000603, NA)
stat.test2$p.adj <- c(0.00366, 0.54, 0.001809, NA)
stat.test2$p.signif <- c("**", "ns", "**", "NA")

figure_1B <- figure_1B + stat_pvalue_manual(stat.test2, tip.length = 0.02, size=2.5, label = "p.signif")

figure_1C <- ggplot(data = data.scores.all.vir, aes(x = NMDS1, y = NMDS2, color=timepoint_type)) + 
  geom_point(size = 1.5, alpha=0.6) + 
  geom_point(data=centroids.vir, aes(fill=timepoint_type),shape=23, size=3, color='black', ) + 
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = timepoint_type, color=timepoint_type), linetype = 2) +
  xlim(0.0041,0.0059) +
  ylim(0.004,0.006) +
  theme_bw()+
  labs(tag="c", color = "Timepoints", title = "Bray-Curtis dissimilarity based on vOTUs") +
  scale_color_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  theme(
    plot.title = element_text(size = 8),
    strip.text = ggtext::element_markdown(size=7),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=13),
    legend.position = "bottom") +
  guides(fill = "none", color = "none") 
figure_1C <- ggMarginal(figure_1C, type="boxplot", groupFill=T)


figure_1D <- ggplot(data = data.scores.all.hosts, aes(x = NMDS1, y = NMDS2, color=timepoint_type)) + 
  geom_point(size = 1.5, alpha=0.6) + 
  geom_point(data=centroids.host, aes(fill=timepoint_type),shape=23, size=3, color='black', ) + 
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = timepoint_type, color=timepoint_type), linetype = 2) +
  xlim(-0.003,0.009) +
  ylim(-0.0075,0.009) +
  theme_bw()+
  labs(tag="d",color = "Timepoints", title = "Bray-Curtis dissimilarity based on host-based \nvOTU aggregates") +
  scale_color_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  theme(
    plot.title = element_text(size = 8),
    strip.text = ggtext::element_markdown(size=7),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=13),
    legend.position = "bottom") +
  guides(fill = "none", color = "none") 
figure_1D <- ggMarginal(figure_1D, type="boxplot", groupFill=T)

combined_plot <- wrap_elements(figure_1A + figure_1B + plot_layout(nrow=1, guides = "collect") & theme(legend.position = "bottom")) /
  (wrap_elements(figure_1C) | wrap_elements(figure_1D))

ggsave("Figure1_rebuttal_edited.pdf", combined_plot, width = 21/2.54, height = 24/2.54)

summarized_df2_frac_stats$ncvssample <- factor(summarized_df2_frac_stats$ncvssample, levels = c("NCs", "SAMPLES"))
supplementary_figure_4A <- ggplot(summarized_df2_frac_stats[summarized_df2_frac_stats$variable != "Unclassified", ], aes(x = variable, y = value, fill = ncvssample, color = ncvssample)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
              size = 1.5, shape = 19, stroke = 0.1, alpha=0.3) +
  geom_boxplot(alpha = 0, outliers = FALSE, position = position_dodge(width = 0.75)) +
  theme_bw() +
  scale_fill_manual(values = c("#619b8a", "#233d4d")) +
  scale_color_manual(values = c("#619b8a", "#233d4d")) +
  scale_y_continuous(limits = c(-0.05, 1.15), breaks = seq(0, 1, by = 0.25), expand = c(0, 0)) +
  facet_grid(. ~ cohort, labeller = labeller(
    cohort = c(
      "garmaeva" = "Garmaeva *et al.* <br>",
      "liang" = "Liang *et al.* <br>",
      "maqsood" = "Maqsood *et al.* <br>",
      "shah" = "Shah *et al.* <br>"))) +
  labs(x = "", y = "Fraction of the viral group richness \n (vOTUs with at least 50% completeness)", fill = "Viral group", tag="a") +
  theme(
    strip.text = ggtext::element_markdown(size=9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 10),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=12),
    legend.position = "none"
  ) +
  guides(color = "none", fill = "none")

summarized_df3_frac_stats$ncvssample <- factor(summarized_df3_frac_stats$ncvssample, levels = c("NCs", "SAMPLES"))
supplementary_figure_4B <- ggplot(summarized_df3_frac_stats[summarized_df3_frac_stats$variable != "Unclassified", ], aes(x = variable, y = value, fill = ncvssample, color = ncvssample)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
              size = 1.5, shape = 19, stroke = 0.1, alpha=0.3) +
  geom_boxplot(alpha = 0, outliers = FALSE, position = position_dodge(width = 0.75)) +
  theme_bw() +
  scale_fill_manual(values = c("#619b8a", "#233d4d")) +
  scale_color_manual(values = c("#619b8a", "#233d4d")) +
  scale_y_continuous(limits = c(-0.05, 1.15), breaks = seq(0, 1, by = 0.25), expand = c(0, 0)) +
  facet_grid(. ~ cohort, labeller = labeller(
    cohort = c(
      "garmaeva" = "Garmaeva *et al.* <br>",
      "liang" = "Liang *et al.* <br>",
      "maqsood" = "Maqsood *et al.* <br>",
      "shah" = "Shah *et al.* <br>"))) +
  labs(x = "", y = "Fraction of the host group richness \n (vOTUs with at least 50% completeness)", fill = "Host group", tag="b")+
  theme(
    strip.text = ggtext::element_markdown(size=9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 10),
    legend.title = element_text(size = 0),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, "cm"),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=12),
    legend.position = "bottom"
  ) +
  guides(fill = "none")

supplementary_figure4 <- supplementary_figure_4A + supplementary_figure_4B + plot_layout(nrow=2, guides = "collect") & theme(legend.position = "bottom")

ggsave("supplementary_figure4.pdf", supplementary_figure4, width = 21/2.54, height = 22/2.54)

meta_all_with_qc_curated$clean_reads_comb_for_plot <- NULL
meta_all_with_qc_curated$richness_for_plot <- NULL
meta_working$ncvssample <- factor(meta_working$ncvssample, levels = c("SAMPLES", "NCs"))
meta_all_with_qc_curated$ncvssample <- factor(meta_all_with_qc_curated$ncvssample, levels = c("SAMPLES", "NCs"))

#######################################################################################################################################

## Stats calculation for part 1
#######################################################################################################################################
summary(lmer(clean_reads_comb ~  ncvssample + (1|cohort/nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$clean_reads_comb > 0, ]))
summary(lmer(clean_reads_comb ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$clean_reads_comb > 0 & 
                                                                                                                meta_all_with_qc_curated$cohort == "liang", ]))
summary(lmer(clean_reads_comb ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$cohort == "maqsood", ]))
summary(lmer(clean_reads_comb ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$cohort == "shah", ]))
p.adjust(c(0.75, 0.0383, 0.725), method="BH")

summary(lmer(contigs_1000 ~ ncvssample + (1|cohort/nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$clean_reads_comb > 0, ]))
summary(lmer(contigs_1000 ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$clean_reads_comb > 0 & 
                                                                                                            meta_all_with_qc_curated$cohort == "liang", ]))
summary(lmer(contigs_1000 ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$cohort == "maqsood", ]))
summary(lmer(contigs_1000 ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$cohort == "shah", ]))
p.adjust(c(0.0321, 0.802, 0.360), method="BH")

summary(lmer(total_viruses_discovered ~ ncvssample + (1|cohort/nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$clean_reads_comb > 0, ]))
summary(lmer(total_viruses_discovered ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$clean_reads_comb > 0 & 
                                                                                                                        meta_all_with_qc_curated$cohort == "liang", ]))
summary(lmer(total_viruses_discovered ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$cohort == "maqsood", ]))
summary(lmer(total_viruses_discovered ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$cohort == "shah", ]))
p.adjust(c(0.000231, 0.645, 0.368), method="BH")

summary(lmer(CHM_LU_richness_discovered_ratio ~ ncvssample + (1|cohort/nc_subject_group), REML = F, data = meta_working))
summary(lmer(CHM_LU_richness_discovered_ratio ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_working[meta_working$cohort == "liang", ]))
summary(lmer(CHM_LU_richness_discovered_ratio ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_working[meta_working$cohort == "maqsood", ]))
summary(lmer(CHM_LU_richness_discovered_ratio ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_working[meta_working$cohort == "shah", ]))
p.adjust(c(0.0742, 0.020469, 0.8654), method="BH")

summary(lmer(richness ~ ncvssample + (1|cohort/nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$clean_reads_comb > 0, ]))
summary(lmer(richness ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$clean_reads_comb > 0 & 
                                                                                                        meta_all_with_qc_curated$cohort == "liang", ]))
summary(lmer(richness ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_working[meta_working$cohort == "maqsood", ]))
summary(lmer(richness ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_working[meta_working$cohort == "shah", ]))
p.adjust(c(6.53e-05, 0.540, 0.000603), method="BH")

summary(lmer(diversity ~ ncvssample + (1|cohort/nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$clean_reads_comb > 0, ]))
summary(lmer(diversity ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_all_with_qc_curated[meta_all_with_qc_curated$clean_reads_comb > 0 & 
                                                                                                         meta_all_with_qc_curated$cohort == "liang", ]))
summary(lmer(diversity ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_working[meta_working$cohort == "maqsood", ]))
summary(lmer(diversity ~  ncvssample + (1|nc_subject_group), REML = F, data = meta_working[meta_working$cohort == "shah", ]))
p.adjust(c(2.79e-06, 0.0696, 0.00723), method="BH")

calculate_mean_ci <- function(data_vector, confidence_level = 0.95) {
  data_vector <- na.omit(data_vector)
  mean_value <- mean(data_vector)
  standard_error <- sd(data_vector) / sqrt(length(data_vector))
  alpha <- 1 - confidence_level
  t_value <- qt(1 - alpha/2, df = length(data_vector) - 1)
  confidence_interval <- t_value * standard_error
  lower_bound <- mean_value - confidence_interval
  upper_bound <- mean_value + confidence_interval
  return(list(mean = mean_value, lower_bound = lower_bound, upper_bound = upper_bound))
}

result <- calculate_mean_ci(meta_working$richness[meta_working$ncvssample == "SAMPLES"])
cat("Mean:", result$mean, "\n")
cat("95% Confidence Interval: [", result$lower_bound, ", ", result$upper_bound, "]\n")

result <- calculate_mean_ci(meta_working$richness[meta_working$ncvssample == "NCs"])
cat("Mean:", result$mean, "\n")
cat("95% Confidence Interval: [", result$lower_bound, ", ", result$upper_bound, "]\n")

# For cohort-dependency in the number of discovered viruses in the NCs
model_richness_all <- lmer(richness ~ 1 + (1|cohort), REML = FALSE, data = meta_working)
variance_components <- VarCorr(model_richness_all)
var_cohort <- as.numeric(variance_components$cohort[1])
var_residual <- attr(variance_components, "sc")^2

icc <- var_cohort / (var_cohort + var_residual)
icc


model_richness_NCs <- lmer(richness ~ 1 + (1|cohort), REML = FALSE, data = meta_working[meta_working$ncvssample == "NCs", ])
variance_components <- VarCorr(model_richness_NCs)
var_cohort <- as.numeric(variance_components$cohort[1])
var_residual <- attr(variance_components, "sc")^2

icc <- var_cohort / (var_cohort + var_residual)
icc

model_richness_samples <- lmer(richness ~ 1 + (1|cohort), REML = FALSE, data = meta_working[meta_working$ncvssample == "SAMPLES", ])
variance_components <- VarCorr(model_richness_samples)
var_cohort <- as.numeric(variance_components$cohort[1])
var_residual <- attr(variance_components, "sc")^2

icc <- var_cohort / (var_cohort + var_residual)
icc

summary(lmer(total_viruses_discovered ~ 1 + (1|cohort), REML = FALSE, data = meta_all_with_qc_curated[meta_all_with_qc_curated$ncvssample == "NCs", ]))


# For age and cohort-dependency in the number of discovered viruses in the NCs
model_richness_samples_age <- lmer(richness ~ 1 + (1|cohort/Timepoint_numeric), REML = FALSE, data = meta_working[meta_working$ncvssample == "SAMPLES", ])
variance_components <- VarCorr(model_richness_samples_age)
var_cohort <- as.numeric(variance_components$cohort[1])
var_residual <- attr(variance_components, "sc")^2

icc <- var_cohort / (var_cohort + var_residual)
icc


#######################################################################################################################################

## Bray-Curtis boxplots: all samples; log scale
#######################################################################################################################################

# bray_dist_matrix_full <- as.matrix(vegdist(t(RPKM), method="bray"))
# bray_dist_matrix_full_rev <- 1 - bray_dist_matrix_full
# write.table(bray_dist_matrix_full_rev, "bray_dist_matrix_full_rev.tsv", sep='\t', row.names=T, col.names=T, quote=F)

bray_dist_matrix_full_rev <- read.delim('bray_dist_matrix_full_rev.tsv')
bray_dist_matrix_full <- 1 - bray_dist_matrix_full_rev

upper_tri <- upper.tri(bray_dist_matrix_full_rev)

# Get row and column indices for upper triangle
row_indices <- row(bray_dist_matrix_full_rev)[upper_tri]
col_indices <- col(bray_dist_matrix_full_rev)[upper_tri]

# Extract distances corresponding to upper triangle indices
distances <- bray_dist_matrix_full_rev[upper_tri]

# Create a dataframe
df_distances <- data.frame(
  Sample1 = rownames(bray_dist_matrix_full_rev)[row_indices],
  Sample2 = colnames(bray_dist_matrix_full_rev)[col_indices],
  Distance = distances
)

# Ensure Sample1 < Sample2 in each row (order does not matter)
df_distances <- df_distances[order(pmin(df_distances$Sample1, df_distances$Sample2)),
]  # Sorting by Sample1, Sample2

meta_working <- as.data.frame(meta_working)
meta_working$Sample1 <- meta_working$Sample_name
meta_working$Sample2 <- meta_working$Sample_name

df_distances <- merge(df_distances, meta_working[c("Sample1", "ncvssample", "cohort")], by="Sample1", all.x=T)
df_distances <- merge(df_distances, meta_working[c("Sample2", "ncvssample", "cohort")], by="Sample2", all.x=T)
colnames(df_distances) <- c("Sample1", "Sample2", "Distance", "Sample1_ncvssample",  "Sample1_cohort", "Sample2_ncvssample",  "Sample2_cohort")

df_distances_nc <- df_distances[!(df_distances$Sample1_ncvssample == "SAMPLES" & df_distances$Sample2_ncvssample == "SAMPLES"), ]
df_distances_nc <- df_distances_nc %>%
  mutate(Sample1_cohort = as.character(Sample1_cohort),
         Sample2_cohort = as.character(Sample2_cohort)) %>%
  mutate(category = ifelse(Sample1_ncvssample == "NCs" & Sample2_ncvssample == "SAMPLES", "Samples", NA),
         category = ifelse(Sample2_ncvssample == "NCs" & Sample1_ncvssample == "SAMPLES", "Samples", category),
         category = ifelse(Sample2_ncvssample == "NCs" & Sample1_ncvssample == "NCs", "NCs", category),
         type_cohort = ifelse(Sample1_cohort == Sample2_cohort, "same", "different"),
         cohort_nc = ifelse(type_cohort == "same" & category == "NCs", Sample1_cohort, NA),
         cohort_nc = ifelse(category == "Samples" & Sample1_ncvssample == "NCs", Sample1_cohort, cohort_nc),
         cohort_nc = ifelse(category == "Samples" & Sample2_ncvssample == "NCs", Sample2_cohort, cohort_nc))

switched_df <- df_distances_nc %>%
  filter(type_cohort == "different" & category == "NCs") %>%
  mutate(Sample3 = Sample1, 
         Sample1 = Sample2,
         Sample2 = Sample3,
         Sample3 = NULL,
         Sample3_cohort = Sample1_cohort,
         Sample1_cohort = Sample2_cohort,
         Sample2_cohort = Sample3_cohort,
         Sample3_cohort = NULL)

df_distances_nc <- rbind(df_distances_nc, switched_df)

df_distances_nc <- df_distances_nc %>%
  mutate(cohort_nc = ifelse(category == "NCs" & type_cohort == "different",  Sample1_cohort, cohort_nc))

# Permutation analysis


# Initialize a list to store the results for each combination
results <- list()

# Define the unique values in cohort_nc
cohort_values <- c("garmaeva", "liang", "maqsood", "shah")

# Loop over each cohort_nc value
for (cohort in cohort_values) {
  
  # Filter the dataframe for the current cohort
  df_filtered <- df_distances_nc[df_distances_nc$cohort_nc == cohort, ]
  
  # Apply the two conditions for Sample1_ncvssample
  conditions <- list(
    "Condition_1" = df_filtered[df_filtered$Sample1_ncvssample == df_filtered$Sample2_ncvssample, ],
    "Condition_2" = df_filtered[df_filtered$Sample1_ncvssample != df_filtered$Sample2_ncvssample, ]
  )
  
  # Loop over each condition
  for (condition_name in names(conditions)) {
    
    # Skip the combination cohort == "garmaeva" and Condition_1
    if (cohort == "garmaeva" && condition_name == "Condition_1") {
      next
    }
    
    # Get the filtered data for the current condition
    df_condition <- conditions[[condition_name]]
    
    # Perform the initial Wilcoxon test to get the baseline p-value
    baseline_result <- wilcox.test(Distance ~ type_cohort, data = df_condition)
    baseline_pvalue <- baseline_result$p.value
    
    # Set up the permutation test
    set.seed(123)  # Set seed for reproducibility
    n_permutations <- 10000
    ppermute <- numeric(n_permutations)
    
    # Perform permutations
    for (i in 1:n_permutations) {
      # Shuffle the "Distance" column
      shuffled_distances <- sample(df_condition$Distance)
      
      # Create a new dataframe with the shuffled "Distance" column
      df_shuffled <- df_condition
      df_shuffled$Distance <- shuffled_distances
      
      # Perform Wilcoxon test on the shuffled data
      result <- wilcox.test(Distance ~ type_cohort, data = df_shuffled)
      
      # Store the p-value in the ppermute vector
      ppermute[i] <- result$p.value
    }
    
    # Calculate the final p-value
    final_pvalue <- sum(ppermute <= baseline_pvalue) / n_permutations
    
    # Save the results in the list
    results[[paste(cohort, condition_name, sep = "_")]] <- list(
      baseline_pvalue = baseline_pvalue,
      final_pvalue = final_pvalue
    )
  }
}

# Output the results
results

#CHECK FOR 0 VALUES FIRST!!!
min_nonzero <- min(df_distances_nc$Distance[df_distances_nc$Distance > 0])
df_distances_nc_for_plot <- df_distances_nc
df_distances_nc_for_plot$Distance <- df_distances_nc_for_plot$Distance + (min_nonzero / 2)

df_distances_nc_for_plot$type_cohort <- factor(df_distances_nc_for_plot$type_cohort, levels = c("same", "different"))
df_distances_nc_for_plot$cohort_nc <- factor(df_distances_nc_for_plot$cohort_nc, levels = c("garmaeva", "liang", "maqsood", "shah"))

write.table(df_distances_nc_for_plot, "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/df_for_figure2c.tsv", sep='\t', row.names=F, col.names=T, quote=F)
#######################################################################################################################################

## Additional analysis - NCs vs samples
#######################################################################################################################################
df_distances_NS_SS <- data.frame(
  Sample1 = rownames(bray_dist_matrix_full_rev)[row_indices],
  Sample2 = colnames(bray_dist_matrix_full_rev)[col_indices],
  Distance = distances
)

df_distances_NS_SS <- merge(df_distances_NS_SS, meta_working[c("Sample1", "ncvssample", "cohort", "Subject_ID", "FAM_ID")], by="Sample1", all.x=T)
df_distances_NS_SS <- merge(df_distances_NS_SS, meta_working[c("Sample2", "ncvssample", "cohort", "Subject_ID", "FAM_ID")], by="Sample2", all.x=T)
colnames(df_distances_NS_SS) <- c("Sample1", "Sample2", "Distance",
                                  "Sample1_ncvssample",  "Sample1_cohort", "Sample1_Subject_ID", "Sample1_FAM_ID",
                                  "Sample2_ncvssample",  "Sample2_cohort", "Sample2_Subject_ID", "Sample2_FAM_ID")

df_distances_NS_SS <- df_distances_NS_SS[df_distances_NS_SS$Sample1_cohort == df_distances_NS_SS$Sample2_cohort, ]
df_distances_NS_SS <- df_distances_NS_SS[!(df_distances_NS_SS$Sample1_ncvssample == "NCs" & df_distances_NS_SS$Sample2_ncvssample == "NCs"), ]
df_distances_NS_SS <- df_distances_NS_SS[!((!is.na(df_distances_NS_SS$Sample1_FAM_ID)) & 
                                             (!is.na(df_distances_NS_SS$Sample2_FAM_ID)) & 
                                             df_distances_NS_SS$Sample1_FAM_ID == df_distances_NS_SS$Sample2_FAM_ID), ]

df_distances_NS_SS <- df_distances_NS_SS %>%
  mutate(category = ifelse(Sample1_ncvssample == Sample2_ncvssample, "samples", "ncs"))

df_distances_NS_SS$Distance_log <- df_distances_NS_SS$Distance + (min(df_distances_NS_SS$Distance[df_distances_NS_SS$Distance > 0]) / 2)

df_distances_NS_SS$category <- as.factor(df_distances_NS_SS$category)

figure_NS_SS_log <- ggplot(df_distances_NS_SS, aes(x = category, y = Distance_log )) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outlier.alpha = 0) +
  scale_y_log10() +
  scale_x_discrete(labels=c("ncs" = "NCs and Samples", "samples" = "Samples and Samples")) +
  labs(y = expression(log[10] * "Similarity index"), x = "", tag="c") +
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
ggsave("figure_NS_SS_log.pdf", figure_NS_SS_log, width = 14/2.54, height = 14/2.54)

# Do by cohort
figure_NS_SS_cohort <- ggplot(df_distances_NS_SS, aes(x = category, y = Distance )) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outlier.alpha = 0) +
  # scale_y_log10() +
  scale_x_discrete(labels=c("ncs" = "NCs and Samples", "samples" = "Samples and Samples")) +
  facet_wrap(~ Sample1_cohort, nrow = 2, ncol = 2, scales = "free", labeller = labeller(
    Sample1_cohort = c(
      "garmaeva" = "Samples Garmaeva *et al.* <br>",
      "liang" = "Samples Liang *et al.* <br>",
      "maqsood" = "Samples Maqsood *et al.* <br>",
      "shah" = "Samples Shah *et al.* <br>"
    )
  )) +
  labs(y = expression("Similarity index"), x = "", tag="c") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=7),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 7), 
    axis.text.y = element_text(size = 7),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = ggtext::element_markdown(face = "bold")
  )
ggsave("figure_NS_SS_cohort.pdf", figure_NS_SS_cohort, width = 14/2.54, height = 14/2.54)


# Permutation analysis: see detailed comments in the previous section
results <- list()
cohort_values <- c("garmaeva", "liang", "maqsood", "shah")

for (cohort in cohort_values) {
  df_filtered <- df_distances_NS_SS[df_distances_NS_SS$Sample1_cohort == cohort, ]
  baseline_result <- wilcox.test(Distance ~ category, data = df_filtered)
  baseline_pvalue <- baseline_result$p.value
  
  set.seed(123)
  n_permutations <- 1000
  ppermute <- numeric(n_permutations)
  
  for (i in 1:n_permutations) {
    print(i)
    print(cohort)
    shuffled_distances <- sample(df_filtered$Distance)
    df_shuffled <- df_filtered
    df_shuffled$Distance <- shuffled_distances
    result <- wilcox.test(Distance ~ category, data = df_shuffled)
    ppermute[i] <- result$p.value
    print(result$p.value)
  }
  
  final_pvalue <- sum(ppermute <= baseline_pvalue) / n_permutations
  
  results[[paste(cohort)]] <- list(
    baseline_pvalue = baseline_pvalue,
    final_pvalue = final_pvalue
  )
}

# Output the results
results

# write.table(df_distances_NS_SS, "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/df_for_figure2d_new.tsv", sep='\t', row.names=F, col.names=T, quote=F)
#######################################################################################################################################

## Analysis for the vOTUs sharedness with the other cohorts: correlation between sharedness with own vs different cohort
#######################################################################################################################################
RPKM_filtered_w_dummy <- RPKM_count

negative_controls_garmaeva <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "garmaeva" & 
                                                         meta_working$Sample_name %in% colnames(RPKM)]

negative_controls_liang <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "liang" &
                                                      meta_working$Sample_name %in% colnames(RPKM)]

negative_controls_maqsood <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "maqsood" &
                                                        meta_working$Sample_name %in% colnames(RPKM)]

negative_controls_shah <- meta_working$Sample_name[meta_working$Type == "Neg_ctrl" & meta_working$cohort == "shah" &
                                                     meta_working$Sample_name %in% colnames(RPKM)]

negative_controls <- c(negative_controls_garmaeva, negative_controls_liang, negative_controls_maqsood, negative_controls_shah)

RPKM_filtered_w_dummy$dummy_non_garmaeva_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% c(negative_controls_liang, negative_controls_maqsood, negative_controls_shah)])
RPKM_filtered_w_dummy$dummy_non_liang_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% c(negative_controls_garmaeva, negative_controls_maqsood, negative_controls_shah)])
RPKM_filtered_w_dummy$dummy_non_maqsood_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% c(negative_controls_garmaeva, negative_controls_liang, negative_controls_shah)])
RPKM_filtered_w_dummy$dummy_non_shah_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% c(negative_controls_garmaeva, negative_controls_liang, negative_controls_maqsood)])

RPKM_filtered_w_dummy$dummy_garmaeva_NC <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) == negative_controls_garmaeva]
RPKM_filtered_w_dummy$dummy_liang_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_liang])
RPKM_filtered_w_dummy$dummy_maqsood_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_maqsood])
RPKM_filtered_w_dummy$dummy_shah_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_shah])

# Save the result for updating eTOF later
writeLines(row.names(RPKM[RPKM_filtered_w_dummy$dummy_garmaeva_NC > 0, ]), con = "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/VIR_DB/sort_n_filter_for_upload/vOTUs_present_in_garmaeva_NCs")
writeLines(row.names(RPKM[RPKM_filtered_w_dummy$dummy_liang_NC > 0, ]), con = "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/VIR_DB/sort_n_filter_for_upload/vOTUs_present_in_liang_NCs")
writeLines(row.names(RPKM[RPKM_filtered_w_dummy$dummy_maqsood_NC > 0, ]), con = "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/VIR_DB/sort_n_filter_for_upload/vOTUs_present_in_maqsood_NCs")
writeLines(row.names(RPKM[RPKM_filtered_w_dummy$dummy_shah_NC > 0, ]), con = "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/VIR_DB/sort_n_filter_for_upload/vOTUs_present_in_shah_NCs")

RPKM_filtered_w_dummy[RPKM_filtered_w_dummy > 0] <- 1
RPKM_filtered_w_dummy <- RPKM_filtered_w_dummy[, !colnames(RPKM_filtered_w_dummy) %in% c(negative_controls)]

dummy_cols <- grep("^dummy", names(RPKM_filtered_w_dummy), value = TRUE)
non_dummy_cols <- setdiff(names(RPKM_filtered_w_dummy), dummy_cols)

# Initialize an empty dataframe for the results
result <- data.frame(matrix(ncol = length(dummy_cols), nrow = length(non_dummy_cols)))
names(result) <- dummy_cols
row.names(result) <- non_dummy_cols

# Create a data frame for present_in_both values
result_present_in_both <- data.frame(matrix(ncol = length(dummy_cols), nrow = length(non_dummy_cols)))
names(result_present_in_both) <- dummy_cols
row.names(result_present_in_both) <- non_dummy_cols

# Calculate the percent similarity and store present_in_both
for (non_dummy_col in non_dummy_cols) {
  for (dummy_col in dummy_cols) {
    # Number of viruses present in the non-dummy column
    present_in_non_dummy <- sum(RPKM_filtered_w_dummy[[non_dummy_col]] == 1)
    
    if (present_in_non_dummy > 0) {
      # Number of viruses present in both non-dummy and dummy column
      present_in_both <- sum(RPKM_filtered_w_dummy[[non_dummy_col]] == 1 & RPKM_filtered_w_dummy[[dummy_col]] == 1)
      
      # Calculate the percent similarity
      similarity <- present_in_both / present_in_non_dummy * 100
    } else {
      similarity <- NA  # or 0 if you prefer
      present_in_both <- NA  # or 0 if you prefer
    }
    
    # Store the similarity result
    result[non_dummy_col, dummy_col] <- similarity
    
    # Store the present_in_both value
    result_present_in_both[non_dummy_col, dummy_col] <- present_in_both
  }
}

result$Sample_name <- row.names(result)
row.names(result) <- NULL
table_for_plot_cor <- merge(result, meta_working[c("Sample_name", "Type", "cohort", "Subject_ID", "Timepoint", "nc_subject_group")], by="Sample_name", all.x = T)
table_for_plot_cor <- table_for_plot_cor %>%
  mutate(same_cohort_NC = ifelse(cohort == "maqsood", dummy_maqsood_NC, NA),
         same_cohort_NC = ifelse(cohort == "liang", dummy_liang_NC, same_cohort_NC),
         same_cohort_NC = ifelse(cohort == "shah", dummy_shah_NC, same_cohort_NC),
         same_cohort_NC = ifelse(cohort == "garmaeva", dummy_garmaeva_NC, same_cohort_NC),
         different_cohort_NC = ifelse(cohort == "maqsood", dummy_non_maqsood_NC, NA),
         different_cohort_NC = ifelse(cohort == "liang", dummy_non_liang_NC, different_cohort_NC),
         different_cohort_NC = ifelse(cohort == "shah", dummy_non_shah_NC, different_cohort_NC),
         different_cohort_NC = ifelse(cohort == "garmaeva", dummy_non_garmaeva_NC, different_cohort_NC)
  )

result_present_in_both$Sample_name <- row.names(result_present_in_both)
result_present_in_both <- merge(result_present_in_both, meta_working[c("Sample_name", "cohort")], by="Sample_name", all.x = T)
result_present_in_both <- result_present_in_both %>%
  mutate(same_cohort_NC = ifelse(cohort == "maqsood", dummy_maqsood_NC, NA),
         same_cohort_NC = ifelse(cohort == "liang", dummy_liang_NC, same_cohort_NC),
         same_cohort_NC = ifelse(cohort == "shah", dummy_shah_NC, same_cohort_NC),
         same_cohort_NC = ifelse(cohort == "garmaeva", dummy_garmaeva_NC, same_cohort_NC),
         different_cohort_NC = ifelse(cohort == "maqsood", dummy_non_maqsood_NC, NA),
         different_cohort_NC = ifelse(cohort == "liang", dummy_non_liang_NC, different_cohort_NC),
         different_cohort_NC = ifelse(cohort == "shah", dummy_non_shah_NC, different_cohort_NC),
         different_cohort_NC = ifelse(cohort == "garmaeva", dummy_non_garmaeva_NC, different_cohort_NC))


RPKM_filtered_w_dummy <- RPKM

RPKM_filtered_w_dummy$dummy_non_garmaeva_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% c(negative_controls_liang, negative_controls_maqsood, negative_controls_shah)])
RPKM_filtered_w_dummy$dummy_non_liang_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% c(negative_controls_garmaeva, negative_controls_maqsood, negative_controls_shah)])
RPKM_filtered_w_dummy$dummy_non_maqsood_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% c(negative_controls_garmaeva, negative_controls_liang, negative_controls_shah)])
RPKM_filtered_w_dummy$dummy_non_shah_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% c(negative_controls_garmaeva, negative_controls_liang, negative_controls_maqsood)])

RPKM_filtered_w_dummy$dummy_garmaeva_NC <- RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) == negative_controls_garmaeva]
RPKM_filtered_w_dummy$dummy_liang_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_liang])
RPKM_filtered_w_dummy$dummy_maqsood_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_maqsood])
RPKM_filtered_w_dummy$dummy_shah_NC <- rowSums(RPKM_filtered_w_dummy[, colnames(RPKM_filtered_w_dummy) %in% negative_controls_shah])

RPKM_filtered_w_dummy <- RPKM_filtered_w_dummy[, !colnames(RPKM_filtered_w_dummy) %in% c(negative_controls)]

dummy_cols <- grep("^dummy", names(RPKM_filtered_w_dummy), value = TRUE)
non_dummy_cols <- setdiff(names(RPKM_filtered_w_dummy), dummy_cols)

# Initialize an empty dataframe for the results
result <- data.frame(matrix(ncol = length(dummy_cols), nrow = length(non_dummy_cols)))
names(result) <- dummy_cols
row.names(result) <- non_dummy_cols

# Calculate the percent similarity
for (non_dummy_col in non_dummy_cols) {
  for (dummy_col in dummy_cols) {
    # Number of viruses present in the non-dummy column
    present_in_non_dummy <- sum(RPKM_filtered_w_dummy[[non_dummy_col]])
    
    if (present_in_non_dummy > 0) {
      # Number of viruses present in both non-dummy and dummy column
      present_in_both <- sum(RPKM_filtered_w_dummy[[non_dummy_col]][RPKM_filtered_w_dummy[[non_dummy_col]] > 0 & RPKM_filtered_w_dummy[[dummy_col]] > 0])
      
      # Calculate the percent similarity
      similarity <- present_in_both / present_in_non_dummy * 100
    } else {
      similarity <- NA  # or 0 if you prefer
    }
    
    # Store the result
    result[non_dummy_col, dummy_col] <- similarity
  }
}


result$Sample_name <- row.names(result)
row.names(result) <- NULL
table_for_plot_cor_a <- merge(result, meta_working[c("Sample_name", "Type", "cohort", "Subject_ID", "Timepoint", "nc_subject_group")], by="Sample_name", all.x = T)
table_for_plot_cor_a <- table_for_plot_cor_a %>%
  mutate(same_cohort_NC = ifelse(cohort == "maqsood", dummy_maqsood_NC, NA),
         same_cohort_NC = ifelse(cohort == "liang", dummy_liang_NC, same_cohort_NC),
         same_cohort_NC = ifelse(cohort == "shah", dummy_shah_NC, same_cohort_NC),
         same_cohort_NC = ifelse(cohort == "garmaeva", dummy_garmaeva_NC, same_cohort_NC),
         different_cohort_NC = ifelse(cohort == "maqsood", dummy_non_maqsood_NC, NA),
         different_cohort_NC = ifelse(cohort == "liang", dummy_non_liang_NC, different_cohort_NC),
         different_cohort_NC = ifelse(cohort == "shah", dummy_non_shah_NC, different_cohort_NC),
         different_cohort_NC = ifelse(cohort == "garmaeva", dummy_non_garmaeva_NC, different_cohort_NC)
  )


table_for_plot_cor_combine <- merge(table_for_plot_cor, table_for_plot_cor_a[c("Sample_name", "same_cohort_NC", "different_cohort_NC")], by="Sample_name", 
                                    suffixes = c("_presence","_abundance"))

# Calculate CI
result <- calculate_mean_ci(table_for_plot_cor_combine$same_cohort_NC_abundance)
cat("Mean:", result$mean, "\n")
cat("95% Confidence Interval: [", result$lower_bound, ", ", result$upper_bound, "]\n")

result <- calculate_mean_ci(table_for_plot_cor_combine$same_cohort_NC_presence)
cat("Mean:", result$mean, "\n")
cat("95% Confidence Interval: [", result$lower_bound, ", ", result$upper_bound, "]\n")

result <- calculate_mean_ci(table_for_plot_cor_combine$different_cohort_NC_presence)
cat("Mean:", result$mean, "\n")
cat("95% Confidence Interval: [", result$lower_bound, ", ", result$upper_bound, "]\n")

spearman_result_garmaeva <- psych::corr.test(table_for_plot_cor_combine[table_for_plot_cor_combine$cohort == "garmaeva", c("different_cohort_NC_presence", "same_cohort_NC_presence")], method = "spearman")
print(spearman_result_garmaeva$r)
print(spearman_result_garmaeva$p)

spearman_result_maqsood <- psych::corr.test(table_for_plot_cor_combine[table_for_plot_cor_combine$cohort == "maqsood", c("different_cohort_NC_presence", "same_cohort_NC_presence")], method = "spearman")
print(spearman_result_maqsood$r)
print(spearman_result_maqsood$p)

spearman_result_liang <- psych::corr.test(table_for_plot_cor_combine[table_for_plot_cor_combine$cohort == "liang", c("different_cohort_NC_presence", "same_cohort_NC_presence")], method = "spearman")
print(spearman_result_liang$r)
print(spearman_result_liang$p)

spearman_result_shah <- psych::corr.test(table_for_plot_cor_combine[table_for_plot_cor_combine$cohort == "shah", c("different_cohort_NC_presence", "same_cohort_NC_presence")], method = "spearman")
print(spearman_result_shah$r)
print(spearman_result_shah$p)


table_for_plot_cor_combine$Type <- as.factor(table_for_plot_cor_combine$Type)
table_for_plot_cor_combine$Type <- factor(table_for_plot_cor_combine$Type, levels = c("Mother", "Infant"))

summary(lmer(same_cohort_NC_abundance ~ cohort + (1|nc_subject_group), REML = F, data = table_for_plot_cor_combine))

model_perc_shared_predic <- lmer(different_cohort_NC_presence ~ same_cohort_NC_presence + Type + cohort + (1 | nc_subject_group), data=table_for_plot_cor_combine)
summary(model_perc_shared_predic)
r_squared <- r.squaredGLMM(model_perc_shared_predic)
print(r_squared)


summary(lmer(same_cohort_NC_presence ~ Type + (1|nc_subject_group), data=table_for_plot_cor_combine[table_for_plot_cor_combine$cohort == "garmaeva", ]))
summary(lm(same_cohort_NC_presence ~ Type, data=table_for_plot_cor_combine[table_for_plot_cor_combine$cohort == "maqsood", ]))

table_for_plot_cor_combine_melt <- melt(table_for_plot_cor_combine[c("Sample_name", "cohort", "same_cohort_NC_presence", "different_cohort_NC_presence", 
                                                                     "same_cohort_NC_abundance", "different_cohort_NC_abundance")], 
                                        id.vars = c("Sample_name", "cohort"))

table_for_plot_cor_combine_melt <- table_for_plot_cor_combine_melt %>%
  mutate(pres_abun = ifelse(grepl("presence", variable), "presence", "abundance"),
         same_diff = ifelse(grepl("same", variable), "same", "different"))

table_for_plot_cor_combine_melt$same_diff <- factor(table_for_plot_cor_combine_melt$same_diff, levels = c("same", "different"))
table_for_plot_cor_combine_melt$pres_abun <- factor(table_for_plot_cor_combine_melt$pres_abun, levels = c("presence", "abundance"))

# write.table(table_for_plot_cor_combine, "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/df_for_figure2d.tsv", sep='\t', row.names=F, col.names=T, quote=F)
# write.table(table_for_plot_cor_combine_melt, "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/df_for_figure9as.tsv", sep='\t', row.names=F, col.names=T, quote=F)

# Permutation analysis


# Initialize a list to store the results for each combination
results <- list()

# Define the unique values in cohort_nc
cohort_values <- c("garmaeva", "liang", "maqsood", "shah")

# Loop over each cohort_nc value
for (cohort in cohort_values) {
  
  # Filter the dataframe for the current cohort
  df_filtered <- table_for_plot_cor_combine_melt[table_for_plot_cor_combine_melt$cohort == cohort, ]
  
  # Apply the two conditions for pres_abun
  conditions <- list(
    "Condition_1" = df_filtered[df_filtered$pres_abun =="presence", ],
    "Condition_2" = df_filtered[df_filtered$pres_abun == "abundance", ]
  )
  
  # Loop over each condition
  for (condition_name in names(conditions)) {
    
    # Get the filtered data for the current condition
    df_condition <- conditions[[condition_name]]
    
    # Perform the initial Wilcoxon test to get the baseline p-value
    baseline_result <- wilcox.test(value ~ same_diff, data = df_condition)
    baseline_pvalue <- baseline_result$p.value
    
    # Set up the permutation test
    set.seed(123)  # Set seed for reproducibility
    n_permutations <- 10000
    ppermute <- numeric(n_permutations)
    
    # Perform permutations
    for (i in 1:n_permutations) {
      # Shuffle the "value" column
      shuffled_values <- sample(df_condition$value)
      
      # Create a new dataframe with the shuffled "value" column
      df_shuffled <- df_condition
      df_shuffled$value <- shuffled_values
      
      # Perform Wilcoxon test on the shuffled data
      result <- wilcox.test(value ~ same_diff, data = df_shuffled)
      
      # Store the p-value in the ppermute vector
      ppermute[i] <- result$p.value
    }
    
    # Calculate the final p-value
    final_pvalue <- sum(ppermute <= baseline_pvalue) / n_permutations
    
    # Save the results in the list
    results[[paste(cohort, condition_name, sep = "_")]] <- list(
      baseline_pvalue = baseline_pvalue,
      final_pvalue = final_pvalue
    )
  }
}

# Output the results
results

table_for_plot_cor_combine_melt_subset <- table_for_plot_cor_combine[table_for_plot_cor_combine$cohort %in% c("garmaeva", "maqsood"), c("cohort", "Type", "same_cohort_NC_presence")]

#CHECK FOR 0 VALUES FIRST!!!
min_nonzero <- min(table_for_plot_cor_combine_melt_subset$same_cohort_NC_presence[table_for_plot_cor_combine_melt_subset$same_cohort_NC_presence > 0])
table_for_plot_cor_combine_melt_subset_log <- table_for_plot_cor_combine_melt_subset
table_for_plot_cor_combine_melt_subset_log$same_cohort_NC_presence <- table_for_plot_cor_combine_melt_subset_log$same_cohort_NC_presence + (min_nonzero / 2)
table_for_plot_cor_combine_melt_subset_log$Type <- factor(table_for_plot_cor_combine_melt_subset_log$Type, levels = c("Infant", "Mother"))
table_for_plot_cor_combine_melt_subset$Type <- factor(table_for_plot_cor_combine_melt_subset$Type, levels = c("Infant", "Mother"))
#######################################################################################################################################

## Analysis for the vOTUs sharedness along the timepoints: with dummy, only Liang and Garmaeva cohorts, no cross-cohort comparison
#######################################################################################################################################
filter_RPKM_table <- function(RPKM_data, samples) {
  filtered_RPKM <- RPKM_data[, colnames(RPKM_data) %in% samples]
  filtered_RPKM <- filtered_RPKM[rowSums(filtered_RPKM) > 0, colSums(filtered_RPKM) > 0]
  return(filtered_RPKM)
}


calculate_similarity <- function(RPKM_data, nc_samples, non_nc_samples) {
  result <- data.frame(matrix(ncol = length(nc_samples), nrow = length(non_nc_samples)))
  names(result) <- nc_samples
  row.names(result) <- non_nc_samples
  
  for (non_dummy_col in non_nc_samples) {
    for (dummy_col in nc_samples) {
      present_in_non_dummy <- sum(RPKM_data[[non_dummy_col]] == 1)
      
      if (present_in_non_dummy > 0) {
        present_in_both <- sum(RPKM_data[[non_dummy_col]] == 1 & RPKM_data[[dummy_col]] == 1)
        similarity <- present_in_both / present_in_non_dummy * 100
      } else {
        similarity <- NA
      }
      
      result[non_dummy_col, dummy_col] <- similarity
    }
  }
  
  result$Sample_name <- row.names(result)
  row.names(result) <- NULL
  return(result)
}

calculate_similarity_a <- function(RPKM_data, nc_samples, non_nc_samples) {
  result <- data.frame(matrix(ncol = length(nc_samples), nrow = length(non_nc_samples)))
  names(result) <- nc_samples
  row.names(result) <- non_nc_samples
  
  for (non_dummy_col in non_nc_samples) {
    for (dummy_col in nc_samples) {
      abundance_in_non_dummy <- sum(RPKM_data[[non_dummy_col]])
      
      if (abundance_in_non_dummy > 0) {
        abundace_present_in_both <- sum(RPKM_data[[non_dummy_col]][RPKM_data[[non_dummy_col]] > 0 & RPKM_data[[dummy_col]] > 0])  
        similarity <- abundace_present_in_both / abundance_in_non_dummy * 100
      } else {
        similarity <- NA
      }
      
      result[non_dummy_col, dummy_col] <- similarity
    }
  }
  
  result$Sample_name <- row.names(result)
  row.names(result) <- NULL
  return(result)
}

negative_controls <- meta_all_with_qc_curated$Sample_name[meta_all_with_qc_curated$Type == "Neg_ctrl" & !is.na(meta_all_with_qc_curated$Type)]
samples_g <- c(meta_working$Sample_name[meta_working$Type %in% c("Infant", "Neg_ctrl") & meta_working$cohort == "garmaeva"])

RPKM_filtered_g <- filter_RPKM_table(RPKM_count, samples_g)

nc_samples_g <- colnames(RPKM_filtered_g)[colnames(RPKM_filtered_g) %in% negative_controls]
non_nc_samples_g <- setdiff(names(RPKM_filtered_g), nc_samples_g)

result_g <- calculate_similarity(RPKM_filtered_g, nc_samples_g, non_nc_samples_g)
table_for_plot_g <- merge(result_g, meta_working[c("Sample_name", "Timepoint", "Subject_ID", "cohort", "total_viruses_discovered", "clean_reads_comb")], by="Sample_name", all.x = TRUE)
table_for_plot_g$Timepoint <- factor(table_for_plot_g$Timepoint, levels = c("M1", "M2", "M3", "M6", "M12"))
table_for_plot_g$Dataset <- "RPKM_count"
colnames(table_for_plot_g) <- c("Sample_name", "value", "Timepoint", "Subject_ID", "cohort", "total_viruses_discovered", "clean_reads_comb", "Dataset")
table_for_plot_g$Subject_ID <- as.factor(table_for_plot_g$Subject_ID)
table_for_plot_g$Timepoint_numeric <- as.integer(gsub("M", "", table_for_plot_g$Timepoint))

summary(lmer(value ~ Timepoint_numeric + (1|Subject_ID), REML = F, data = table_for_plot_g))

RPKM_filtered_ga <- filter_RPKM_table(RPKM, samples_g)

result_ga <- calculate_similarity_a(RPKM_filtered_ga, nc_samples_g, non_nc_samples_g)
table_for_plot_ga <- merge(result_ga, meta_working[c("Sample_name", "Timepoint", "Subject_ID", "cohort", "total_viruses_discovered", "clean_reads_comb")], by="Sample_name", all.x = TRUE)
table_for_plot_ga$Timepoint <- factor(table_for_plot_ga$Timepoint, levels = c("M1", "M2", "M3", "M6", "M12"))
table_for_plot_ga$Dataset <- "RPKM"
colnames(table_for_plot_ga) <- c("Sample_name", "value", "Timepoint", "Subject_ID", "cohort", "total_viruses_discovered", "clean_reads_comb", "Dataset")
table_for_plot_ga$Subject_ID <- as.factor(table_for_plot_ga$Subject_ID)
table_for_plot_ga$Timepoint_numeric <- as.integer(gsub("M", "", table_for_plot_ga$Timepoint))

summary(lmer(value ~ Timepoint_numeric + (1|Subject_ID), REML = F, data = table_for_plot_ga))

samples_l <- c(meta_working$Sample_name[(meta_working$Subcohort == "Discovery" & !is.na(meta_working$Subcohort) 
                                         & meta_working$cohort == "liang" & meta_working$Type == "Infant") | 
                                          (meta_working$cohort == "liang" & meta_working$Type == "Neg_ctrl")])

RPKM_filtered_l <- filter_RPKM_table(RPKM_count, samples_l)

RPKM_filtered_l$liang_nc_dummy <- rowSums(RPKM_filtered_l[, colnames(RPKM_filtered_l) %in% negative_controls], na.rm = TRUE)
RPKM_filtered_l$liang_nc_dummy[RPKM_filtered_l$liang_nc_dummy > 0] <- 1
RPKM_filtered_l <- RPKM_filtered_l[, !colnames(RPKM_filtered_l) %in% negative_controls]

nc_samples_l <- "liang_nc_dummy"
non_nc_samples_l <- setdiff(names(RPKM_filtered_l), nc_samples_l)

result_l <- calculate_similarity(RPKM_filtered_l, nc_samples_l, non_nc_samples_l)

table_for_plot_l <- merge(result_l, meta_working[c("Sample_name", "Sample_ID", "Timepoint", "Subject_ID", "cohort", "total_viruses_discovered", "clean_reads_comb")], by="Sample_name", all.x = T)
table_for_plot_l <- table_for_plot_l %>%
  mutate(nucleic_acid = ifelse(grepl("_D", Sample_ID), "DNA", "RNA"))
table_for_plot_l <- table_for_plot_l[table_for_plot_l$nucleic_acid == "DNA", ]
table_for_plot_l$nucleic_acid <- NULL
table_for_plot_l$Sample_ID <- NULL
table_for_plot_l$Dataset <- "RPKM_count"
table_for_plot_l$Timepoint <- factor(table_for_plot_l$Timepoint, levels = c("M0", "M1", "M4"))
colnames(table_for_plot_l) <- c("Sample_name", "value", "Timepoint", "Subject_ID", "cohort", "total_viruses_discovered", "clean_reads_comb", "Dataset")
table_for_plot_l$Subject_ID <- as.factor(table_for_plot_l$Subject_ID)
table_for_plot_l$Timepoint_numeric <- as.integer(gsub("M", "", table_for_plot_l$Timepoint))

summary(lmer(value ~ Timepoint_numeric + (1|Subject_ID), REML = F, data = table_for_plot_l))

RPKM_filtered_la <- filter_RPKM_table(RPKM, samples_l)

RPKM_filtered_la$liang_nc_dummy <- rowSums(RPKM_filtered_la[, colnames(RPKM_filtered_la) %in% negative_controls], na.rm = TRUE)
RPKM_filtered_la$liang_nc_dummy[RPKM_filtered_la$liang_nc_dummy > 0] <- 1
RPKM_filtered_la <- RPKM_filtered_la[, !colnames(RPKM_filtered_la) %in% negative_controls]

result_la <- calculate_similarity_a(RPKM_filtered_la, nc_samples_l, non_nc_samples_l)
table_for_plot_la <- merge(result_la, meta_working[c("Sample_name", "Sample_ID", "Timepoint", "Subject_ID", "cohort", "total_viruses_discovered", "clean_reads_comb")], by="Sample_name", all.x = T)
table_for_plot_la <- table_for_plot_la %>%
  mutate(nucleic_acid = ifelse(grepl("_D", Sample_ID), "DNA", "RNA"))
table_for_plot_la <- table_for_plot_la[table_for_plot_la$nucleic_acid == "DNA", ]
table_for_plot_la$nucleic_acid <- NULL
table_for_plot_la$Sample_ID <- NULL
table_for_plot_la$Dataset <- "RPKM"
table_for_plot_la$Timepoint <- factor(table_for_plot_la$Timepoint, levels = c("M0", "M1", "M4"))
colnames(table_for_plot_la) <- c("Sample_name", "value", "Timepoint", "Subject_ID", "cohort", "total_viruses_discovered", "clean_reads_comb", "Dataset")
table_for_plot_la$Subject_ID <- as.factor(table_for_plot_la$Subject_ID)
table_for_plot_la$Timepoint_numeric <- as.integer(gsub("M", "", table_for_plot_la$Timepoint))

summary(lmer(value ~ Timepoint_numeric + (1|Subject_ID), REML = F, data = table_for_plot_la))

combined_table_for_plot <- rbind(table_for_plot_g, table_for_plot_ga, table_for_plot_l, table_for_plot_la)
combined_table_for_plot$Timepoint <- factor(combined_table_for_plot$Timepoint, levels = c("M0", "M1", "M2", "M3", "M4", "M6", "M12"))
combined_table_for_plot$Dataset <- factor(combined_table_for_plot$Dataset, levels = c("RPKM_count", "RPKM"))

#CHECK FOR 0 VALUES FIRST!!!
min_nonzero <- min(combined_table_for_plot$value[combined_table_for_plot$value > 0])
# combined_table_for_plot_log <- combined_table_for_plot
# combined_table_for_plot_log$value <- combined_table_for_plot_log$value + (min_nonzero / 2)
# 
write.table(combined_table_for_plot, "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/df_for_figure3c.tsv", sep='\t', row.names=F, col.names=T, quote=F)
#######################################################################################################################################

## Additional analysis: absolute number calculation
#######################################################################################################################################
result_present_in_both_mi <- result_present_in_both[result_present_in_both$cohort %in% c("maqsood", "garmaeva"), c("Sample_name", "cohort", "same_cohort_NC")]
result_present_in_both_mi <- merge(result_present_in_both_mi, meta_working[c("Sample_name", "Type", "Subject_ID", "Timepoint", "nc_subject_group")], by="Sample_name", all.x = T)

summary(lmer(same_cohort_NC ~ Type + (1|nc_subject_group), data=result_present_in_both_mi[result_present_in_both_mi$cohort == "garmaeva", ]))
summary(lm(same_cohort_NC ~ Type, data=result_present_in_both_mi[result_present_in_both_mi$cohort == "maqsood", ]))

stat.testaaan <- result_present_in_both_mi %>%
  group_by(cohort) %>%
  t_test(same_cohort_NC ~ Type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.testaaan <- stat.testaaan %>% add_xy_position(x = "same_cohort_NC")
stat.testaaan$xmin <- 1
stat.testaaan$xmax <- 2
stat.testaaan$p.signif <- c("***", "***")

pdf('absolute_count_contaminants_mom_vs_inf.pdf', width=12/2.54, height=8/2.54)
ggplot(result_present_in_both_mi, aes(x = Type, y = same_cohort_NC)) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outlier.alpha = 0) +
  facet_grid(. ~ cohort, scales = "free", labeller = labeller(
    cohort = c(
      "garmaeva" = "Samples Garmaeva *et al.* <br>",
      "maqsood" = "Samples Maqsood *et al.* <br>"
    )
  )) +
  labs(y = "N vOTUs shared with NCs") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=8),
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    strip.background = element_rect(fill = "transparent")
  ) + 
  stat_pvalue_manual(stat.testaaan, tip.length = 0.02, size=2.5, label = "p.signif")
dev.off()
#######################################################################################################################################

## Strain check
#######################################################################################################################################
strains_df <- strains_df_ini
strains_df <- strains_df_ini[!is.na(strains_df_ini$popANI), ]

popANI_distribution <- ggplot(data=strains_df, aes(x = popANI)) +
  geom_histogram(bins = 200, fill = "red") +
  labs(x = "population ANI", y = "N viral ") +
  xlim(0.998, 1) +
  ylim(0, 7000) +
  theme_bw()
#######################################################################################################################################

## Diversity association
#######################################################################################################################################
meta_working <- meta_working %>%
  mutate(Timepoint_numeric = as.integer(gsub("M", "", Timepoint)),
         Timepoint_numeric = ifelse(grepl("Y2-5", Timepoint), 42, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother", 384, Timepoint_numeric))

table_div_percshared_asoc <- merge(table_for_plot_cor_combine_melt[table_for_plot_cor_combine_melt$pres_abun == "presence" 
                                                                   & table_for_plot_cor_combine_melt$same_diff == "same",],
                                   meta_working[c("Sample_name", "diversity", "Type", "Timepoint_numeric", "nc_subject_group", "timepoint_type", "Timepoint")], by="Sample_name", all.x=T)

summary(lmer(value ~ diversity + (1|cohort/nc_subject_group), REML = F, data = table_div_percshared_asoc))
summary(lmer(value ~ diversity + Type + (1|cohort/nc_subject_group), REML = F, data = table_div_percshared_asoc))
summary(lmer(value ~ diversity + Type + Timepoint_numeric + (1|cohort/nc_subject_group), REML = F, data = table_div_percshared_asoc))

table_div_percshared_asoc <- table_div_percshared_asoc %>%
  mutate(timepoint_type = ifelse(grepl("Infant", timepoint_type), Timepoint, timepoint_type))
table_div_percshared_asoc$timepoint_type <- as.factor(table_div_percshared_asoc$timepoint_type)
table_div_percshared_asoc$timepoint_type <- factor(table_div_percshared_asoc$timepoint_type, 
                                                   levels = c("M0", "M1", "M2", "M3", "M4", "M6", "M12", "Y2-5", "Mother"))
#######################################################################################################################################

## Supplementary information
#######################################################################################################################################
filtered_extended_tof <- extended_tof[extended_tof$New_CID %in% row.names(RPKM), ]

figure_2A_supp <- ggplot() + 
  geom_histogram(data = filtered_extended_tof, aes(POST_CHV_length, color="All", fill="All"), alpha = 0.2, bins=60) + 
  geom_histogram(data = filtered_extended_tof[filtered_extended_tof$checkv_quality=="Not-determined",], aes(POST_CHV_length, color="Not-determined", fill="Not-determined"), alpha = 0.2, bins=60) +
  geom_histogram(data = filtered_extended_tof[filtered_extended_tof$checkv_quality=="Low-quality",], aes(POST_CHV_length, color="Low-quality", fill="Low-quality"), alpha = 0.2, bins=60) +
  geom_histogram(data = filtered_extended_tof[filtered_extended_tof$checkv_quality=="Medium-quality",], aes(POST_CHV_length, color="Medium-quality", fill="Medium-quality"), alpha = 0.2, bins=60) +
  geom_histogram(data = filtered_extended_tof[filtered_extended_tof$checkv_quality=="High-quality",], aes(POST_CHV_length, color="High-quality", fill="High-quality"), alpha = 0.2, bins=60) +
  geom_histogram(data = filtered_extended_tof[filtered_extended_tof$checkv_quality=="Complete",], aes(POST_CHV_length, color="Complete", fill="Complete"), alpha = 0.2, bins=60) +
  labs(x="Virus contig length, bp", y="log10(N virus contigs \n or fragments)", fill="Genome Quality", color="Genome Quality", tag="a") +
  scale_color_manual(breaks=c("All", "Not-determined", "Low-quality", "Medium-quality", "High-quality", "Complete"), 
                     values=c("#adadad", "#FDE725FF", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF")) + 
  scale_fill_manual(breaks=c("All", "Not-determined", "Low-quality", "Medium-quality", "High-quality", "Complete"), 
                    values=c("#DDDDDD", "#FDE725FF", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF")) +
  scale_x_log10(breaks=c(100, 1000, 10000, 100000, 1000000), labels=c("100", "1,000", "10,000", "100,000", "1,000,000")) +
  scale_y_log10() +
  theme_bw() + 
  theme(axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size=8),
        legend.text = element_text(size = 6),
        legend.position = "bottom",
        plot.tag = element_text(face="bold", size=6),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-3,-10)) +
  guides(color = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1, keywidth = 1, keyheight = 1)) + 
  guides(fill = guide_legend(title.position = "top",label.position = "left", title.hjust = 0.5, nrow = 1, keywidth = 1, keyheight = 1))

filtered_extended_tof$virus_host_ratio <- (filtered_extended_tof$viral_genes) / (filtered_extended_tof$host_genes)

figure_2B_supp <- ggplot(data=filtered_extended_tof, aes(x = virus_host_ratio)) +
  geom_histogram(aes(y = after_stat(count)), bins = 30, fill = "red", alpha = 0.2) +
  labs(x = "log10(virus to host gene ratio)", y = "N viruses", tag="b") +
  scale_x_log10() +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.position = "bottom",
    plot.tag = element_text(face="bold", size=6),
    strip.background = element_rect(fill = "transparent"))

figure_2_supp <- figure_2A_supp + figure_2B_supp + plot_layout(guides = 'collect') & theme(legend.position = "bottom")
ggsave("Supplementary_figure2.pdf", figure_2_supp, width = 24/2.54, height = 10/2.54)
ggsave("Supplementary_figure2.png", figure_2_supp, width = 15.92/2.54, height = 6.63/2.54)

meta_working$ncvssample <- factor(meta_working$ncvssample, levels = c("NCs", "SAMPLES"))

figure_3A_supp <- ggplot(meta_working, aes(x=ncvssample, y=contigs_1000)) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=1, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(alpha=0, outliers = FALSE) +
  facet_grid(. ~ cohort, labeller = labeller(
    cohort = c(
      "garmaeva" = "Garmaeva *et al.* <br>",
      "liang" = "Liang *et al.* <br>",
      "maqsood" = "Maqsood *et al.* <br>",
      "shah" = "Shah *et al.* <br>"))) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  scale_y_log10() +
  labs(y = "Number of assembled contigs \n longer than 1kb (log10)", fill = "Timepoints", tag="a") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.position = "bottom",
    plot.tag = element_text(face="bold", size=6),
    strip.background = element_rect(fill = "transparent")) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 3)
    )
  )

stat.test3a <- meta_working[meta_working$cohort!="garmaeva",] %>%
  group_by(cohort) %>%
  t_test(contigs_1000 ~ ncvssample) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test3a <- stat.test3a %>% add_xy_position(x = "contigs_1000")
stat.test3a[4,] <- stat.test3a[3,]
stat.test3a[4,"cohort"] <- "garmaeva"
stat.test3a$y.position <- log10(stat.test3a$y.position)
stat.test3a[4,"y.position"] <- log10(max(meta_working[meta_working$cohort=="garmaeva",]$contigs_1000 + 50000) )

stat.test3a$xmin <- 1
stat.test3a$xmax <- 2

stat.test3a$p.signif <- c("ns", "ns", "ns", "NA")
figure_3A_supp <- figure_3A_supp + stat_pvalue_manual(stat.test3a, tip.length = 0.02, size=2.5, label = "p.signif")

figure_3B_supp <- ggplot(meta_working, aes(x=ncvssample, y=total_viruses_discovered)) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=1, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(alpha=0, outliers = FALSE) +
  facet_grid(. ~ cohort, labeller = labeller(
    cohort = c(
      "garmaeva" = "Garmaeva *et al.* <br>",
      "liang" = "Liang *et al.* <br>",
      "maqsood" = "Maqsood *et al.* <br>",
      "shah" = "Shah *et al.* <br>"))) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  scale_y_log10() +
  labs(y = "Number of the discovered putative \n viral sequences (log10)", fill = "Timepoints", tag="b") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.position = "bottom",
    plot.tag = element_text(face="bold", size=6),
    strip.background = element_rect(fill = "transparent")) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 3)
    )
  )

stat.test3b <- meta_working[meta_working$cohort!="garmaeva",] %>%
  group_by(cohort) %>%
  t_test(total_viruses_discovered ~ ncvssample) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test3b <- stat.test3b %>% add_xy_position(x = "total_viruses_discovered")
stat.test3b[4,] <- stat.test3b[3,]
stat.test3b[4,"cohort"] <- "garmaeva"
stat.test3b$y.position <- log10(stat.test3b$y.position)
stat.test3b[4,"y.position"] <- log10(max(meta_working[meta_working$cohort=="garmaeva",]$total_viruses_discovered + 20000) )

stat.test3b$xmin <- 1
stat.test3b$xmax <- 2

stat.test3b$p.signif <- c("***", "ns", "ns", "NA")
figure_3B_supp <- figure_3B_supp + stat_pvalue_manual(stat.test3b, tip.length = 0.02, size=2.5, label = "p.signif")


figure_3C_supp <- ggplot(meta_working, aes(x=ncvssample, y=CHM_LU_richness_discovered_ratio)) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=1, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(alpha=0, outliers = FALSE) +
  facet_grid(. ~ cohort, labeller = labeller(
    cohort = c(
      "garmaeva" = "Garmaeva *et al.* <br>",
      "liang" = "Liang *et al.* <br>",
      "maqsood" = "Maqsood *et al.* <br>",
      "shah" = "Shah *et al.* <br>"))) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  labs(y = "Ratio of the discovered vOTUs \n with sufficient quality to all discovered vOTUs", fill = "Timepoints", tag="c") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.position = "bottom",
    plot.tag = element_text(face="bold", size=6),
    strip.background = element_rect(fill = "transparent")) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 3)
    )
  )

stat.test3c <- meta_working[meta_working$cohort!="garmaeva",] %>%
  group_by(cohort) %>%
  t_test(CHM_LU_richness_discovered_ratio ~ ncvssample) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test3c <- stat.test3c %>% add_xy_position(x = "CHM_LU_richness_discovered_ratio")
stat.test3c[4,] <- stat.test3c[3,]
stat.test3c[4,"cohort"] <- "garmaeva"
stat.test3c[4,"y.position"] <- max(meta_working[meta_working$cohort=="garmaeva",]$CHM_LU_richness_discovered_ratio + 0.2)

stat.test3c$xmin <- 1
stat.test3c$xmax <- 2

stat.test3c$p.signif <- c("ns", "ns", "ns", "NA")
figure_3C_supp <- figure_3C_supp + stat_pvalue_manual(stat.test3c, tip.length = 0.02, size=2.5, label = "p.signif")


figure_3D_supp <- ggplot(meta_all_with_qc_curated[meta_all_with_qc_curated$clean_reads_comb > 0, ], aes(x=ncvssample, y=diversity)) +
  geom_jitter(width = 0.3, aes(fill = timepoint_type), size=1, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(alpha=0, outliers = FALSE) +
  facet_grid(. ~ cohort, labeller = labeller(
    cohort = c(
      "garmaeva" = "Garmaeva *et al.* <br>",
      "liang" = "Liang *et al.* <br>",
      "maqsood" = "Maqsood *et al.* <br>",
      "shah" = "Shah *et al.* <br>"))) +
  scale_fill_manual(values = c("#fcbf49", "#f77f00", "#d62828", "#4C6E7F")) +
  labs(y = "Shannon diversity", fill = "Timepoints", tag="d") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.position = "bottom",
    plot.tag = element_text(face="bold", size=6),
    strip.background = element_rect(fill = "transparent")) +
  guides(
    fill = guide_legend(
      override.aes = list(size = 3)
    )
  )

stat.test3d <- meta_working[meta_working$cohort!="garmaeva",] %>%
  group_by(cohort) %>%
  t_test(diversity ~ ncvssample) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test3d <- stat.test3d %>% add_xy_position(x = "diversity")
stat.test3d[4,] <- stat.test3d[3,]
stat.test3d[4,"cohort"] <- "garmaeva"
stat.test3d[4,"y.position"] <- max(meta_working[meta_working$cohort=="garmaeva",]$diversity + 0.2)

stat.test3d$xmin <- 1
stat.test3d$xmax <- 2

stat.test3d$p.signif <- c("*", "ns", "*", "NA")
figure_3D_supp <- figure_3D_supp + stat_pvalue_manual(stat.test3d, tip.length = 0.02, size=2.5, label = "p.signif")

figure_3_supp <- ((figure_3A_supp + figure_3B_supp) / (figure_3C_supp + figure_3D_supp)) + plot_layout(guides = 'collect') & theme(legend.position = "bottom")

# Save the combined plot as a PDF
ggsave("Supplementary_figure3.pdf", figure_3_supp, width = 21/2.54, height = 21/2.54)
ggsave("Supplementary_figure3.png", figure_3_supp, width = 15.92/2.54, height = 15.92/2.54)

#######################################################################################################################################

## Additional analysis: Venn diagram
#######################################################################################################################################
RPKM_NC <- RPKM[, colnames(RPKM) %in% negative_controls]
RPKM_SMPLS <- RPKM[, !(colnames(RPKM) %in% negative_controls)]
RPKM_NC <- RPKM_NC[rowSums(RPKM_NC) > 0, ]
RPKM_SMPLS <- RPKM_SMPLS[rowSums(RPKM_SMPLS) > 0, ]

data <- list(
  vOTUs_NC = row.names(RPKM_NC),
  vOTUs_SMPLS = row.names(RPKM_SMPLS)
)

pdf('venn_vOTUs_NCsandSAMPLES_adjusted.pdf', width=10/2.54, height=10/2.54)
venn.plot <- venn.diagram(
  x = data,
  category.names = c("Detected \n in NCs", "Detected \n in samples"),
  filename = NULL,
  output = TRUE,
  col = "transparent",
  fill = c("skyblue", "pink"),
  alpha = 0.5,
  cex = 1,
  cat.cex = 1,
  cat.fontface = "bold",
  main = "vOTUs detected in the study"
)
grid::grid.draw(venn.plot)
dev.off()

#######################################################################################################################################

## Additional analysis: differential abundance of the vOTUs between samples and NCs (IN CONSTRUCTION)
#######################################################################################################################################
vOTUs_samples_NCs <- intersect(data$vOTUs_NC, data$vOTUs_SMPLS)

sample_garmaeva <- meta_working$Sample_name[meta_working$Type != "Neg_ctrl" & meta_working$cohort == "garmaeva" & meta_working$Sample_name %in% colnames(RPKM)]
sample_liang <- meta_working$Sample_name[meta_working$Type != "Neg_ctrl" & meta_working$cohort == "liang" & meta_working$Sample_name %in% colnames(RPKM)]
sample_maqsood <- meta_working$Sample_name[meta_working$Type != "Neg_ctrl" & meta_working$cohort == "maqsood" & meta_working$Sample_name %in% colnames(RPKM)]
sample_shah <- meta_working$Sample_name[meta_working$Type != "Neg_ctrl" & meta_working$cohort == "shah" & meta_working$Sample_name %in% colnames(RPKM)]

RPKM_shared <- RPKM[row.names(RPKM) %in% vOTUs_samples_NCs, ]
RPKM_shared <- RPKM_shared[rowSums(RPKM_shared) > 0, colSums(RPKM_shared) > 0]
RPKM_shared_count <- RPKM_shared
RPKM_shared_count[RPKM_shared_count > 0] <- 1

RPKM_shared_count$rowsumsnc_g <- RPKM_shared_count$LN_7C08_VL_405
RPKM_shared_count$rowsumssamples_g <- rowSums(RPKM_shared_count[, colnames(RPKM_shared_count) %in% sample_garmaeva])

RPKM_shared_count$rowsumsnc_l <- rowSums(RPKM_shared_count[, colnames(RPKM_shared_count) %in% negative_controls_liang])
RPKM_shared_count$rowsumssamples_l <- rowSums(RPKM_shared_count[, colnames(RPKM_shared_count) %in% sample_liang])

RPKM_shared_count$rowsumsnc_m <- rowSums(RPKM_shared_count[, colnames(RPKM_shared_count) %in% negative_controls_maqsood])
RPKM_shared_count$rowsumssamples_m <- rowSums(RPKM_shared_count[, colnames(RPKM_shared_count) %in% sample_maqsood])

RPKM_shared_count$rowsumsnc_s <- rowSums(RPKM_shared_count[, colnames(RPKM_shared_count) %in% negative_controls_shah])
RPKM_shared_count$rowsumssamples_s <- rowSums(RPKM_shared_count[, colnames(RPKM_shared_count) %in% sample_shah])

vOTUs_interest_garmaeva <- row.names(RPKM_shared_count)[RPKM_shared_count$rowsumsnc_g > 0 & RPKM_shared_count$rowsumssamples_g > 1]
vOTUs_interest_liang <- row.names(RPKM_shared_count)[RPKM_shared_count$rowsumsnc_l > 1 & RPKM_shared_count$rowsumssamples_l > 1]
vOTUs_interest_maqsood <- row.names(RPKM_shared_count)[RPKM_shared_count$rowsumsnc_m > 1 & RPKM_shared_count$rowsumssamples_m > 1]
vOTUs_interest_shah <- row.names(RPKM_shared_count)[RPKM_shared_count$rowsumsnc_s > 1 & RPKM_shared_count$rowsumssamples_s > 1]

result_diff_ab <- RPKM_shared_count[, grep("^rowsums", names(RPKM_shared_count))]
result_diff_ab$vOTU <- row.names(result_diff_ab)
row.names(result_diff_ab) <- NULL

result_diff_ab_g <- result_diff_ab[result_diff_ab$vOTU %in% vOTUs_interest_garmaeva, c("vOTU", "rowsumssamples_g", "rowsumsnc_g")]
result_diff_ab_l <- result_diff_ab[result_diff_ab$vOTU %in% vOTUs_interest_liang, c("vOTU", "rowsumssamples_l", "rowsumsnc_l")]
result_diff_ab_m <- result_diff_ab[result_diff_ab$vOTU %in% vOTUs_interest_maqsood, c("vOTU", "rowsumssamples_m", "rowsumsnc_m")]
result_diff_ab_s <- result_diff_ab[result_diff_ab$vOTU %in% vOTUs_interest_shah, c("vOTU", "rowsumssamples_s", "rowsumsnc_s")]

RPKM_INT$Sample_name <- row.names(RPKM_INT)
row.names(RPKM_INT) <- NULL
RPKM_INT <-  merge(RPKM_INT, meta_working[c("Sample_name", "ncvssample", "cohort", "nc_subject_group")], by="Sample_name", all.x=T)
RPKM_INT$ncvssample <- factor(RPKM_INT$ncvssample, levels = c("SAMPLES", "NCs"))
RPKM_INT <- RPKM_INT[, colnames(RPKM_INT) %in% result_diff_ab_m$vOTU | colnames(RPKM_INT) %in% result_diff_ab_l$vOTU |
                            colnames(RPKM_INT) %in% result_diff_ab_s$vOTU | colnames(RPKM_INT) %in% c("Sample_name", "ncvssample", "cohort", "nc_subject_group")]


wilcox.test(value ~ ncvssample, data=RPKM_INT[RPKM_INT$vOTU == "maqsood_C0132iv_N38_L10939_K9.5_E0_P0_F0" & RPKM_INT$cohort == "maqsood", ])$p.value
summary(lmer(garmaeva_LN_4E07_VL_260_N107_L4930_K24.8_E1_P0_F0 ~ ncvssample + (1|nc_subject_group), REML = F, data = RPKM_INT[RPKM_INT$cohort == "liang" , ]))
coef(summary(lmer(liang_SRR8653029_N31_L12224_K2.6_E0_P0_F0 ~ ncvssample + (1|nc_subject_group), REML = F, data = RPKM_INT[RPKM_INT$cohort == "maqsood" , ])))["ncvssampleNCs", "Pr(>|t|)"]


result_diff_ab_m$pval <- sapply(result_diff_ab_m$vOTU, function(current_vOTU) {
  subset_data <- RPKM_INT[RPKM_INT$cohort == "maqsood", ]
  if (nrow(subset_data) > 1) {
    formula <- as.formula(paste(current_vOTU, "~ ncvssample + (1|nc_subject_group)"))
    model <- lmer(formula, REML = FALSE, data = subset_data)
    p_value <- coef(summary(model))["ncvssampleNCs", "Pr(>|t|)"]
    return(p_value)
  } else {
    return(NA) # Return NA if there's not enough data
  }
})


result_diff_ab_l$pval <- sapply(result_diff_ab_l$vOTU, function(current_vOTU) {
  subset_data <- RPKM_INT[RPKM_INT$cohort == "liang", ]
  if (nrow(subset_data) > 1) {
    formula <- as.formula(paste(current_vOTU, "~ ncvssample + (1|nc_subject_group)"))
    model <- lmer(formula, REML = FALSE, data = subset_data)
    p_value <- coef(summary(model))["ncvssampleNCs", "Pr(>|t|)"]
    return(p_value)
  } else {
    return(NA) # Return NA if there's not enough data
  }
})

result_diff_ab_s$pval <- sapply(result_diff_ab_s$vOTU, function(current_vOTU) {
  subset_data <- RPKM_INT[RPKM_INT$cohort == "shah", ]
  if (nrow(subset_data) > 1) {
    formula <- as.formula(paste(current_vOTU, "~ ncvssample + (1|nc_subject_group)"))
    model <- lmer(formula, REML = FALSE, data = subset_data)
    p_value <- coef(summary(model))["ncvssampleNCs", "Pr(>|t|)"]
    return(p_value)
  } else {
    return(NA) # Return NA if there's not enough data
  }
})


RPKM_t <- as.data.frame(t(RPKM))
RPKM_t <- RPKM_t[, colnames(RPKM_t) %in% result_diff_ab_m$vOTU | colnames(RPKM_t) %in% result_diff_ab_l$vOTU |
                   colnames(RPKM_t) %in% result_diff_ab_s$vOTU]
min_nonzero <- min(RPKM_t[RPKM_t > 0], na.rm = TRUE)
RPKM_t$Sample_name <- row.names(RPKM_t)
row.names(RPKM_t) <- NULL
RPKM_t <-  merge(RPKM_t, meta_working[c("Sample_name", "ncvssample", "cohort", "nc_subject_group")], by="Sample_name", all.x=T)
RPKM_t$ncvssample <- factor(RPKM_t$ncvssample, levels = c("SAMPLES", "NCs"))



result_diff_ab_m$log2fold_change <- sapply(result_diff_ab_m$vOTU, function(current_vOTU) {
  subset_data <- RPKM_t[RPKM_t$cohort == "maqsood", ]
  if (nrow(subset_data) > 1) {
    values_NCs <- subset_data[[current_vOTU]][subset_data$ncvssample == "NCs"]
    values_SAMPLES <- subset_data[[current_vOTU]][subset_data$ncvssample == "SAMPLES"]
    test_result <- log2(mean(values_NCs)/mean(values_SAMPLES))
    return(test_result)
  } else {
    return(NA) # Return NA if there's not enough data
  }
})

result_diff_ab_l$log2fold_change <- sapply(result_diff_ab_l$vOTU, function(current_vOTU) {
  subset_data <- RPKM_t[RPKM_t$cohort == "liang", ]
  if (nrow(subset_data) > 1) {
    values_NCs <- subset_data[[current_vOTU]][subset_data$ncvssample == "NCs"]
    values_SAMPLES <- subset_data[[current_vOTU]][subset_data$ncvssample == "SAMPLES"]
    test_result <- log2(mean(values_NCs)/mean(values_SAMPLES))
    return(test_result)
  } else {
    return(NA) # Return NA if there's not enough data
  }
})

result_diff_ab_s$log2fold_change <- sapply(result_diff_ab_s$vOTU, function(current_vOTU) {
  subset_data <- RPKM_t[RPKM_t$cohort == "shah", ]
  if (nrow(subset_data) > 1) {
    values_NCs <- subset_data[[current_vOTU]][subset_data$ncvssample == "NCs"]
    values_SAMPLES <- subset_data[[current_vOTU]][subset_data$ncvssample == "SAMPLES"]
    test_result <- log2(mean(values_NCs)/mean(values_SAMPLES))
    return(test_result)
  } else {
    return(NA) # Return NA if there's not enough data
  }
})

write.table(result_diff_ab_m, "../result_diff_ab_m_upd.tsv", sep='\t', row.names=F, col.names=T, quote=F)
write.table(result_diff_ab_l, "../result_diff_ab_l_upd.tsv", sep='\t', row.names=F, col.names=T, quote=F)
write.table(result_diff_ab_s, "../result_diff_ab_s_upd.tsv", sep='\t', row.names=F, col.names=T, quote=F)

#######################################################################################################################################

## Additional analysis: vOTU and strain sharedness correlation
#######################################################################################################################################
strains_df$name1 <- gsub("_filtered.sorted.bam", "", strains_df$name1)
strains_df$name2 <- gsub("_filtered.sorted.bam", "", strains_df$name2)
colnames(strains_df)[colnames(strains_df) %in% c("name1", "name2")] <- c("Sample1", "Sample2")
strains_df <- merge(strains_df, meta_working[c("Sample1", "ncvssample", "cohort")], by="Sample1", all.x=T)
strains_df <- merge(strains_df, meta_working[c("Sample2", "ncvssample", "cohort")], by="Sample2", all.x=T, suffixes=c("_sample1", "_sample2"))
strains_df <- strains_df[strains_df$cohort_sample1 == strains_df$cohort_sample2 & strains_df$ncvssample_sample1 != strains_df$ncvssample_sample2, ]
strains_df <- strains_df %>%
  mutate(SAMPLE = ifelse(ncvssample_sample1 == "SAMPLES", Sample1, Sample2),
         NC = ifelse(ncvssample_sample1 == "NCs", Sample1, Sample2))

strains_df_popani_99999 <- strains_df[strains_df$popANI > 0.99999, ]
shared_list <- unique(strains_df_popani_99999[c("SAMPLE", "scaffold")])
presence_strain <- as.data.frame(table(shared_list$SAMPLE))
colnames(presence_strain) <- c("Sample_name", "N_strain_shared")

rpkm_sums <- merge(as.data.frame(colSums(RPKM)), as.data.frame(colSums(RPKM > 0)), by="row.names")
colnames(rpkm_sums) <-c("Sample_name", "total_abundance", "total_presence")
strains_vOTU_shared <- merge(table_for_plot_cor_combine, rpkm_sums, all.x=T, by="Sample_name")
strains_vOTU_shared <- merge(strains_vOTU_shared, presence_strain, all.x=T, by="Sample_name")
strains_vOTU_shared$N_strain_shared[is.na(strains_vOTU_shared$N_strain_shared)] <- 0
strains_vOTU_shared$abundance_strain_shared <- 0

for (i in shared_list$SAMPLE) {
  strains_vOTU_shared$abundance_strain_shared[strains_vOTU_shared$Sample_name == i] <- 
    sum(RPKM[row.names(RPKM) %in% shared_list$scaffold[shared_list$SAMPLE == i], i], na.rm = TRUE)
}

strains_vOTU_shared$presence_strain_shared_perc <- (strains_vOTU_shared$N_strain_shared / strains_vOTU_shared$total_presence)*100
strains_vOTU_shared$abundance_strain_shared_perc <- (strains_vOTU_shared$abundance_strain_shared / strains_vOTU_shared$total_abundance)*100

# Strain sharedness to different cohorts

strains_df_diff <- strains_df_ini[!is.na(strains_df_ini$popANI), ]
strains_df_diff$name1 <- gsub("_filtered.sorted.bam", "", strains_df_diff$name1)
strains_df_diff$name2 <- gsub("_filtered.sorted.bam", "", strains_df_diff$name2)
colnames(strains_df_diff)[colnames(strains_df_diff) %in% c("name1", "name2")] <- c("Sample1", "Sample2")
strains_df_diff <- merge(strains_df_diff, meta_working[c("Sample1", "ncvssample", "cohort")], by="Sample1", all.x=T)
strains_df_diff <- merge(strains_df_diff, meta_working[c("Sample2", "ncvssample", "cohort")], by="Sample2", all.x=T, suffixes=c("_sample1", "_sample2"))
strains_df_diff <- strains_df_diff[strains_df_diff$cohort_sample1 != strains_df_diff$cohort_sample2 &
                                     strains_df_diff$ncvssample_sample1 != strains_df_diff$ncvssample_sample2, ]

strains_df_diff <- strains_df_diff %>%
  mutate(SAMPLE = ifelse(ncvssample_sample1 == "SAMPLES", Sample1, Sample2),
         NC = ifelse(ncvssample_sample1 == "NCs", Sample1, Sample2),
         cohort_sample = ifelse(ncvssample_sample1 == "SAMPLES", as.character(cohort_sample1), as.character(cohort_sample2)))

strains_df_popani_99999_diff <- strains_df_diff[strains_df_diff$popANI > 0.99999, ]
shared_list_diff <- unique(strains_df_popani_99999_diff[c("SAMPLE", "scaffold")])
presence_strain_diff <- as.data.frame(table(shared_list_diff$SAMPLE))
colnames(presence_strain_diff) <- c("Sample_name", "N_strain_shared_diff_study")

strains_vOTU_shared <- merge(strains_vOTU_shared, presence_strain_diff, by="Sample_name", all.x=T)
strains_vOTU_shared$N_strain_shared_diff_study[is.na(strains_vOTU_shared$N_strain_shared_diff_study)] <- 0
strains_vOTU_shared$presence_strain_shared_perc_diff_study <- (strains_vOTU_shared$N_strain_shared_diff_study / strains_vOTU_shared$total_presence)*100

# for plot correlation

spearman_result_4a <- psych::corr.test(strains_vOTU_shared[c("presence_strain_shared_perc", "presence_strain_shared_perc_diff_study")], method = "spearman")
print(spearman_result_4a$r)  # 0.3677397
print(spearman_result_4a$p)  # 1.92456e-41

spearman_result_garmaeva_4b <- psych::corr.test(strains_vOTU_shared[strains_vOTU_shared$cohort == "garmaeva", c("presence_strain_shared_perc", "different_cohort_NC_presence")], method = "spearman")
print(spearman_result_garmaeva_4b$r)
print(spearman_result_garmaeva_4b$p)

spearman_result_maqsood_4b <- psych::corr.test(strains_vOTU_shared[strains_vOTU_shared$cohort == "maqsood", c("presence_strain_shared_perc", "different_cohort_NC_presence")], method = "spearman")
print(spearman_result_maqsood_4b$r)
print(spearman_result_maqsood_4b$p)

spearman_result_liang_4b <- psych::corr.test(strains_vOTU_shared[strains_vOTU_shared$cohort == "liang", c("presence_strain_shared_perc", "different_cohort_NC_presence")], method = "spearman")
print(spearman_result_liang_4b$r)
print(spearman_result_liang_4b$p)

spearman_result_shah_4b <- psych::corr.test(strains_vOTU_shared[strains_vOTU_shared$cohort == "shah", c("presence_strain_shared_perc", "different_cohort_NC_presence")], method = "spearman")
print(spearman_result_shah_4b$r)
print(spearman_result_shah_4b$p)

correlations_presence_strain_vs_diffvotu <- strains_vOTU_shared %>%
  group_by(cohort) %>%
  summarize(cor = cor(presence_strain_shared_perc, different_cohort_NC_presence, method = "spearman"))

correlations_presence_strain_vs_diffstrain <- strains_vOTU_shared %>%
  summarize(cor = cor(presence_strain_shared_perc, presence_strain_shared_perc_diff_study, method = "spearman"))
#######################################################################################################################################

## Figure 4
#######################################################################################################################################

Figure4a <- ggplot(strains_vOTU_shared, aes(x=presence_strain_shared_perc, y=presence_strain_shared_perc_diff_study)) +
  geom_point(size = 0.7, color="#2E236C", alpha=0.75) +  # Adjusted point size and added transparency
  geom_smooth(method="lm", color="#2E236C", fill="#C8ACD6", se=TRUE, linewidth=0.5) +  # Use linewidth instead of size
  labs(x = "% shared strains with NCs from same study", y = "% shared strains with NCs from different studies", tag="a") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = ggtext::element_markdown(face = "bold")
  ) +
  geom_text(data = correlations_presence_strain_vs_diffstrain, aes(x = Inf, y = Inf, label = paste("rho = ", round(cor, 2))),
            hjust = 1.1, vjust = 1.5, size = 3)


Figure4b <- ggplot(strains_vOTU_shared, aes(x=presence_strain_shared_perc, y=different_cohort_NC_presence)) +
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
  labs(x = "% shared strains with NCs from same study", y = "% shared vOTUs with NCs from different studies", tag="b") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=7),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = ggtext::element_markdown(face = "bold")
  ) +
  geom_text(data = correlations_presence_strain_vs_diffvotu, aes(x = Inf, y = Inf, label = paste("rho = ", round(cor, 2))),
            hjust = 1.1, vjust = 1.5, size = 3)

Figure4 <- Figure4a + Figure4b

ggsave("Figure4.pdf", Figure4, width = 21/2.54, height = 12/2.54)
#######################################################################################################################################

## Part of the figure 3d & saving the df for the figure 3abe
#######################################################################################################################################

summary(lmer(N_strain_shared ~ Type + (1|nc_subject_group), data=strains_vOTU_shared[strains_vOTU_shared$cohort == "garmaeva", ]))
summary(lm(N_strain_shared ~ Type, data=strains_vOTU_shared[strains_vOTU_shared$cohort == "maqsood", ]))


df_figure3abe <- merge(strains_vOTU_shared[strains_vOTU_shared$cohort %in% c("garmaeva", "maqsood"), c("Sample_name", "Type", "cohort", "same_cohort_NC_presence", "N_strain_shared")],
                       result_present_in_both_mi[c("Sample_name", "same_cohort_NC")], by="Sample_name")

write.table(df_figure3abe, "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/df_figure3abe.tsv", sep='\t', row.names=F, col.names=T, quote=F)
#######################################################################################################################################

## Metadata tuning for public repository
#######################################################################################################################################
meta_all_with_qc_curated <- merge(meta_all_with_qc_curated, table_for_plot_cor_combine[c("Sample_name", "same_cohort_NC_presence", "different_cohort_NC_presence",
                                                                                         "same_cohort_NC_abundance", "different_cohort_NC_abundance")], by="Sample_name", all.x = T)

meta_all_with_qc_curated <- merge(meta_all_with_qc_curated, result_present_in_both[c("Sample_name", "same_cohort_NC", "different_cohort_NC")], by="Sample_name", all.x = T)

meta_all_with_qc_curated <- merge(meta_all_with_qc_curated, strains_vOTU_shared[c("Sample_name", "presence_strain_shared_perc", "abundance_strain_shared_perc", 
                                                                                  "N_strain_shared", "N_strain_shared_diff_study", "presence_strain_shared_perc_diff_study")], by="Sample_name", all.x = T)

meta_all_with_qc_curated_clean <- meta_all_with_qc_curated %>%
  mutate(in_reads_se = NULL,
         clean_reads_se = NULL,
         to_1kb_contigs = NULL,
         to_all_vir = NULL,
         to_ext_vir = NULL,
         to_ext_prun_vir = NULL,
         DeepVirFinder = NULL,
         geNomad = NULL,
         VIBRANT = NULL,
         VirSorter2 = NULL,
         Subject_ID = NULL,
         Sample_ID = NULL,
         Subcohort = NULL,
         COMMENTS = NULL,
         incl_longitudinal = ifelse(Sample_name %in% combined_table_for_plot$Sample_name, "yes", "no"),
         richness = ifelse(clean_reads_comb == 0, NA, richness),
         diversity = ifelse(clean_reads_comb == 0, NA, diversity)) %>%
  rename(status = ncvssample,
         Study = cohort,
         perc_reads_mapped_to_all_contigs = to_all_contigs,
         clean_reads_mapped_to_vOTUs = reads_mapped,
         Subject_ID = nc_subject_group,
         perc_shared_own_NC = same_cohort_NC_presence,
         perc_shared_other_NC = different_cohort_NC_presence,
         N_shared_own_NC = same_cohort_NC,
         N_shared_other_NC = different_cohort_NC,
         rna_dna_liang = rna_dna,
         timepoint_custom_category = timepoint_type,       
         perc_abundance_shared_own_NC = same_cohort_NC_abundance,
         perc_abundance_shared_other_NC = different_cohort_NC_abundance,
         N_strains_shared_own_NC = N_strain_shared,
         N_strains_shared_other_NC = N_strain_shared_diff_study,
         perc_strains_shared_own_NC = presence_strain_shared_perc,
         perc_strains_shared_other_NC = presence_strain_shared_perc_diff_study,
         perc_abundance_strains_shared_own_NC = abundance_strain_shared_perc)

meta_all_with_qc_curated_clean <- meta_all_with_qc_curated_clean[c("Sample_name", "Subject_ID", "FAM_ID", "Type", "status", "Timepoint", "timepoint_custom_category",
                                                                   "Study", "input_files", "in_reads_comb", "clean_reads_comb", "dedup_efficiency", "contigs_total",
                                                                   "contigs_1000", "N50", "total_viruses_discovered", "perc_reads_mapped_to_all_contigs",
                                                                   "clean_reads_mapped_to_vOTUs", "richness", "diversity", "N_shared_own_NC", "N_shared_other_NC",
                                                                   "perc_shared_own_NC", "perc_shared_other_NC", "perc_abundance_shared_own_NC", "perc_abundance_shared_other_NC",
                                                                   "N_strains_shared_own_NC", "N_strains_shared_other_NC", "perc_strains_shared_own_NC", 
                                                                   "perc_strains_shared_other_NC", "perc_abundance_strains_shared_own_NC", "rna_dna_liang", 
                                                                   "incl_longitudinal")]




write.table(meta_all_with_qc_curated_clean, "/scratch/hb-llnext/VLP_public_data/nc_project/for_upload/Sample_metadata_v2.tsv", sep='\t', row.names=F, col.names=T, quote=F)
write.table(meta_all_with_qc_curated_clean, "/scratch/hb-llnext/VLP_public_data/nc_project/Sample_metadata_upd.tsv", sep='\t', row.names=F, col.names=T, quote=F)
#######################################################################################################################################