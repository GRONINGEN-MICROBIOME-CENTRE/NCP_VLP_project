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
#######################################################################################################################################

## Set the working directory
#######################################################################################################################################
setwd("/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/plots")
#######################################################################################################################################

## Load the  data
#######################################################################################################################################
extended_tof <- read.delim('../../VIR_DB/virus_contigs/MERGED_Extended_TOF_NCP')
row.names(extended_tof) <- extended_tof$New_CID
dim(extended_tof)  # 971583     41
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
dim(RPKM_initial)  # 306275   1301

RPKM <- RPKM_initial[row.names(RPKM_initial) %in% RPKM_contigs_keep, ]
dim(RPKM)  # 193970   1301
RPKM <- RPKM[rowSums(RPKM) > 0, colSums(RPKM) > 0]
dim(RPKM)  # 193970   1291

RPKM_count <- RPKM
RPKM_count[RPKM_count > 0] <- 1

# RPKM_cleaned_98_cov75 <- read.delim("../RPKM_counts_VLP_NCP_fin_98_cov75.txt")
# RPKM_cleaned_99_cov75 <- read.delim("../RPKM_counts_VLP_NCP_fin_99_cov75.txt")
# RPKM_cleaned_95_cov75 <- read.delim("../RPKM_counts_VLP_NCP_fin_cleaned_95_cov75.txt")

meta_all_with_qc_curated <- as.data.frame(read_tsv('../../metadata_with_qc_NCPv2.tsv'))
dim(meta_all_with_qc_curated)  # 1376   28
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
         timepoint_type = ifelse(Type == "Mother", "Mother", timepoint_type)
  )

meta_all_with_qc_curated$nc_subject_group <- as.factor(meta_all_with_qc_curated$nc_subject_group)
meta_all_with_qc_curated$cohort <- as.factor(meta_all_with_qc_curated$cohort)
meta_all_with_qc_curated$ncvssample <- as.factor(meta_all_with_qc_curated$ncvssample)
meta_all_with_qc_curated$ncvssample <- factor(meta_all_with_qc_curated$ncvssample, levels = c("SAMPLES", "NCs"))

meta_working <- meta_all_with_qc_curated[meta_all_with_qc_curated$Sample_name %in% colnames(RPKM), ]
meta_working$ncvssample <- factor(meta_working$ncvssample, levels = c("SAMPLES", "NCs"))

host_prediction <- read.csv('../../VIR_DB/host_prediction_w_neg_der95_NCP/results/MERGED_Host_prediction_to_genus_m90_v2.csv')
dim(host_prediction)  # 216679      5 

filtered_host_prediction <- host_prediction %>%
  group_by(Virus) %>%
  slice_max(order_by = Confidence.score, with_ties = FALSE) %>%
  ungroup()
filtered_host_prediction <- as.data.frame(filtered_host_prediction)
dim(filtered_host_prediction)  # 188122      5

strains_df_ini <- read.delim('../inStrain_70_vOTUr_pairwise_popANI.txt')
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

# Applying the function to create the new columns
split_results <- t(sapply(extended_tof$New_CID, get_study_and_sample_of_origin))

# Adding the new columns to the dataframe
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
summarized_df2_frac$Sample_name <- NULL

summarized_df2_frac_stats <- melt(summarized_df2_frac, id.vars = c("ncvssample", "cohort", "nc_subject_group"))

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


summarized_df2_frac <- summarized_df2_frac %>%
  group_by(ncvssample, cohort) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

summarized_df2_frac <- as.data.frame(summarized_df2_frac)

## Redo for long format
summarized_df2_frac_melt <- melt(summarized_df2_frac, id.vars = c("ncvssample", "cohort"))
summarized_df2_frac_melt$variable <- as.factor(summarized_df2_frac_melt$variable)
summarized_df2_frac_melt$variable <- factor(summarized_df2_frac_melt$variable, levels = c("RNA", "ssDNA", "dsDNA", "Unclassified"))
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

summarized_df3_frac_stats <- summarized_df3_frac

summarized_df3_frac_stats$Prokaryotes_log_trnsfrmd <- summarized_df3_frac_stats$Prokaryotes + (min(summarized_df3_frac_stats$Prokaryotes[summarized_df3_frac_stats$Prokaryotes > 0])/2)
summarized_df3_frac_stats$Eukaryotes_log_trnsfrmd <- summarized_df3_frac_stats$Eukaryotes + (min(summarized_df3_frac_stats$Eukaryotes[summarized_df3_frac_stats$Eukaryotes > 0])/2)

summarized_df3_frac_stats$Prokaryotes_log_trnsfrmd <- log(summarized_df3_frac_stats$Prokaryotes_log_trnsfrmd)
summarized_df3_frac_stats$Eukaryotes_log_trnsfrmd <- log(summarized_df3_frac_stats$Eukaryotes_log_trnsfrmd)

prok_vs_euk_combined <- summary(lmer(Prokaryotes_log_trnsfrmd ~ Eukaryotes_log_trnsfrmd + (1|cohort/nc_subject_group), REML = F, data = summarized_df3_frac_stats))
prok_vs_euk_combined$coefficients
prok_vs_euk_liang <- summary(lmer(Prokaryotes_log_trnsfrmd ~ Eukaryotes_log_trnsfrmd + (1|nc_subject_group), REML = F, data = summarized_df3_frac_stats[summarized_df3_frac_stats$cohort == "liang", ]))
prok_vs_euk_liang$coefficients
prok_vs_euk_maqsood <- summary(lmer(Prokaryotes_log_trnsfrmd ~ Eukaryotes_log_trnsfrmd + (1|nc_subject_group), REML = F, data = summarized_df3_frac_stats[summarized_df3_frac_stats$cohort == "maqsood", ]))
prok_vs_euk_maqsood$coefficients
prok_vs_euk_shah <- summary(lmer(Prokaryotes_log_trnsfrmd ~ Eukaryotes_log_trnsfrmd + (1|nc_subject_group), REML = F, data = summarized_df3_frac_stats[summarized_df3_frac_stats$cohort == "shah", ]))
prok_vs_euk_shah$coefficients
p.adjust(c(5.363049e-31, 1.41e-06, 7.386918e-112))

summarized_df3_frac <- summarized_df3_frac %>%
  group_by(ncvssample, cohort) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

summarized_df3_frac <- as.data.frame(summarized_df3_frac)

## Redo for long format
summarized_df3_frac_melt <- melt(summarized_df3_frac, id.vars = c("ncvssample", "cohort"))
summarized_df3_frac_melt$variable <- as.factor(summarized_df3_frac_melt$variable)
summarized_df3_frac_melt$variable <- factor(summarized_df3_frac_melt$variable, levels = c("Prokaryotes", "Eukaryotes", "Unclassified"))
#######################################################################################################################################

## Analysis for the phages host richness
#######################################################################################################################################
BU_viruses <- row.names(extended_tof[extended_tof$host_group %in% c("Prokaryotes", "Unclassified"), ])
RPKM_count_ch4 <- RPKM_count[row.names(RPKM_count) %in% CHM_viruses & row.names(RPKM_count) %in% BU_viruses, ]


row.names(filtered_host_prediction) <- filtered_host_prediction$Virus
filtered_host_prediction$Genus <- sub(".*;g__", "", filtered_host_prediction$Host.genus)
filtered_host_prediction$Genus[filtered_host_prediction$Genus == ""] <- "Unclassified"

RPKM_count_ch4 <- merge(RPKM_count_ch4, filtered_host_prediction[colnames(filtered_host_prediction) %in% c("Genus")], all.x = T, by="row.names")
RPKM_count_ch4$Genus[is.na(RPKM_count_ch4$Genus)] <- "Unclassified"
RPKM_count_ch4$Row.names <- NULL

summarized_df4 <- RPKM_count_ch4 %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

summarized_df4 <- as.data.frame(t(summarized_df4))
colnames(summarized_df4) <- c(t(summarized_df4[1,]))
summarized_df4 <- summarized_df4[row.names(summarized_df4) != "Genus", ]

summarized_df4 <- summarized_df4 %>%
  mutate(across(where(is.character), ~ as.numeric(.)))

summarized_df4 <- summarized_df4[rowSums(summarized_df4) > 0, ]

summarized_df4_frac <- as.data.frame(t(apply(summarized_df4, 1, fraction_per_sample)))

summarized_df4_frac$Sample_name <- row.names(summarized_df4_frac)
summarized_df4_frac <- merge(summarized_df4_frac, meta_working[, colnames(meta_working) %in% c("Sample_name", "ncvssample", "cohort")], all.x = T, by="Sample_name")
summarized_df4_frac$Sample_name <- NULL

summarized_df4_frac <- summarized_df4_frac %>%
  group_by(ncvssample, cohort) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

summarized_df4_frac <- as.data.frame(summarized_df4_frac)

columns_other <- c()
for (col in names(summarized_df4_frac[3:ncol(summarized_df4_frac)])) {
  if (max(summarized_df4_frac[[col]], na.rm = TRUE) < 0.01) {
    columns_other <- c(columns_other, col)
  }
}

summarized_df4_frac <- summarized_df4_frac %>%
  mutate(Other = rowSums(select(., all_of(columns_other)), na.rm = TRUE))

summarized_df4_frac <- summarized_df4_frac %>%
  select(-all_of(columns_other))


## Redo for long format
summarized_df4_frac_melt <- melt(summarized_df4_frac, id.vars = c("ncvssample", "cohort"))
summarized_df4_frac_melt$variable <- as.factor(summarized_df4_frac_melt$variable)

summarized_df4_frac_melt$variable <- factor(summarized_df4_frac_melt$variable, 
                                            levels = levels(summarized_df4_frac_melt$variable)[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 31)])

#######################################################################################################################################

## Patching the Figure 1
#######################################################################################################################################
write.table(meta_working, "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/df_for_figures1ab.tsv", sep='\t', row.names=F, col.names=T, quote=F)
write.table(summarized_df2_frac_melt, "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/df_for_figure1c.tsv", sep='\t', row.names=F, col.names=T, quote=F)
write.table(summarized_df3_frac_melt, "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/df_for_figure1d.tsv", sep='\t', row.names=F, col.names=T, quote=F)
write.table(summarized_df4_frac_melt, "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/df_for_figures1e.tsv", sep='\t', row.names=F, col.names=T, quote=F)

meta_working$ncvssample <- factor(meta_working$ncvssample, levels = c("NCs", "SAMPLES"))
meta_all_with_qc_curated$ncvssample <- factor(meta_all_with_qc_curated$ncvssample, levels = c("NCs", "SAMPLES"))

figure_1A <- ggplot(meta_working, aes(x=ncvssample, y=clean_reads_comb)) +
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
  labs(y = "Number of clean reads (log10)", tag="a", fill = "Timepoints") +
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

figure_1B <- ggplot(meta_all_with_qc_curated[meta_all_with_qc_curated$clean_reads_comb > 0, ], aes(x=ncvssample, y=richness)) +
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
  labs(y = "Richness (log10)", tag="b", fill = "Timepoints") +
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


summarized_df2_frac_melt$ncvssample <- factor(summarized_df2_frac_melt$ncvssample, levels = c("NCs", "SAMPLES"))
figure_1C <- ggplot(summarized_df2_frac_melt, aes(x = ncvssample, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("#a1c181", "#fcca46", "#fe7f2d", "#233d4d")) +
  facet_grid(. ~ cohort, labeller = labeller(
    cohort = c(
      "garmaeva" = "Garmaeva *et al.* <br>",
      "liang" = "Liang *et al.* <br>",
      "maqsood" = "Maqsood *et al.* <br>",
      "shah" = "Shah *et al.* <br>"))) +
  labs(tag="c") +
  theme(
    strip.text = ggtext::element_markdown(size=7),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=13),
    legend.position = "bottom"
  ) +
  labs(x = "", y = "Mean fraction of the viral group richness \n (vOTUs with at least 50% completeness)", fill = "Viral group")


summarized_df3_frac_melt$ncvssample <- factor(summarized_df3_frac_melt$ncvssample, levels = c("NCs", "SAMPLES"))
figure_1D <- ggplot(summarized_df3_frac_melt, aes(x = ncvssample, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("#619b8a", "#f4a261", "#233d4d")) +
  facet_grid(. ~ cohort, labeller = labeller(
    cohort = c(
      "garmaeva" = "Garmaeva *et al.* <br>",
      "liang" = "Liang *et al.* <br>",
      "maqsood" = "Maqsood *et al.* <br>",
      "shah" = "Shah *et al.* <br>"))) +
  labs(tag="d")+
  theme(
    strip.text = ggtext::element_markdown(size=7),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=13),
    legend.position = "bottom"
  ) +
  labs(x = "", y = "Mean fraction of the host group richness \n (vOTUs with at least 50% completeness)", fill = "Host group")

summarized_df4_frac_melt$ncvssample <- factor(summarized_df4_frac_melt$ncvssample, levels = c("NCs", "SAMPLES"))
figure_1E <- ggplot(summarized_df4_frac_melt, aes(x = ncvssample, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("cornflowerblue", "gold2", "turquoise3", "orchid", "steelblue", "burlywood", "purple3", "slategray4","yellow", "blue", "skyblue","darkgoldenrod1",
                               "darksalmon", "green", "pink", "tomato2", "rosybrown3", "snow2","plum2", "aquamarine3",
                               "palegreen3", "seagreen", "aquamarine", "yellow3", "chartreuse", "thistle", "orange2", 
                               "navajowhite", "firebrick", "purple", "violetred1", "red", "darkblue", "black")) +
  facet_grid(. ~ cohort, labeller = labeller(
    cohort = c(
      "garmaeva" = "Garmaeva *et al.* <br>",
      "liang" = "Liang *et al.* <br>",
      "maqsood" = "Maqsood *et al.* <br>",
      "shah" = "Shah *et al.* <br>"))) +
  labs(tag="e") +
  theme(
    strip.text = ggtext::element_markdown(size=7),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7, face = "italic"),
    legend.key.size = unit(0.3, "cm"),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=13),
    legend.position = "bottom"
  ) +
  labs(x = "", y = "Mean fraction of the host genus richness \n (vOTUs with at least 50% completeness)", fill = "Host genus") +
  guides(fill=guide_legend(nrow=5, byrow=F, title.position = 'top', title.hjust = 0.5,label.theme = element_text(face = "italic", size=7)), 
         alpha=guide_legend(nrow=2, byrow=TRUE, title.position = 'top', title.hjust = 0.5))


figure_1AB <- figure_1A + figure_1B + plot_layout(nrow=4, guides = "collect")
# Combine the plots using patchwork
combined_plot <- (figure_1A + figure_1B + plot_layout(nrow=1, guides = "collect") & theme(legend.position = "bottom")) / (figure_1C + figure_1D) / figure_1E

# Save the combined plot as a PDF
ggsave("combined_figure.pdf", combined_plot, width = 21/2.54, height = 29.7/2.54)

meta_working$ncvssample <- factor(meta_working$ncvssample, levels = c("SAMPLES", "NCs"))
meta_all_with_qc_curated$ncvssample <- factor(meta_all_with_qc_curated$ncvssample, levels = c("SAMPLES", "NCs"))

summarized_df2_frac_melt$ncvssample <- factor(summarized_df2_frac_melt$ncvssample, levels = c("SAMPLES", "NCs"))
summarized_df3_frac_melt$ncvssample <- factor(summarized_df3_frac_melt$ncvssample, levels = c("SAMPLES", "NCs"))
summarized_df4_frac_melt$ncvssample <- factor(summarized_df4_frac_melt$ncvssample, levels = c("SAMPLES", "NCs"))

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

# Calculate the ICC
icc <- var_cohort / (var_cohort + var_residual)
icc


model_richness_NCs <- lmer(richness ~ 1 + (1|cohort), REML = FALSE, data = meta_working[meta_working$ncvssample == "NCs", ])
variance_components <- VarCorr(model_richness_NCs)
var_cohort <- as.numeric(variance_components$cohort[1])
var_residual <- attr(variance_components, "sc")^2

# Calculate the ICC
icc <- var_cohort / (var_cohort + var_residual)
icc

model_richness_samples <- lmer(richness ~ 1 + (1|cohort), REML = FALSE, data = meta_working[meta_working$ncvssample == "SAMPLES", ])
variance_components <- VarCorr(model_richness_samples)
var_cohort <- as.numeric(variance_components$cohort[1])
var_residual <- attr(variance_components, "sc")^2

# Calculate the ICC
icc <- var_cohort / (var_cohort + var_residual)
icc



summary(lmer(total_viruses_discovered ~ 1 + (1|cohort), REML = FALSE, data = meta_all_with_qc_curated[meta_all_with_qc_curated$ncvssample == "NCs", ]))

#######################################################################################################################################

## Bray-Curtis boxplots: all samples; log scale
#######################################################################################################################################

bray_dist_matrix_full <- as.matrix(vegdist(t(RPKM), method="bray"))
bray_dist_matrix_full_rev <- 1 - bray_dist_matrix_full

bray_dist_matrix_full_rev <- read.delim('bray_dist_matrix_full_rev.tsv')
#write.table(bray_dist_matrix_full_rev, "bray_dist_matrix_full_rev.tsv", sep='\t', row.names=T, col.names=T, quote=F)

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

# Save the rsult for updating eTOF later
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
         different_cohort_NC = ifelse(cohort == "garmaeva", dummy_non_garmaeva_NC, different_cohort_NC),
         Timepoint_numeric = as.integer(gsub("M", "", Timepoint)),
         Timepoint_numeric = ifelse(grepl("Y2-5", Timepoint), 24, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother" & Timepoint == "Mtrim3", 0, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother" & Timepoint == "M0", 3, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother" & Timepoint == "M1", 4, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother" & Timepoint == "M2", 5, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother" & Timepoint == "M3", 6, Timepoint_numeric))

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
         different_cohort_NC = ifelse(cohort == "garmaeva", dummy_non_garmaeva_NC, different_cohort_NC),
         Timepoint_numeric = as.integer(gsub("M", "", Timepoint)),
         Timepoint_numeric = ifelse(grepl("Y2-5", Timepoint), 24, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother" & Timepoint == "Mtrim3", 0, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother" & Timepoint == "M0", 3, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother" & Timepoint == "M1", 4, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother" & Timepoint == "M2", 5, Timepoint_numeric),
         Timepoint_numeric = ifelse(Type == "Mother" & Timepoint == "M3", 6, Timepoint_numeric))


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

write.table(table_for_plot_cor_combine, "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/df_for_figure2d.tsv", sep='\t', row.names=F, col.names=T, quote=F)

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
combined_table_for_plot_log <- combined_table_for_plot
combined_table_for_plot_log$value <- combined_table_for_plot_log$value + (min_nonzero / 2)

write.table(combined_table_for_plot_log, "/scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R/df_for_figure2e.tsv", sep='\t', row.names=F, col.names=T, quote=F)
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


figure_4A_supp <- ggplot(table_for_plot_cor_combine_melt, aes(x = same_diff, y = value)) +
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
    strip.text = ggtext::element_markdown(size=6),
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 8), 
    axis.text.y = element_text(size = 6),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=6)
  )

stat.test4 <- table_for_plot_cor_combine_melt %>%
  group_by(cohort, pres_abun) %>%
  t_test(value ~ same_diff) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test4 <- stat.test4 %>% 
  add_xy_position(x = "same_diff", dodge = 0.8)

stat.test4$xmin <- 1
stat.test4$xmax <- 2

stat.test4$p.signif <- c("****", "****", "ns", "ns", "ns", "ns", "****", "****")

figure_4A_supp <- figure_4A_supp + 
  stat_pvalue_manual(stat.test4, tip.length = 0.02, size=2.5, label = "p.signif")

figure_4B_supp <- ggplot(table_for_plot_cor_combine_melt_subset_log, aes(x = Type, y = same_cohort_NC_presence)) +
  geom_jitter(width = 0.3, fill = "#2E236C", size = 1.3, shape = 21, stroke = 0.1, color = "white") +
  geom_boxplot(fill = "#C8ACD6", alpha=0.3, outlier.alpha = 0) +
  scale_y_log10() +
  facet_grid(. ~ cohort, scales = "free", labeller = labeller(
    cohort = c(
      "garmaeva" = "Samples Garmaeva *et al.* <br>",
      "maqsood" = "Samples Maqsood *et al.* <br>"
    )
  )) +
  labs(y = "log10(% vOTUs shared with NCs)", tag="b") +
  theme_bw() +
  theme(
    strip.text = ggtext::element_markdown(size=6),
    plot.title = element_text(size = 10),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    strip.background = element_rect(fill = "transparent"),
    plot.tag = element_text(face="bold", size=6)
  )

stat.test5 <- table_for_plot_cor_combine_melt_subset_log %>%
  group_by(cohort) %>%
  t_test(same_cohort_NC_presence ~ Type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test5 <- stat.test5 %>% add_xy_position(x = "same_cohort_NC_presence")
stat.test5$y.position <- log10(stat.test5$y.position)

stat.test5$xmin <- 1
stat.test5$xmax <- 2

stat.test5$p.signif <- c("****", "****")
figure_4B_supp <- figure_4B_supp + stat_pvalue_manual(stat.test5, tip.length = 0.02, size=2.5, label = "p.signif")

figure_4_supp <- figure_4A_supp | figure_4B_supp + 
  plot_layout( 
    heights = c(1, 0.6))

ggsave("Supplementary_figure4.pdf", figure_4_supp, width = 21/2.54, height = 21/2.54)
ggsave("Supplementary_figure4.png", figure_4_supp, width = 15.92/2.54, height = 15.92/2.54)
#######################################################################################################################################

## Metadata tuning for public repository
#######################################################################################################################################
meta_all_with_qc_curated <- merge(meta_all_with_qc_curated, table_for_plot_cor_combine[c("Sample_name", "same_cohort_NC_presence", "different_cohort_NC_presence")],
                                  all.x = T)

meta_all_with_qc_curated <- merge(meta_all_with_qc_curated, result_present_in_both[c("Sample_name", "same_cohort_NC", "different_cohort_NC")], all.x = T)

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
         timepoint_custom_category = timepoint_type)

meta_all_with_qc_curated_clean <- meta_all_with_qc_curated_clean[c("Sample_name", "Subject_ID", "FAM_ID", "Type", "status", "Timepoint", "timepoint_custom_category",
                                                                   "Study", "input_files", "in_reads_comb", "clean_reads_comb", "dedup_efficiency", "contigs_total",
                                                                   "contigs_1000", "N50", "total_viruses_discovered", "perc_reads_mapped_to_all_contigs",
                                                                   "clean_reads_mapped_to_vOTUs", "richness", "diversity", "N_shared_own_NC", "N_shared_other_NC",
                                                                   "perc_shared_own_NC", "perc_shared_other_NC", "rna_dna_liang", "incl_longitudinal")]

write.table(meta_all_with_qc_curated_clean, "/scratch/hb-llnext/VLP_public_data/nc_project/for_upload/metadata_with_qc_clean.tsv", sep='\t', row.names=F, col.names=T, quote=F)
#######################################################################################################################################
