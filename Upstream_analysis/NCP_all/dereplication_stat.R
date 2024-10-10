##########################################
# Parsing VC data and VC size after
# dereplication
##########################################
args = commandArgs(trailingOnly=TRUE)
##############################
# Loading libraries
##############################
.libPaths("/home1/p309176/R_libraries")
library(tidyverse)
##############################
# Functions
##############################

##############################
# Input data
##############################
viral_clusters <- read_delim("./NCP_viral_clusters_clean.tsv", col_names = F)
##############################
# ANALYSIS
##############################

# Transforming into long format to merge with Extended_TOF
long_df <- viral_clusters %>%
  separate_rows(X2, sep = ",\\s*") %>%
  arrange(X2, X1) %>%
  rename(Representative = X1, Cluster_member = X2)

MERGED_Extended_TOF_NCP <- read.delim("./MERGED_Extended_TOF_NCP")
MERGED_Extended_TOF_NCP$Cluster_member <- MERGED_Extended_TOF_NCP$New_CID
MERGED_Extended_TOF_NCP <- merge(MERGED_Extended_TOF_NCP, long_df, by="Cluster_member", all.y = T)
MERGED_Extended_TOF_NCP$Cluster_member <- NULL

MERGED_Extended_TOF_NCP_filt <- MERGED_Extended_TOF_NCP[MERGED_Extended_TOF_NCP$POST_CHV_length >= 1000 & (MERGED_Extended_TOF_NCP$viral_genes >= MERGED_Extended_TOF_NCP$host_genes) & MERGED_Extended_TOF_NCP$plasmid=="No", ]

long_df <- long_df[long_df$Cluster_member %in% MERGED_Extended_TOF_NCP_filt$New_CID, ]
write.table(long_df, './NCP_viral_clusters_clean_long_format.txt', sep='\t', row.names = F, quote=F)

# Getting the size of dereplication clusters
cluster_counts <- long_df %>%
  count(Representative) %>%
  rename(Cluster_size = n)

MERGED_Extended_TOF_NCP_filt <- merge(MERGED_Extended_TOF_NCP_filt, cluster_counts, by="Representative", all.x = T)
lines_vector <- readLines("./virus_representatives_NCs_filtered.txt")
garmaeva_ncs_vector <-  ("./vOTUs_present_in_garmaeva_NCs")
liang_ncs_vector <- readLines("./vOTUs_present_in_liang_NCs")
maqsood_ncs_vector <- readLines("./vOTUs_present_in_maqsood_NCs")
shah_ncs_vector <- readLines("./vOTUs_present_in_shah_NCs")

MERGED_Extended_TOF_NCP_filt <- MERGED_Extended_TOF_NCP_filt %>% 
	mutate(vOTUr_present_in_NC = ifelse(Representative %in% lines_vector, "Yes", "No"),
	       Study_of_sample_origin = str_extract(New_CID, "^[^_]+"),
	       plasmid = NULL,
	       conjugation_genes = NULL,
	       amr_genes = NULL,
	       plasmid_score = NULL,
	       plasmid_fdr = NULL) %>%
	rename(vOTU_cluster_size = Cluster_size,
	       vOTUr_ID = Representative,
	       virus_ID = New_CID,
	       final_length = POST_CHV_length) %>%
	select(1:2, vOTU_cluster_size, everything())

Extended_TOF_NCP_clean_representative_all <- MERGED_Extended_TOF_NCP_filt[MERGED_Extended_TOF_NCP_filt$vOTUr_ID == MERGED_Extended_TOF_NCP_filt$virus_ID, ]
Extended_TOF_NCP_clean_representative_all <- Extended_TOF_NCP_clean_representative_all %>%
	mutate(virus_ID = NULL,
	       Study_of_sample_origin = NULL)

Extended_TOF_NCP_clean_redundant_NCs <- MERGED_Extended_TOF_NCP_filt[MERGED_Extended_TOF_NCP_filt$vOTUr_present_in_NC == "Yes", ]
Extended_TOF_NCP_clean_redundant_NCs <- Extended_TOF_NCP_clean_redundant_NCs %>%
	mutate(cluster_status = ifelse(vOTUr_ID == virus_ID, "representative", "member"),
	       vOTUr_present_in_NC_garmaeva = ifelse(vOTUr_ID %in% garmaeva_ncs_vector, "Yes", "No"),
	       vOTUr_present_in_NC_liang = ifelse(vOTUr_ID %in% liang_ncs_vector, "Yes", "No"),
	       vOTUr_present_in_NC_maqsood = ifelse(vOTUr_ID %in% maqsood_ncs_vector, "Yes", "No"),
	       vOTUr_present_in_NC_shah = ifelse(vOTUr_ID %in% shah_ncs_vector, "Yes", "No")) %>%
	select(1:2, cluster_status, everything())


write.table(Extended_TOF_NCP_clean_representative_all, "./Extended_TOF_NCP_clean_representative_all.tsv", sep='\t', row.names = F, quote=F)
write.table(MERGED_Extended_TOF_NCP_filt, "./Extended_TOF_NCP_clean_redundant_all.tsv", sep='\t', row.names = F, quote=F)
write.table(Extended_TOF_NCP_clean_redundant_NCs, "./Extended_TOF_NCP_clean_redundant_NCs.tsv", sep='\t', row.names = F, quote=F)

writeLines(MERGED_Extended_TOF_NCP_filt$New_CID, con = "./viruses_redundant_all_clean.txt")
writeLines(Extended_TOF_NCP_clean_redundant_NCs$New_CID, con = "./viruses_redundant_NCs_clean.txt")

