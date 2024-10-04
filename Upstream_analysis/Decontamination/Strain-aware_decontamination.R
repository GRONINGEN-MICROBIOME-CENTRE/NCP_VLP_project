setwd("~/Desktop/Projects_2024/AMG_paper/")

##########################################
# Strain-aware decontamination
##########################################

##############################
# Loading libraries
##############################
library(readr)
library(vegan)
library(lme4)
library(RLRsim)
library(lmerTest)
library(ggplot2)
library(ggforce)
library(tidyr)
##############################
# Functions
##############################

##############################
# INPUT
##############################
sample_metadata = read.table("metadata_with_qc_NCPv2.tsv", header=T, sep="\t")
sample_metadata$status <- ifelse(sample_metadata$Type == "Neg_ctrl", "NCs", "SAMPLES")

# list of NCs
NC_samples = sample_metadata[sample_metadata$Type == "Neg_ctrl","Sample_name"]

contig_metadata = read.table("MERGED_Extended_TOF_NCP.txt",header=T, sep="\t")
contig_metadata <- contig_metadata[contig_metadata$viral_genes >= contig_metadata$host_genes & 
                                     contig_metadata$POST_CHV_length >= 1000 &
                                     contig_metadata$plasmid == "No",]

# rpkm table:
rpkm <- read_delim("RPKM_counts_VLP_NCP.txt")
rpkm <- as.data.frame(rpkm)
row.names(rpkm) <- rpkm$V1
rpkm$V1 <- NULL

# filtering the rpkm table & etof:
rpkm <- rpkm[row.names(rpkm) %in% contig_metadata$New_CID,]
contig_metadata <- contig_metadata[contig_metadata$New_CID %in% row.names(rpkm),]

sample_metadata$richness <- colSums(rpkm>0)[match(sample_metadata$Sample_name, colnames(rpkm))]
sample_metadata[sample_metadata$clean_reads_comb >0 & is.na(sample_metadata$richness),]$richness <- 0
sample_metadata$diversity <- diversity(t(rpkm))[match(sample_metadata$Sample_name, colnames(rpkm))]

# popANI:
compare <- read.table('instrain_compare_ALL_comparisonsTable.tsv', sep='\t', header=T)
##############################
# ANALYSIS
##############################
# create dummies for NCs
studies <- unique(sample_metadata$cohort)

NCs_p_study <- list()

sample_metadata$shared_oNC <- NA

for (study in studies) {
  
  # getting the list of samples per study
  
  NCs_p_study[[study]] <- sample_metadata[sample_metadata$Type=="Neg_ctrl" & sample_metadata$cohort==study,]$Sample_name
  
  # creating a dummy NC per study
  dummy_name <- paste0( toupper(substr(study, 1,1)), 'NC' )
  if ( length(NCs_p_study[[study]]) > 1 ) {
    rpkm[,dummy_name] <- as.numeric(rowSums(rpkm[, colnames(rpkm) %in% NCs_p_study[[study]]]) > 0)
  } else {
    rpkm[,dummy_name] <- as.numeric((rpkm[,NCs_p_study[[study]]] > 0))
  }
  
  # calculating the N vOTUs shared to NCs:
  bs <- sample_metadata[sample_metadata$Type!="Neg_ctrl" & sample_metadata$cohort==study & sample_metadata$Sample_name %in% colnames(rpkm), ]$Sample_name
  NC_vOTUs <- row.names(  rpkm[rpkm[,dummy_name]!=0,]  )
  sample_metadata[sample_metadata$Sample_name %in% bs, ]$shared_oNC <- colSums( rpkm[NC_vOTUs, bs] > 0 )[match( sample_metadata[sample_metadata$Sample_name %in% bs, ]$Sample_name, bs   )]
  
}

# CLEANING THE RPKM TABLE:
DF <- compare

DF$name1 <- gsub("_filtered.sorted.bam", "", DF$name1)
DF$name2 <- gsub("_filtered.sorted.bam", "", DF$name2)

swap_indices <- DF$name1 %in% NC_samples & !DF$name2 %in% NC_samples

DF[swap_indices, c("name1", "name2")] <- DF[swap_indices, c("name2", "name1")]

# KEEPING ONLY STRAINS IDENTICAL TO NCs
DF <- DF[DF$name2 %in% NC_samples & !DF$name1 %in% NC_samples, ]
DF <- DF[DF$popANI>=0.99999 ,]

clean_rpkm <- rpkm

# REMOVING STRAINS IDENTICAL TO NCs:
for (study in studies) {
  bs <- sample_metadata[sample_metadata$Type!="Neg_ctrl" & sample_metadata$cohort==study & sample_metadata$Sample_name %in% colnames(clean_rpkm), ]$Sample_name
  
  per_virus <- DF[(DF$name2 %in% NCs_p_study[[study]]) & (DF$name1 %in% bs) ,]

  for (virus in unique(per_virus$scaffold)) {
    samples <- per_virus[per_virus$scaffold==virus,]$name1
    clean_rpkm [row.names(clean_rpkm)==virus, colnames(clean_rpkm) %in% samples] <- 0
    
  }
  
}

# CALCULATING NEW RICHNESS & DIVERSITY
sample_metadata$richness_clean <- colSums(clean_rpkm>0)[match(sample_metadata$Sample_name, colnames(clean_rpkm))]
sample_metadata$diversity_clean <- diversity(t(clean_rpkm))[match(sample_metadata$Sample_name, colnames(clean_rpkm))]

# CALCULATING NEW SHAREDNESS TO NCs
sample_metadata$shared_oNC_clean <- NA

# CALCULATING NEW SHAREDNESS TO NCs
for (study in studies) {
  
  # creating a dummy NC per study
  dummy_name <- paste0( toupper(substr(study, 1,1)), 'NC' )
  # calculating the N vOTUs shared to NCs after cleaning:
  bs <- sample_metadata[sample_metadata$Type!="Neg_ctrl" & sample_metadata$cohort==study & sample_metadata$Sample_name %in% colnames(clean_rpkm), ]$Sample_name
  NC_vOTUs <- row.names(  clean_rpkm[clean_rpkm[,dummy_name]!=0,]  )
  sample_metadata[sample_metadata$Sample_name %in% bs, ]$shared_oNC_clean <- colSums( clean_rpkm[NC_vOTUs, bs] > 0 )[match( sample_metadata[sample_metadata$Sample_name %in% bs, ]$Sample_name, bs   )]
  
}

# RPKM TABLE FOR NCs
cols_select <- (colnames(rpkm) %in%NC_samples) | (colnames(rpkm) %in% paste0( toupper(substr(studies, 1,1)), 'NC' ))
rpkm_NC <- rpkm[unique(compare$scaffold), cols_select]
rpkm_NC <- rpkm_NC[,colSums(rpkm_NC)>0]
#shared_in_NCs <- names(which(rowSums(rpkm_NC[,grep('NC', colnames(rpkm_NC))] > 0) > 1 ) )

# STRAIN IDENTITY FOR vOTUs IDENTIFIED IN NCs:
NC_stid <- compare
NC_stid$name1 <- gsub("_filtered.sorted.bam", "", NC_stid$name1)
NC_stid$name2 <- gsub("_filtered.sorted.bam", "", NC_stid$name2)
NC_stid <- NC_stid[NC_stid$name1 %in% NC_samples & NC_stid$name2 %in% NC_samples,]

# samples that are the same in the NCs from different cohorts
unique(c(NC_stid[NC_stid$percent_genome_compared >= 0.75 & NC_stid$popANI >= 0.99999, ]$name1, NC_stid[NC_stid$percent_genome_compared >= 0.75 & NC_stid$popANI >= 0.99999, ]$name2))

# HOW MANY STRAINS ARE SHARED OWN vs OTHER COHORTS
str_sharing_NCs <- as.data.frame(matrix(NA, nrow=length(studies), ncol=3), row.names = studies)
colnames(str_sharing_NCs) <- c("N_sp_detected", "N_st_shared_own", "N_st_shared_other")

for (study in studies) {
  dummy_name <- paste0( toupper(substr(study, 1,1)), 'NC' )
  
  str_sharing_NCs[study,"N_sp_detected"] <- sum(rpkm_NC[,dummy_name]==1)
if (length(NCs_p_study[[study]]) > 1) {

  
  shared_in_NCs <- unique(NC_stid[NC_stid$percent_genome_compared >= 0.75 & 
                                    NC_stid$popANI >= 0.99999 & 
                                    NC_stid$name1 %in% NCs_p_study[[study]] & 
                                    NC_stid$name2 %in% NCs_p_study[[study]],]$scaffold)
  
  str_sharing_NCs[study,"N_st_shared_own"] <- length(shared_in_NCs)
}
  
  shared_in_NCs <- unique(NC_stid[NC_stid$percent_genome_compared >= 0.75 & 
                                    NC_stid$popANI >= 0.99999 & 
                                    (NC_stid$name1 %in% NCs_p_study[[study]] & (
                                    !NC_stid$name2 %in% NCs_p_study[[study]] ) | 
                                    (!NC_stid$name1 %in% NCs_p_study[[study]] & 
                                       NC_stid$name2 %in% NCs_p_study[[study]] )),]$scaffold)
  
  str_sharing_NCs[study,"N_st_shared_other"] <- length(shared_in_NCs)
}

str_sharing_NCs$N_st_shared_own_perc <- str_sharing_NCs$N_st_shared_own/str_sharing_NCs$N_sp_detected

# STRAIN-SHARING IN SAMPLES, STRICT CUT-OFFs:
DF_cov_included <- DF[DF$coverage_overlap >= 0.75,]
length(unique(DF_cov_included$scaffold))

# IS THERE AN ENRICHMENT FOR STRAINS IDENTICAL TO 

# IF CLEANED AT THE SPECIES-LEVEL:
sp_clean <- rpkm
for (study in studies) {
  
  # creating a dummy NC per study
  dummy_name <- paste0( toupper(substr(study, 1,1)), 'NC' )
  # calculating the N vOTUs shared to NCs after cleaning:
  bs <- sample_metadata[sample_metadata$Type!="Neg_ctrl" & sample_metadata$cohort==study & sample_metadata$Sample_name %in% colnames(sp_clean), ]$Sample_name
  NC_vOTUs <- row.names(  sp_clean[sp_clean[,dummy_name]!=0,]  )
  sp_clean[NC_vOTUs, bs] <- 0
}

sample_metadata$richness_sp_clean <- colSums(sp_clean>0)[match(sample_metadata$Sample_name, colnames(sp_clean))]

# CALCULATING DELTA
sample_metadata$delta_st <- sample_metadata$richness - sample_metadata$richness_clean
sample_metadata$delta_sp <- sample_metadata$richness - sample_metadata$richness_sp_clean
sample_metadata$delta_st_sp <- sample_metadata$richness_clean - sample_metadata$richness_sp_clean

sample_metadata$delta_st_perc <- (sample_metadata$richness - sample_metadata$richness_clean)/sample_metadata$richness
sample_metadata$delta_sp_perc <- (sample_metadata$richness - sample_metadata$richness_sp_clean)/sample_metadata$richness

# CHANGE IN RICHNESS:
summary(sample_metadata[sample_metadata$Type %in% c('Infant', 'Mother')  & !is.na(sample_metadata$shared_oNC) & sample_metadata$shared_oNC>0,]$delta_st_perc*100)
# median is 1.5%, (IQR: 0.5-5)
summary(sample_metadata[sample_metadata$Type %in% c('Infant', 'Mother')  & !is.na(sample_metadata$shared_oNC) & sample_metadata$shared_oNC>0,]$delta_sp_perc*100)
# median is 4.9%, (IQR: 1.9-12.9)

# CHANGE IN N SHARED TO NCs:
sample_metadata$shared_oNC_clean_perc <- sample_metadata$shared_oNC_clean/sample_metadata$richness_clean
sample_metadata$shared_oNC_perc <- sample_metadata$shared_oNC/sample_metadata$richness
sample_metadata$shared_oNC_decrease_perc <- (sample_metadata$shared_oNC - sample_metadata$shared_oNC_clean)/sample_metadata$shared_oNC
summary(sample_metadata[sample_metadata$Type %in% c('Infant', 'Mother')  & !is.na(sample_metadata$shared_oNC) & sample_metadata$shared_oNC>0,,]$shared_oNC_clean_perc*100)
# 2.9%, (IQR: 1.0-7.8)

summary(sample_metadata[sample_metadata$Type %in% c('Infant', 'Mother')  & !is.na(sample_metadata$shared_oNC) & sample_metadata$shared_oNC>0,,]$shared_oNC_decrease_perc*100)
# 33.3%, (IQR: 14.9-50.0)

# how many samples shared at least 1 strain to NC?
sum(sample_metadata$delta_st_perc>0, na.rm=T)
sum(sample_metadata$shared_oNC_perc>0, na.rm=T)

##############################
# OUTPUT
##############################

