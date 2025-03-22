##########################################
# Preparing files for strain-aware decon-
# tamination
##########################################

##############################
# Loading libraries
##############################
library(readr)
##############################
# Functions
##############################

##############################
# INPUT
##############################
sample_metadata = read.table("../metadata_with_qc_NCPv1.tsv",header=T,sep="\t")
print(dim(sample_metadata))
# list of NCs
NC_samples = sample_metadata[sample_metadata$Type == "Neg_ctrl","Sample_name"]

contig_metadata = read.table("../MERGED_Extended_TOF_NCP",header=T, sep="\t")
print(dim(contig_metadata))

rpkm <- read.delim("../RPKM_counts_VLP_NCP.txt")
print(dim(rpkm))

contig_metadata <- contig_metadata[contig_metadata$viral_genes >= contig_metadata$host_genes & 
                                     contig_metadata$POST_CHV_length >= 1000 &
                                     contig_metadata$plasmid == "No",]
print(dim(contig_metadata))

##############################
# ANALYSIS
##############################
rpkm <- rpkm[row.names(rpkm) %in% contig_metadata$New_CID,]
print(dim(rpkm))
contig_metadata <- contig_metadata[contig_metadata$New_CID %in% row.names(rpkm),]
print(dim(rpkm))
rpkm_NC <- rpkm[, colnames(rpkm) %in% NC_samples]

NC_vOTUs <- row.names(rpkm_NC[rowSums(rpkm_NC) > 0,])

rpkm <- rpkm[row.names(rpkm) %in% NC_vOTUs,]

notNC <- rpkm[,colnames(rpkm) %in% sample_metadata[sample_metadata$Type!="Neg_ctrl",]$Sample_name]
notNC <- notNC[,colSums(notNC)==0] #256 samples do not share vOTUs to NCs
ncol(notNC)/(ncol(rpkm) - ncol(rpkm_NC)) # 20.3%

rpkm <- rpkm[rowSums(rpkm>0) > 1,]
rpkm <- rpkm[,colSums(rpkm) >0]

print(dim(rpkm))
##############################
# OUTPUT
##############################
for (i in colnames(rpkm)) {
  votus_keep <- row.names(rpkm[rpkm[,i] > 0,])
  
  write.table(votus_keep, paste0("./keep_vOTUs_per_sample/", i, "_vOTUs_to_keep"), sep='\t', row.names=F, col.names=F, quote=F)
}
write.table(colnames(rpkm), "./SAMPLES_and_NCs_run_decon", sep = "\t", row.names = F, col.names = F, quote=F)

write.table(row.names(rpkm), "./vOTUs_shared_by_samples_and_ncs", sep = "\t", row.names = F, col.names = F, quote=F)
