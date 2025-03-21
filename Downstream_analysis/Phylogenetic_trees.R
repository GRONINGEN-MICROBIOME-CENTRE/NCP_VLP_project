setwd("~/Desktop/Projects_2024/AMG_paper/")

##########################################
# Calculating N shared strains between NCs 
# and biological samples;
# Reconstructing the trees based on popANI
# Authors: Alex Kurilshikov & Sana Garmaeva
##########################################

##############################
# Loading libraries
##############################
library(ggtree)
library(ggplot2)
library(dplyr)
library(ape)
library(stringr)
library(MetBrewer)
library(readr)
library(circlize)
##############################
# Functions
##############################
Kandinsky <- met.brewer('Kandinsky')
Egypt <- met.brewer('Egypt')

# fill_symmetric_nas function
# This function fills NA values in a symmetrical matrix by using the corresponding symmetrical values.
# If mat[i, j] is NA and mat[j, i] is not NA, it sets mat[i, j] to mat[j, i], and vice versa.
# The function returns the matrix with NA values filled based on symmetry.
fill_symmetric_nas <- function(mat) {
  n <- nrow(mat)
  for (i in 1:n) {
    for (j in i:n) {
      if (is.na(mat[i, j]) && !is.na(mat[j, i])) {
        mat[i, j] <- mat[j, i]
      } else if (!is.na(mat[i, j]) && is.na(mat[j, i])) {
        mat[j, i] <- mat[i, j]
      }
    }
  }
  return(mat)
}

# make_Dmat function
# This function transforms popANI values into distances between samples and reconstructs a symmetrical matrix,
# using these distances

make_Dmat = function(x){
# making a matrix of pre-distances
  mat = matrix(0,nrow = length(unique(c(x$name1,x$name2))),ncol = length(unique(c(x$name1,x$name2))))
  
  rownames(mat) = sub("_w_neg_der95_NCP.sorted.bam","",unique(c(x$name1,x$name2)))
  colnames(mat) = sub("_w_neg_der95_NCP.sorted.bam","",unique(c(x$name1,x$name2)))
  
  mat[lower.tri(mat)] = 1-x[,11]
  mat[upper.tri(mat)] = 1-x[,11]
  
# filling NAs if a corresponding pair has non-NA value
  mat <- fill_symmetric_nas(mat)
  
# removal of rows containing NAs yet keeping the most of the matrix:
  while (any(is.na(mat))) {
    
    # Count non-NA values in each row and column
    row_non_na <- rowSums(!is.na(mat))
    col_non_na <- colSums(!is.na(mat))
    
    # Find the minimum non-NA count in rows and columns
    min_row_non_na <- min(row_non_na[row_non_na > 0])
    min_col_non_na <- min(col_non_na[col_non_na > 0])
    
    # Remove the row or column with the least non-NA values
    if (min_row_non_na <= min_col_non_na) {
      row_to_remove <- which(row_non_na == min_row_non_na)[1]
      mat <- mat[-row_to_remove, -row_to_remove]
    } else {
      col_to_remove <- which(col_non_na == min_col_non_na)[1]
      mat <- mat[-col_to_remove, -col_to_remove]
    }
    
  }
  return(as.dist(mat))
  
}

make_plot = function(x,offset = 0.0001,size = 2.5,xlimMulti = 3.5){
  phylo1 <- as.phylo(hclust(as.dist(x[[1]]),method="complete"))
  
  metadata_subset = sample_metadata[match(phylo1$tip.label,sample_metadata$Sample_name),]
  
  colnames(metadata_subset)[grep('Sample_name$', colnames(metadata_subset))] = 'label'
  
  phylo2 = full_join(phylo1, metadata_subset, by = 'label')
  
  ggtree(phylo2,options(ignore.negative.edge = T),layout = "rectangular") + 
    labs(tag="a") +
    geom_tippoint(aes(shape=Type, fill=cohort),size = 2) +
    scale_shape_manual(values=c("Infant"=21,"NC"=24)) +
    scale_fill_manual(values = c(`Garmaeva *et al*., 2024`=Kandinsky[4],
                                   `Liang *et al*., 2020`=Kandinsky[2],
                                   `Maqsood *et al*., 2019`=Kandinsky[3],
                                   `Shah *et al*., 2023`=Kandinsky[1]), drop=T)  +
    labs(shape="Sample type", fill="Study") +
    ggtitle(gsub("maqsood_bctrl4633v_N1_L5513_K159.5_E0_P0_F0", "phi X 174", names(x))) +
    #geom_treescale(offset=2, fontsize=3) +
    geom_treescale(fontsize=3) + # for phix
    theme(legend.text = ggtext::element_markdown(size=7), 
          legend.title = element_text(size=7),
          plot.tag = element_text(face="bold"),
          plot.title = element_text(size=7),
          plot.title.position = "plot",
          legend.position = c(0.1,0.7), # for phix
          #legend.position = c(0.7,0.7),
          legend.background = element_rect(fill=NA),
          legend.direction = "vertical", legend.box = "vertical") + # for phix
          #legend.direction = "vertical", legend.box = "horizontal") +
    guides(fill=guide_legend(override.aes=list(shape=21))) + 
    xlim(0,max(phylo1$edge.length) * xlimMulti) #+ 
    #scale_x_reverse() + 
    #coord_flip() 
  
}


##############################
# INPUT
##############################
contig_metadata = read.table("MERGED_Extended_TOF_NCP.txt",header=T, sep="\t")

sample_metadata = read.table("metadata_with_qc_NCPv1.tsv",header=T,sep="\t")
sample_metadata[sample_metadata$Type=="Neg_ctrl",]$Type <- "NC"
sample_metadata$Type = factor(sample_metadata$Type)
sample_metadata[sample_metadata$cohort=="garmaeva",]$cohort <- "Garmaeva *et al*., 2024"
sample_metadata[sample_metadata$cohort=="shah",]$cohort <- "Shah *et al*., 2023"
sample_metadata[sample_metadata$cohort=="liang",]$cohort <- "Liang *et al*., 2020"
sample_metadata[sample_metadata$cohort=="maqsood",]$cohort <- "Maqsood *et al*., 2019"
sample_metadata$cohort <- factor(sample_metadata$cohort)

# list of NCs
NC_samples = sample_metadata[sample_metadata$Type == "NC","Sample_name"]

# InStrain output (cut-offs: 1X coverage):
compare = read.table("instrain_compare2_comparisonsTable.tsv",header=T)
compare2 <- read.table("add_instrain_compare3_comparisonsTable.tsv", header=T)
compare <- rbind(compare, compare2)
rm(compare2)

compare$name1 <- gsub("_w_neg_der95_NCP.sorted.bam", "", compare$name1)
compare$name2 <- gsub("_w_neg_der95_NCP.sorted.bam", "", compare$name2)

# excluding sequences comparisons that produced less than 75% genome length overlap to denoise trees:
compare[compare$percent_genome_compared < 0.75,"popANI"] <- NA

# predicted prokaryotic hosts:
hosts <- read.table('MERGED_Host_prediction_to_genus_m90_v2.csv', sep=',', header=T)
##############################
# ANALYSIS
##############################


# splitting data per sequence:
compare.split = split(compare,compare$scaffold)

working_contigs = compare.split

# filtering contigs metadata to focus on the sequences of interest:
check_contig_metadata <- contig_metadata[contig_metadata$New_CID %in% names(compare.split),]

# filtering out sequences with low number of virus genes:
check_contig_metadata <- check_contig_metadata[check_contig_metadata$host_genes <= check_contig_metadata$viral_genes,]
working_contigs <- working_contigs[names(working_contigs) %in% check_contig_metadata$New_CID]

# calculating distances from (1 - popANI)
dmats = lapply(working_contigs,make_Dmat)

dmats <- Filter(function(x) length(x) >= 2, dmats)

check_contig_metadata <- check_contig_metadata[check_contig_metadata$New_CID %in% names(dmats),]

# filtering host prediction table:
hosts <- hosts[hosts$Virus %in% check_contig_metadata$New_CID,]

# filtering compare table:
compare <- compare[compare$scaffold %in% check_contig_metadata$New_CID,]

make_plot(dmats["maqsood_bctrl4633v_N1_L5513_K159.5_E0_P0_F0"], offset = 0.000005, xlimMulti = 1) 

# adding same vs diff to NCs strains info to the sample metadata to place the tree highlights properly:
phages <- c('maqsood_bctrl4633v_N1_L5513_K159.5_E0_P0_F0', #phix
            'maqsood_C0251iv_N4_L41225_K60.7_E0_P0_F0', # burkholderia phage
            'shah_kid122_N774_L6428_K6.1_E0_P0_F0' # Microviridae, Bacteroides phage
            )

for (i in phages) {
  
  # filtering the compare data.frame
  DF <- compare[compare$scaffold==i,]
  swap_indices <- DF$name1 %in% NC_samples & !DF$name2 %in% NC_samples
  
  DF[swap_indices, c("name1", "name2")] <- DF[swap_indices, c("name2", "name1")]
  
  DF <- DF[DF$name2 %in% NC_samples & !DF$name1 %in% NC_samples, ]
  
  # creating factor for coloring:
  sample_metadata[, stringr::str_extract(i, 'L[0-9]+')] <- 'Different'
  
  if (i == "maqsood_C0251iv_N4_L41225_K60.7_E0_P0_F0") {
    
    sample_metadata[sample_metadata$Sample_name %in% DF[DF$popANI >= 0.9999,]$name1, stringr::str_extract(i, 'L[0-9]+')] <- 'Same'
  
    } else {
    
    sample_metadata[sample_metadata$Sample_name %in% DF[DF$popANI >= 0.99999,]$name1, stringr::str_extract(i, 'L[0-9]+')] <- 'Same'
    
  }
  
  
}

# phiX tree:

x <- dmats["maqsood_bctrl4633v_N1_L5513_K159.5_E0_P0_F0"]

offset =  0.000005
xlimMulti = 1

phylo1 <- as.phylo(hclust(as.dist(x[[1]]),method="complete"))

metadata_subset = sample_metadata[match(phylo1$tip.label,sample_metadata$Sample_name),]
colnames(metadata_subset)[grep('Sample_name$', colnames(metadata_subset))] = 'label'
phylo2 = full_join(phylo1, metadata_subset, by = 'label')

phix <- ggtree(phylo2,options(ignore.negative.edge = T),layout = "rectangular") + 
  labs(tag="a", fill="Strains") +
  geom_highlight(node=16, fill="#810955", alpha=0.4, type="auto") +
  geom_rect(aes(xmin=0, 
                xmax=0,
                ymin=0, ymax=0,
                fill = "Same as in NCs"),alpha=0.4) + 
  scale_fill_manual(values = c("Same as in NCs" = "#810955"), drop=T) +
  scale_shape_manual(values=c("Infant"=21,"NC"=24)) +
  ggnewscale::new_scale_fill() +
  #geom_tippoint(aes(shape=Type, fill=phiX),size = 2) +
  geom_tippoint(aes(shape=Type, fill=cohort),size = 2) +
  scale_fill_manual(values = c(`Garmaeva *et al*., 2024`=Kandinsky[4],
                               `Liang *et al*., 2020`=Kandinsky[2],
                               `Maqsood *et al*., 2019`=Kandinsky[3],
                               `Shah *et al*., 2023`=Kandinsky[1]), drop=T)  +
  labs(shape="Sample type", fill="Study") +
  ggtitle(gsub("maqsood_bctrl4633v_N1_L5513_K159.5_E0_P0_F0", "phi X 174", names(x))) +
  geom_treescale(fontsize=3, offset = 0.5) + # for phix
  theme(legend.text = ggtext::element_markdown(size=8), 
        legend.title = element_text(size=9),
        plot.tag = element_text(face="bold"),
        plot.title = element_text(size=9),
        plot.title.position = "plot",
        #legend.position = c(0.1,0.7), # for phix
        legend.background = element_rect(fill=NA),
        legend.direction = "vertical", legend.box = "vertical",
        #legend.box.background = element_rect(
        #  fill = '#FFFFFF', size = 0.0, linetype = "solid"),
        #legend.spacing.y = unit(0.5, "lines"),
        #legend.margin = margin(0, 0, 0, 0)
        ) +
  guides(fill=guide_legend(order=1, override.aes=list(shape=21))) + 
  xlim(0,max(phylo1$edge.length) * xlimMulti)


# burkholeria tree:
x <- dmats["maqsood_C0251iv_N4_L41225_K60.7_E0_P0_F0"]

offset = 0.5
xlimMulti = 1.4
 
phylo1 <- as.phylo(hclust(as.dist(x[[1]]),method="complete"))
metadata_subset = sample_metadata[match(phylo1$tip.label,sample_metadata$Sample_name),]
colnames(metadata_subset)[grep('Sample_name$', colnames(metadata_subset))] = 'label'
phylo2 = full_join(phylo1, metadata_subset, by = 'label')
 
burkholderia <- ggtree(phylo2,options(ignore.negative.edge = T),layout = "rectangular") + 
   labs(tag="a", fill="Strains") +
  geom_highlight(node=66, fill="#810955", alpha=0.4) +
  geom_highlight(node=69, fill="#810955", alpha=0.4) +
  geom_highlight(node=72, fill="#810955", alpha=0.4) +
  geom_highlight(node=61, fill="#810955", alpha=0.4) +
  geom_highlight(node=63, fill="#810955", alpha=0.4) +
  geom_highlight(node=67, fill="#810955", alpha=0.4) +
   geom_rect(aes(xmin=0, 
                 xmax=0,
                 ymin=0, ymax=0,
                 fill = "Same as in NCs"),alpha=0.4) + 
   scale_fill_manual(values = c("Same as in NCs" = "#810955"), drop=T) +
   scale_shape_manual(values=c("Infant"=21,"NC"=24)) +
   ggnewscale::new_scale_fill() +
   #geom_tippoint(aes(shape=Type, fill=L41225),size = 2) +
   geom_tippoint(aes(shape=Type, fill=cohort),size = 2) +
     scale_fill_manual(values = c(`Garmaeva *et al*., 2024`=Kandinsky[4],
                                  `Liang *et al*., 2020`=Kandinsky[2],
                                  `Maqsood *et al*., 2019`=Kandinsky[3],
                                  `Shah *et al*., 2023`=Kandinsky[1]), drop=T)  +
   labs(shape="Sample type", fill="Study") +
   ggtitle(names(x)) +
   geom_treescale(offset=2, fontsize=3, x = 0.00005) +
   theme(legend.text = ggtext::element_markdown(size=8), 
         legend.title = element_text(size=9),
         plot.tag = element_text(face="bold"),
         plot.title = element_text(size=9),
         plot.title.position = "plot",
         legend.position = c(0.6,0.7),
         legend.background = element_rect(fill=NA),
         legend.direction = "vertical", legend.box = "horizontal") +
   guides(fill=guide_legend(order=1, override.aes=list(shape=21))) + 
   scale_x_reverse() + 
   coord_flip() 
 

# Bacteroides phage, Microviridae:
x <- dmats["shah_kid122_N774_L6428_K6.1_E0_P0_F0"]

offset = 1
xlimMulti = 1

phylo1 <- as.phylo(hclust(as.dist(x[[1]]),method="complete"))
metadata_subset = sample_metadata[match(phylo1$tip.label,sample_metadata$Sample_name),]
colnames(metadata_subset)[grep('Sample_name$', colnames(metadata_subset))] = 'label'
phylo2 = full_join(phylo1, metadata_subset, by = 'label')

micro_bacteroides <- ggtree(phylo2,options(ignore.negative.edge = T),layout = "rectangular") + 
  labs(tag="a", fill="Strains") +
  geom_highlight(node=37, fill="#810955", alpha=0.4) +
  geom_highlight(node=44, fill="#810955", alpha=0.4) +
  geom_rect(aes(xmin=0, 
                xmax=0,
                ymin=0, ymax=0,
                fill = "Same as in NCs"),alpha=0.4) + 
  scale_fill_manual(values = c("Same as in NCs" = "#810955"), drop=T) +
  scale_shape_manual(values=c("Infant"=21,"NC"=24)) +
  ggnewscale::new_scale_fill() +
  #geom_tippoint(aes(shape=Type, fill=L6428),size = 2) +
  geom_tippoint(aes(shape=Type, fill=cohort),size = 2) +
   scale_fill_manual(values = c(`Garmaeva *et al*., 2024`=Kandinsky[4],
                                 `Liang *et al*., 2020`=Kandinsky[2],
                                 `Maqsood *et al*., 2019`=Kandinsky[3],
                                 `Shah *et al*., 2023`=Kandinsky[1]), drop=T)  +
  labs(shape="Sample type", fill="Study") +
  ggtitle(names(x)) +
  geom_treescale(offset=2, fontsize=3, x = 0.00005) +
  theme(legend.text = ggtext::element_markdown(size=8), 
        legend.title = element_text(size=9),
        plot.tag = element_text(face="bold"),
        plot.title = element_text(size=9),
        plot.title.position = "plot",
        legend.position = c(0.6,0.7),
        legend.background = element_rect(fill=NA),
        legend.direction = "vertical", legend.box = "horizontal") +
  guides(fill=guide_legend(order=1, override.aes=list(shape=21))) + 
  scale_x_reverse() + 
  coord_flip() 



##############################
# OUTPUT
##############################
# save for patching:
saveRDS(phix, file = "phix.rds")

# save for patching:
saveRDS(burkholderia, file = "burkholderia.rds")

# save for patching:
saveRDS(micro_bacteroides, file = "micro_bacteroides.rds")


