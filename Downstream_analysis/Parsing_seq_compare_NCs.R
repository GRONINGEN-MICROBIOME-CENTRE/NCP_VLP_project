library(openxlsx)
library(tidyverse)

##### On the Habrok: the result table is too heavy to work with it on PC ####
# cd /scratch/p309176/amg_paper/raw_data/NCP_studies_vir/downstream_R
# ml R
# R

.libPaths("/home1/p309176/R_libraries")
library(readr)
library(dplyr)

vir_list_NT_ext <- read_tsv("nt_hits.outfmt6.txt", col_names=F)
colnames(vir_list_NT_ext) <- c("qseqid", "Accession", "pident", "length", "qlen", "slen", "evalue", "qstart", "qend", "sstart", "send", "stitle", "qcovs")

filtered_vir_list <- vir_list_NT_ext %>%
  group_by(qseqid) %>%
  filter(evalue == min(evalue)) %>%
  filter(qcovs == max(qcovs)) %>%
  filter(pident == max(pident)) %>%
  ungroup()

write.table(filtered_vir_list, "nt_hits.outfmt6.prefiltered.txt", sep='\t', row.names=F, col.names=T, quote=F)

##### End of the script that was run on the Habrok ####
vir_list <- read.xlsx("../data/mmc3.xlsx")
vir_list_NT_ext <- read_tsv("../data/nt_hits.outfmt6.prefiltered.txt", col_names=T)


# blast_results <- read.delim("../data/result_NCs_compare_w_other_study.txt", header = F)
# colnames(blast_results) <- c("qseqid", "Accession", "pident", "length", "qlen", "slen", "evalue", "qstart", "qend", "sstart", "send", "stitle", "qcovs")

# set.seed(123)  # Set seed for reproducibility
assingnments <- vir_list_NT_ext %>%
  group_by(qseqid) %>%
  filter(length == max(length)) %>%
  ungroup()

assingnments <- assingnments %>%
  mutate(Accession = str_extract(Accession, "[^|]+(?=\\|$)"))

assingnments <- merge(assingnments, vir_list, by="Accession", all.x = T)

write.table(assingnments, "../data/assingments_NT_db.tsv", sep='\t', row.names=F, col.names=T, quote=F)
