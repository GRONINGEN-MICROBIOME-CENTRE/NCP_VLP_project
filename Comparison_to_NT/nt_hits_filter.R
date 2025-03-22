##########################################
# filtering blastn nt hits
##########################################
args = commandArgs(trailingOnly=TRUE)
##############################
# Loading libraries
##############################
library(xlsx)
library(readr)
#library(tidyverse)
library(janitor)
##############################
# Functions
##############################

##############################
# Input data
##############################
asplund <- read.xlsx2("mmc3.xlsx", 1)
colnames(asplund) <- make_clean_names(colnames(asplund))

nt_hits <- read_delim("./nt_hits.outfmt6.txt", col_names=F)
colnames(nt_hits) <- c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "qstart", "qend", "sstart", "send", "stitle", "qcovs")
nt_hits$accession <- gsub(".*\\|([A-Z0-9]+\\.\\d+)\\|.*", "\\1", nt_hits$sseqid)
##############################
# ANALYSIS
##############################
nt_hits$match_to_asplund <- "NO"
nt_hits[nt_hits$accession %in% asplund$accession,]$match_to_asplund <- "YES"
nt_hits <- nt_hits[nt_hits$qcovs >= 30,]
df <- nt_hits[nt_hits$match_to_asplund=="YES",]
nt_hits <- nt_hits[with(nt_hits, order(qseqid, -qcovs, -pident, evalue)),]
nt_hits <- nt_hits[!duplicated(nt_hits$qseqid),]

nt_hits[unique(c(grep('vir', nt_hits$stitle),
                 grep('MAG', nt_hits$stitle),
                 grep('phage', nt_hits$stitle))),]$Virus <- "YES"
##############################
# OUTPUT
##############################
write.table(nt_hits, "sorted_unique_nt_hits.txt", sep='\t', quote=F, row.names=F)
write.table(df, "redundant_hits_asplund.txt", sep='\t', quote=F, row.names=F)
