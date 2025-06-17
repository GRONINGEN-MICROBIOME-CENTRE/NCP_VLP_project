setwd("~/Desktop/Projects_2024/AMG_paper/")

d1 = read.delim("Sample_metadata_for_cor.tsv",header=T,sep="\t")

#Liang
invrank= function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
liang = d1[d1$cohort=="liang",]
liang_d = liang[liang$Type != "Mother",]
liang_d_trans = liang_d
liang_d_trans = liang_d_trans[liang_d_trans$rna_dna=="DNA",]
liang_d_trans = liang_d_trans[!(liang_d_trans$Timepoint_numeric %in% 42),]
liang_d_trans$clean_reads_comb = invrank(liang_d_trans$clean_reads_comb)
liang_d_trans$contigs_1000 = invrank(liang_d_trans$contigs_1000)
liang_d_trans$total_viruses_discovered = invrank(liang_d_trans$total_viruses_discovered)
liang_d_trans$richness = invrank(liang_d_trans$richness)
liang_d_trans$diversity = invrank(liang_d_trans$diversity)
liang_d_trans$CHM_LU_richness_discovered_ratio = invrank(liang_d_trans$CHM_LU_richness_discovered_ratio)

l5 = lm(Timepoint_numeric ~ #diversity +
          contigs_1000 +
          #clean_reads_comb+
          total_viruses_discovered +
          CHM_LU_richness_discovered_ratio +
          richness , data = liang_d_trans[liang_d_trans$Timepoint_numeric %in% c(0,1,3,4),])
p5 = predict(l5, newdata = liang_d_trans[liang_d_trans$Type == "Neg_ctrl",])
mean(p5, na.rm = T)
sd(p5, na.rm = T)


garmaeva = d1[d1$cohort=="garmaeva",]
garmaeva_d = garmaeva[garmaeva$Type != "Mother",]
garmaeva_d_trans = garmaeva_d
garmaeva_d_trans$clean_reads_comb = invrank(garmaeva_d_trans$clean_reads_comb)
garmaeva_d_trans$contigs_1000 = invrank(garmaeva_d_trans$contigs_1000)
garmaeva_d_trans$total_viruses_discovered = invrank(garmaeva_d_trans$total_viruses_discovered)
garmaeva_d_trans$richness = invrank(garmaeva_d_trans$richness)
garmaeva_d_trans$diversity = invrank(garmaeva_d_trans$diversity)
garmaeva_d_trans$CHM_LU_richness_discovered_ratio = invrank(garmaeva_d_trans$CHM_LU_richness_discovered_ratio)


l6 = lm(Timepoint_numeric ~ #diversity +
          contigs_1000 +
          #clean_reads_comb+
          total_viruses_discovered +
          CHM_LU_richness_discovered_ratio +
          richness
        , data = garmaeva_d_trans[garmaeva_d_trans$Type == "Infant",])
p6 = predict(l6, newdata = garmaeva_d_trans[garmaeva_d_trans$Type == "Neg_ctrl",])

p6
