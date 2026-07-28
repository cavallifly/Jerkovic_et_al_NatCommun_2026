rm(list=ls())
library(ggpubr)

setwd("C://Users/pc/Documents/Thèse/Temp from Archive/")

#Import WT data
NPC_WT <- rbind(read.csv2("Merge_summary/All_replicates/NPC_wt_rep1_merge.csv", header = T, stringsAsFactors = F),
                read.csv2("Merge_summary/All_replicates/NPC_WT_rep3_merge.csv", header = T, stringsAsFactors = F))
NPC_WT[,-1]<- apply(NPC_WT[,-1],2, as.numeric)

#Calculate distances then remove dist>2.2µm for NPC WT
NPC_WT$F1F2 <- sqrt(((NPC_WT$xFC2-NPC_WT$xFC1)*40)**2+((NPC_WT$yFC2-NPC_WT$yFC1)*40)**2+((NPC_WT$zFC2-NPC_WT$zFC1)*140)**2)
NPC_WT$F2F3 <- sqrt(((NPC_WT$xFC3-NPC_WT$xFC2)*40)**2+((NPC_WT$yFC3-NPC_WT$yFC2)*40)**2+((NPC_WT$zFC3-NPC_WT$zFC2)*140)**2)
NPC_WT$F1F3 <- sqrt(((NPC_WT$xFC3-NPC_WT$xFC1)*40)**2+((NPC_WT$yFC3-NPC_WT$yFC1)*40)**2+((NPC_WT$zFC3-NPC_WT$zFC1)*140)**2)

cut_off_mydist<- function(a){
  a <- a[a$F1F2 %in% subset(a$F1F2, a$F1F2 < 2200),]
  a <- a[a$F2F3 %in% subset(a$F2F3, a$F2F3 < 2200),]
  a <- a[a$F1F3 %in% subset(a$F1F3, a$F1F3 < 2200),]
}
NPC_WT <- cut_off_mydist(NPC_WT)


#Import dA data
NPC_dA <- rbind(read.csv2("Merge_summary/All_replicates/NPC_dA_rep1_merge.csv", header = T, stringsAsFactors = F),
                read.csv2("Merge_summary/All_replicates/NPC_dA_rep2_merge.csv", header = T, stringsAsFactors = F),
                read.csv2("Merge_summary/All_replicates/NPC_dA_rep3_merge.csv", header = T, stringsAsFactors = F))
NPC_dA[,-1]<- apply(NPC_dA[,-1],2, as.numeric)

#Calculate distances then remove dist>2.2µm for NPC dA
NPC_dA$F1F2 <- sqrt(((NPC_dA$xFC2-NPC_dA$xFC1)*40)**2+((NPC_dA$yFC2-NPC_dA$yFC1)*40)**2+((NPC_dA$zFC2-NPC_dA$zFC1)*140)**2)
NPC_dA$F2F3 <- sqrt(((NPC_dA$xFC3-NPC_dA$xFC2)*40)**2+((NPC_dA$yFC3-NPC_dA$yFC2)*40)**2+((NPC_dA$zFC3-NPC_dA$zFC2)*140)**2)
NPC_dA$F1F3 <- sqrt(((NPC_dA$xFC3-NPC_dA$xFC1)*40)**2+((NPC_dA$yFC3-NPC_dA$yFC1)*40)**2+((NPC_dA$zFC3-NPC_dA$zFC1)*140)**2)
NPC_dA <- cut_off_mydist(NPC_dA)


#Calculate area for WT and dA
NPC_dA_484 <- NPC_dA[1:484,]

#Perimeter
NPC_WT$perimeter <- (NPC_WT$F1F2+ NPC_WT$F1F3+ NPC_WT$F2F3)
NPC_dA_484$perimeter <- (NPC_dA_484$F1F2+ NPC_dA_484$F1F3+ NPC_dA_484$F2F3)

#Semi-Perimeter
NPC_WT$semiperimeter <- (NPC_WT$perimeter/2)
NPC_dA_484$semiperimeter <- (NPC_dA_484$perimeter/2)

#Area
NPC_WT$area <- sqrt(NPC_WT$semiperimeter*(NPC_WT$semiperimeter - NPC_WT$F1F2)*(NPC_WT$semiperimeter - NPC_WT$F2F3)*(NPC_WT$semiperimeter - NPC_WT$F1F3))
NPC_dA_484$area <- sqrt(NPC_dA_484$semiperimeter*(NPC_dA_484$semiperimeter - NPC_dA_484$F1F2)*(NPC_dA_484$semiperimeter - NPC_dA_484$F2F3)*(NPC_dA_484$semiperimeter - NPC_dA_484$F1F3))

NPC_dA_484$condition <- "dA"
NPC_WT$condition <- "_WT"

AllCondition<- rbind(NPC_WT, NPC_dA_484)
AllCondition$area <- AllCondition$area/1000

pdf("R_outpout/area_all_condition.pdf",width = 5, height = 5)
ggviolin(AllCondition, x="condition", y="area", fill="condition", palette = "Dark2",
         add=c("boxplot"), add.params=list(fill="white", width=0.1),
         width = 0.5, order=c("_WT", "dA"), title = "All points")+
  stat_compare_means( method = "wilcox.test" )
dev.off()
