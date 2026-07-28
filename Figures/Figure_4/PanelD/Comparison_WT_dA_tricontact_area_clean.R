rm(list=ls())
library(ggplot2)
library(ggpubr)
setwd("C://Users/pc/Documents/Thèse/Temp from Archive/")

#Set up the threshold for the Tricontacts (ttc)
My_Threshold <- 400

#Import AfterAPC for WT and dA
NPC_WT_AfterAPC <- read.csv2("R_outpout/NPC_WT/NPC_WT_AfterAPC.csv", stringsAsFactors = F)
NPC_dA_AfterAPC <- read.csv2("R_outpout/NPC_dA//NPC_dA_AfterAPC.csv", stringsAsFactors = F)
NPC_dA_AfterAPC[,2] <- NPC_dA_AfterAPC[,2]*-1
NPC_dA_AfterAPC[,3] <- NPC_dA_AfterAPC[,3]*-1


#Calculate the distances between the points after PCA : AfterAPC

#For NPC_WT
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

tAfterAPC <- data.frame(ID=rep(NA,NROW(NPC_WT)),
                        Xprom=rep(NA,NROW(NPC_WT)),
                        Yprom=rep(NA,NROW(NPC_WT)),
                        XEnhA=rep(NA,NROW(NPC_WT)),
                        YEnhA=rep(NA,NROW(NPC_WT)),
                        XEnhb=rep(NA,NROW(NPC_WT)),
                        YEnhb=rep(NA,NROW(NPC_WT)))
u <- 1
for (i in seq(1,NROW(NPC_WT_AfterAPC),by=3)) {
  tAfterAPC[u,1] <- NPC_WT_AfterAPC[i,1]
  tAfterAPC[u,2:3] <- NPC_WT_AfterAPC[i,2:3]
  tAfterAPC[u,4:5] <- NPC_WT_AfterAPC[i+1,2:3]
  tAfterAPC[u,6:7] <- NPC_WT_AfterAPC[i+2,2:3]
  u <- u+1
}

tAfterAPC$d_Prom_EnhA <- sqrt((tAfterAPC$XEnhA-tAfterAPC$Xprom)**2+(tAfterAPC$YEnhA-tAfterAPC$Yprom)**2)
tAfterAPC$d_EnhA_EnhB <- sqrt((tAfterAPC$XEnhb-tAfterAPC$XEnhA)**2+(tAfterAPC$YEnhb-tAfterAPC$YEnhA)**2)
tAfterAPC$d_Prom_Enhb <- sqrt((tAfterAPC$XEnhb-tAfterAPC$Xprom)**2+(tAfterAPC$YEnhb-tAfterAPC$Yprom)**2)

tAfterAPC_NPC_WT <- tAfterAPC
NPC_WT_Tricontact <- tAfterAPC_NPC_WT[tAfterAPC_NPC_WT$d_Prom_EnhA< My_Threshold &
                                        tAfterAPC_NPC_WT$d_EnhA_EnhB<My_Threshold &
                                        tAfterAPC_NPC_WT$d_Prom_Enhb<My_Threshold, 1]
NPC_WT_Tricontact <- NPC_WT_AfterAPC[NPC_WT_AfterAPC$ID %in% NPC_WT_Tricontact,]


#For NPC_dA
NPC_dA <- rbind(read.csv2("Merge_summary/All_replicates/NPC_dA_rep1_merge.csv", header = T, stringsAsFactors = F),
                read.csv2("Merge_summary/All_replicates/NPC_dA_rep2_merge.csv", header = T, stringsAsFactors = F),
                read.csv2("Merge_summary/All_replicates/NPC_dA_rep3_merge.csv", header = T, stringsAsFactors = F))
NPC_dA[,-1]<- apply(NPC_dA[,-1],2, as.numeric)

#Calculate distances then remove dist>2.2µm for NPC dA
NPC_dA$F1F2 <- sqrt(((NPC_dA$xFC2-NPC_dA$xFC1)*40)**2+((NPC_dA$yFC2-NPC_dA$yFC1)*40)**2+((NPC_dA$zFC2-NPC_dA$zFC1)*140)**2)
NPC_dA$F2F3 <- sqrt(((NPC_dA$xFC3-NPC_dA$xFC2)*40)**2+((NPC_dA$yFC3-NPC_dA$yFC2)*40)**2+((NPC_dA$zFC3-NPC_dA$zFC2)*140)**2)
NPC_dA$F1F3 <- sqrt(((NPC_dA$xFC3-NPC_dA$xFC1)*40)**2+((NPC_dA$yFC3-NPC_dA$yFC1)*40)**2+((NPC_dA$zFC3-NPC_dA$zFC1)*140)**2)
NPC_dA <- cut_off_mydist(NPC_dA)

tAfterAPC <- data.frame(ID=rep(NA,NROW(NPC_dA)),
                        Xprom=rep(NA,NROW(NPC_dA)),
                        Yprom=rep(NA,NROW(NPC_dA)),
                        XEnhA=rep(NA,NROW(NPC_dA)),
                        YEnhA=rep(NA,NROW(NPC_dA)),
                        XEnhb=rep(NA,NROW(NPC_dA)),
                        YEnhb=rep(NA,NROW(NPC_dA)))
u <- 1
for (i in seq(1,NROW(NPC_dA_AfterAPC),by=3)) {
  tAfterAPC[u,1] <- NPC_dA_AfterAPC[i,1]
  tAfterAPC[u,2:3] <- NPC_dA_AfterAPC[i,2:3]
  tAfterAPC[u,4:5] <- NPC_dA_AfterAPC[i+1,2:3]
  tAfterAPC[u,6:7] <- NPC_dA_AfterAPC[i+2,2:3]
  u <- u+1
}

tAfterAPC$d_Prom_EnhA <- sqrt((tAfterAPC$XEnhA-tAfterAPC$Xprom)**2+(tAfterAPC$YEnhA-tAfterAPC$Yprom)**2)
tAfterAPC$d_EnhA_EnhB <- sqrt((tAfterAPC$XEnhb-tAfterAPC$XEnhA)**2+(tAfterAPC$YEnhb-tAfterAPC$YEnhA)**2)
tAfterAPC$d_Prom_Enhb <- sqrt((tAfterAPC$XEnhb-tAfterAPC$Xprom)**2+(tAfterAPC$YEnhb-tAfterAPC$Yprom)**2)

#Keep the same number points for both conditions, same as WT
NPC_dA_AfterAPC_1452 <- NPC_dA_AfterAPC[1:1452,]
tAfterAPC_484 <- tAfterAPC[1:484,]

NPC_dA_Tricontact_1452 <- tAfterAPC_484[tAfterAPC_484$d_Prom_EnhA< My_Threshold &
                                          tAfterAPC_484$d_EnhA_EnhB<My_Threshold &
                                          tAfterAPC_484$d_Prom_Enhb<My_Threshold, 1]
NPC_dA_Tricontact_1452 <- NPC_dA_AfterAPC_1452[NPC_dA_AfterAPC_1452$ID %in% NPC_dA_Tricontact_1452,]


#Calculate the area of tricontacts for WT and dA
area_tricontact_NPC <- tAfterAPC_NPC_WT[tAfterAPC_NPC_WT$ID %in% unique(NPC_WT_Tricontact$ID),]
area_tricontact_dA <-  tAfterAPC_484[tAfterAPC_484$ID %in% unique(NPC_dA_Tricontact_1452$ID),]

#Perimeter
area_tricontact_NPC$perimeter <- (area_tricontact_NPC$d_Prom_EnhA+ area_tricontact_NPC$d_EnhA_EnhB+ area_tricontact_NPC$d_Prom_Enhb)
area_tricontact_dA$perimeter <- (area_tricontact_dA$d_Prom_EnhA+ area_tricontact_dA$d_EnhA_EnhB+ area_tricontact_dA$d_Prom_Enhb)

#semi perimeter
area_tricontact_NPC$semiperimeter <- (area_tricontact_NPC$perimeter/2)
area_tricontact_dA$semiperimeter <- (area_tricontact_dA$perimeter/2)

#area
area_tricontact_NPC$area <- sqrt(area_tricontact_NPC$semiperimeter*(area_tricontact_NPC$semiperimeter - area_tricontact_NPC$d_Prom_EnhA)*(area_tricontact_NPC$semiperimeter - area_tricontact_NPC$d_EnhA_EnhB)*(area_tricontact_NPC$semiperimeter - area_tricontact_NPC$d_Prom_Enhb))
area_tricontact_dA$area <- sqrt(area_tricontact_dA$semiperimeter*(area_tricontact_dA$semiperimeter - area_tricontact_dA$d_Prom_EnhA)*(area_tricontact_dA$semiperimeter - area_tricontact_dA$d_EnhA_EnhB)*(area_tricontact_dA$semiperimeter - area_tricontact_dA$d_Prom_Enhb))

area_tricontact_NPC$condition <- "_WT"
area_tricontact_dA$condition <- "dA"
All_area_tricontact <- rbind(area_tricontact_NPC, area_tricontact_dA)
All_area_tricontact$area <- All_area_tricontact$area/1000

pdf(paste0("R_outpout/area_tricontact_T",My_Threshold,".pdf"), width = 5, height = 5)
ggviolin(All_area_tricontact, x="condition", y="area", fill="condition", palette = "Dark2",
         add=c("boxplot"), add.params=list(fill="white", width=0.1),
         width = 0.5, order=c("_WT", "dA"), title = "Tricontact t=400nm")+
  stat_compare_means( method = "wilcox.test", label.y = 75)+
  ylab("Area (µm2)")
dev.off()
