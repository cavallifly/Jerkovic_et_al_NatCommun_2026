rm(list=ls())
library(paleomorph)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(reshape2)

# Import and select the coordinates of triangles made by Prom/EnhA/EnhB
setwd("../Temp from Archive/")
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

#Select the data corresponding to the coordinates of the 3 different points
B <- NPC_WT
Bx <- B[, c(13, 16, 19)]
By <- B[, c(14, 17, 20)]
Bz <- B[, c(15, 18, 21)]
B <- cbind(Bx, By, Bz)
B <- B[!B$xFC1=="NaN",]
B[,1:3] <- B[,1:3]*40
B[,4:6] <- B[,4:6]*40
B[,7:9] <- B[,7:9]*140

#Generate array with the coordinates for Procrustes Superimposition
num.coords <- dim(B)[2]
coor3D <- 3
ns <- dim(B)[1]
ma <- array(dim=c(num.coords/coor3D, coor3D, ns))
xi <- c(1,2,3); yi <- xi+3; zi <- yi+3
for (i in 1:ns) {
  ma[, 1, i] <- unlist(B[i, xi])
  ma[, 2, i] <- unlist(B[i, yi])
  ma[, 3, i] <- unlist(B[i, zi])
}

# Align the data
aligned <- procrustes(ma, scale=F)

#Create a new data frame with the coordinates of the 3 points after Procrustes superimposition
molecular_ID = c("Prom", "EnhA", "EnhB")
PreAPC <- data.frame(
  ID = rep(NA, dim(aligned)[3]*3),
  X = rep(NA, dim(aligned)[3]*3),
  Y = rep(NA, dim(aligned)[3]*3),
  Z = rep(NA, dim(aligned)[3]*3),
  molecular_ID = rep(molecular_ID,dim(aligned)[3])
)
u <- 1; v <- 2; w <- 3
for (i in 1:dim(aligned)[3]) {
  PreAPC[u:w,1] <- i
  PreAPC[u,2:4] <- aligned[1,1:3,i]
  PreAPC[v,2:4] <- aligned[2,1:3,i]
  PreAPC[w,2:4] <- aligned[3,1:3,i]
  PreAPC[u:w,5] <- molecular_ID
  u <- u+3; v <- v+3; w <- w+3
}

res.pca <- PCA(PreAPC[,2:4],scale.unit = F, ncp=2, graph = FALSE)
ind <- get_pca_ind(res.pca)
AfterAPC <- data.frame(
  ID = PreAPC$ID,
  Dim1 = ind$coord[, 1],
  Dim2 = ind$coord[, 2],
  molecular_ID = PreAPC$molecular_ID
)

#Calculate the distances between the points after PCA
tAfterAPC <- data.frame(ID=rep(NA,NROW(NPC_WT)),
                        Xprom=rep(NA,NROW(NPC_WT)),
                        Yprom=rep(NA,NROW(NPC_WT)),
                        XEnhA=rep(NA,NROW(NPC_WT)),
                        YEnhA=rep(NA,NROW(NPC_WT)),
                        XEnhb=rep(NA,NROW(NPC_WT)),
                        YEnhb=rep(NA,NROW(NPC_WT)))
u <- 1
for (i in seq(1,NROW(AfterAPC),by=3)) {
  tAfterAPC[u,1] <- AfterAPC[i,1]
  tAfterAPC[u,2:3] <- AfterAPC[i,2:3]
  tAfterAPC[u,4:5] <- AfterAPC[i+1,2:3]
  tAfterAPC[u,6:7] <- AfterAPC[i+2,2:3]
  u <- u+1
}

tAfterAPC$d_Prom_EnhA <- sqrt((tAfterAPC$XEnhA-tAfterAPC$Xprom)**2+(tAfterAPC$YEnhA-tAfterAPC$Yprom)**2)
tAfterAPC$d_EnhA_EnhB <- sqrt((tAfterAPC$XEnhb-tAfterAPC$XEnhA)**2+(tAfterAPC$YEnhb-tAfterAPC$YEnhA)**2)
tAfterAPC$d_Prom_Enhb <- sqrt((tAfterAPC$XEnhb-tAfterAPC$Xprom)**2+(tAfterAPC$YEnhb-tAfterAPC$Yprom)**2)

#Calculate the cumulative number of bi- and tri-contacts for each threshold
q <- seq(0,2200,by=50)
threshold <- data.frame(threshold=rep(NA,NROW(q)),
                        tricontact=rep(NA,NROW(q)),
                        d_Prom_EnhA=rep(NA,NROW(q)),
                        d_EnhA_EnhB=rep(NA,NROW(q)),
                        d_Prom_Enhb=rep(NA,NROW(q)))
u <- 1
for (i in q) {
  threshold[u,1] <- i
  threshold[u,2]<- NROW(tAfterAPC[tAfterAPC$d_Prom_EnhA < i &
                                    tAfterAPC$d_EnhA_EnhB < i &
                                    tAfterAPC$d_Prom_Enhb < i,])
  threshold[u,3] <- NROW(tAfterAPC[tAfterAPC$d_Prom_EnhA < i,])
  threshold[u,4]<- NROW(tAfterAPC[tAfterAPC$d_EnhA_EnhB < i,])
  threshold[u,5] <- NROW(tAfterAPC[tAfterAPC$d_Prom_Enhb < i,])
  u <- u + 1
}

#Normalize the contact frequencies
thresholdnorm_noTricontact <- data.frame(threshold=threshold$threshold,
                                         tricontact=threshold$tricontact,
                                         d_Prom_EnhA=threshold$d_Prom_EnhA,
                                         d_EnhA_EnhB=threshold$d_EnhA_EnhB,
                                         d_Prom_Enhb=threshold$d_Prom_Enhb)

for (i in 1:NROW(thresholdnorm_noTricontact)) {
  thresholdnorm_noTricontact[i,3:5] <- threshold[i,3:5]/107
}
thresholdnorm_noTricontact$Observed_TriContact <- thresholdnorm_noTricontact$tricontact/107

t_thresholdnorm_noTricontact <- melt(
  thresholdnorm_noTricontact[,c(1,3,4,5,6)],
  id.vars = "threshold",
  measure.vars=c("d_Prom_EnhA", "d_EnhA_EnhB", "d_Prom_Enhb", "Observed_TriContact")
)

pdf("Cumulative_frequency_Bi_Tri_Contact_from_0_to_2200nm.pdf")
ggplot(t_thresholdnorm_noTricontact, aes(x=threshold, y=value, colour=variable))+
  geom_point()+
  geom_smooth(se=F, span=0.2)+
  ylab("Bi / Tri-contact frequency")+xlab("Threshold (nm)")+
  scale_color_discrete(name="",
                       labels=c("Promoter - Enh. A",
                                "Enh. A - Enh. B ",
                                "Promoter - Enh. B",
                                "Tri-Contact") )
dev.off()

pdf("Cumulative_frequency_Bi_Tri_Contact_from_0_to_600nm.pdf")
ggplot(t_thresholdnorm_noTricontact[t_thresholdnorm_noTricontact$threshold<650,], aes(x=threshold, y=value, colour=variable))+
  geom_point()+
  geom_smooth(se=F)+
  ylab("Bi / Tri-contact frequency")+xlab("Threshold (nm)")+
  scale_color_discrete(name="",
                       labels=c("Promoter - Enh. A",
                                "Enh. A - Enh. B ",
                                "Promoter - Enh. B",
                                "Tri-Contact") )
dev.off()
