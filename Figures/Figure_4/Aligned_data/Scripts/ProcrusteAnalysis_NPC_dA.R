rm(list=ls())
library(paleomorph)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(reshape2)

# Import and select the coordinates of triangles made by Prom/EnhA/EnhB
setwd("../Desktop/Temp from Archive/")
NPC_dA <- rbind(read.csv2("Merge_summary/All_replicates/NPC_dA_rep1_merge.csv", header = T, stringsAsFactors = F),
                read.csv2("Merge_summary/All_replicates/NPC_dA_rep2_merge.csv", header = T, stringsAsFactors = F),
                read.csv2("Merge_summary/All_replicates/NPC_dA_rep3_merge.csv", header = T, stringsAsFactors = F))
NPC_dA[,-1]<- apply(NPC_dA[,-1],2, as.numeric)
#str(NPC_dA)
#Calculate distances then remove dist>2.2µm for NPC dA
#dist FC1FC2
NPC_dA$F1F2 <- sqrt(((NPC_dA$xFC2-NPC_dA$xFC1)*40)**2+((NPC_dA$yFC2-NPC_dA$yFC1)*40)**2+((NPC_dA$zFC2-NPC_dA$zFC1)*140)**2)
#dist FC2FC3
NPC_dA$F2F3 <- sqrt(((NPC_dA$xFC3-NPC_dA$xFC2)*40)**2+((NPC_dA$yFC3-NPC_dA$yFC2)*40)**2+((NPC_dA$zFC3-NPC_dA$zFC2)*140)**2)
#dist FC1FC3
NPC_dA$F1F3 <- sqrt(((NPC_dA$xFC3-NPC_dA$xFC1)*40)**2+((NPC_dA$yFC3-NPC_dA$yFC1)*40)**2+((NPC_dA$zFC3-NPC_dA$zFC1)*140)**2)
cut_off_mydist<- function(a){
  a <- a[a$F1F2 %in% subset(a$F1F2, a$F1F2 < 2200),]
  a <- a[a$F2F3 %in% subset(a$F2F3, a$F2F3 < 2200),]
  a <- a[a$F1F3 %in% subset(a$F1F3, a$F1F3 < 2200),]
}
NPC_dA <- cut_off_mydist(NPC_dA)

#Select the data corresponding to the coordinates of the 3 different points
B <- NPC_dA
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
ma[,,1]
B[1,]
# Align the data 
aligned <- procrustes(ma, scale=F)

plotSpecimens(aligned,size= 5)


##Create a new data frame with the coordinate of the 3 points after 
##procrutres superimposition
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

eig.val <- get_eigenvalue(res.pca)
eig.val
fviz_eig(res.pca, addlabels = TRUE)
fviz_pca_ind (res.pca)
fviz_pca_ind(res.pca,
             geom.ind = "point", # Montre les points seulement (mais pas le "text")
             col.ind = PreAPC$molecular_ID, # colorer by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Ellipses de concentration
             legend.title = "Groups"
)
#To get the coordinates
ind <- get_pca_ind(res.pca)
ind
head(ind$coord)
dim(ind$coord)
str(ind)
AfterAPC <- data.frame(
  ID = PreAPC$ID,
  Dim1 = ind$coord[, 1],
  Dim2 = ind$coord[, 2],
  molecular_ID = PreAPC$molecular_ID
)
write.csv2(AfterAPC, "R_outpout/NPC_dA/NPC_dA_AfterAPC.csv",row.names = F)


#First graph with the 2D coordinates after APC
#pdf("R_outpout/2Dcoorrdintes_afterAPC.pdf", width = 8, height = 8)
ggplot(AfterAPC, aes(x=Dim1, y=Dim2, group=ID))+
  geom_polygon(aes(col=ID), fill=NA, colour="grey51")+
  geom_point(aes(col=molecular_ID), size=2, alpha=0.7)+
  #geom_polygon(aes(col=ID), fill=NA, colour="black")
  #stat_density_2d(aes(group=molecular_ID), geom = "polygon", alpha=0.5)
  theme_minimal()
#dev.off()
#Zoom btw -500 and 500 of the first graph

ggplot(AfterAPC, aes(x=Dim1, y=Dim2, group=ID))+
  geom_polygon(aes(col=ID), fill=NA, colour="grey51")+
  geom_point(aes(col=molecular_ID), size=2, alpha=0.7)+
  #geom_polygon(aes(col=ID), fill=NA, colour="black")
  #stat_density_2d(aes(group=molecular_ID), geom = "polygon", alpha=0.5)
  theme_minimal()+
  xlim(-500,500)+ylim(-500,500)+
  xlab("Distance (nm)")+ylab("Distance (nm)")

ggplot(AfterAPC[1:321,], aes(x=Dim1, y=Dim2, group=ID))+
  geom_polygon(aes(col=ID), fill=NA, colour="grey51")+
  geom_point(aes(col=molecular_ID), size=2, alpha=0.7)

#pdf("Density_2Dcoordinates_aFterAPC.pdf",width = 8,height = 8)  
ggplot(AfterAPC, aes(x=Dim1, y=Dim2))+
  geom_point(aes(col=molecular_ID), size=2)+
  stat_density2d(
    data=AfterAPC, aes( x=Dim1,y=Dim2, fill=molecular_ID), alpha=0.2, geom="polygon", bins=15)
#dev.off()
ggplot(AfterAPC, aes(x=Dim1, y=Dim2))+
  geom_bin2d(bins=70)+
  scale_fill_continuous(type = "viridis")

ggplot(AfterAPC, aes(x=Dim1, y=Dim2, col=ID))+
  geom_point(aes(col=molecular_ID), size=2)+
  stat_density2d(
    data=AfterAPC, aes( x=Dim1,y=Dim2, col=molecular_ID), alpha=0.2, geom="polygon")+
  xlim(-500,500)+ylim(-500,500)


#Calculate the distances between the points after PCA : AfterAPC
#dist FC1FC2
#rm(tAfterAPC)
tAfterAPC <- data.frame(ID=rep(NA,NROW(NPC_dA)),
                        Xprom=rep(NA,NROW(NPC_dA)),
                        Yprom=rep(NA,NROW(NPC_dA)),
                        XEnhA=rep(NA,NROW(NPC_dA)),
                        YEnhA=rep(NA,NROW(NPC_dA)),
                        XEnhb=rep(NA,NROW(NPC_dA)),
                        YEnhb=rep(NA,NROW(NPC_dA)))
u <- 1
for (i in seq(1,NROW(AfterAPC),by=3)) {
  tAfterAPC[u,1] <- AfterAPC[i,1]
  tAfterAPC[u,2:3] <- AfterAPC[i,2:3]
  tAfterAPC[u,4:5] <- AfterAPC[i+1,2:3]
  tAfterAPC[u,6:7] <- AfterAPC[i+2,2:3]
  u <- u+1
  
}
#dist Prom EnhA
tAfterAPC$d_Prom_EnhA <- sqrt((tAfterAPC$XEnhA-tAfterAPC$Xprom)**2+(tAfterAPC$YEnhA-tAfterAPC$Yprom)**2)
#dist EnhA EnhB
tAfterAPC$d_EnhA_EnhB <- sqrt((tAfterAPC$XEnhb-tAfterAPC$XEnhA)**2+(tAfterAPC$YEnhb-tAfterAPC$YEnhA)**2)
#dist Prom EnhB
tAfterAPC$d_Prom_Enhb <- sqrt((tAfterAPC$XEnhb-tAfterAPC$Xprom)**2+(tAfterAPC$YEnhb-tAfterAPC$Yprom)**2)
summary(tAfterAPC[,8:10])
tAfterAPC$dd_EnhA_EnhB <- rep(NA, dim(tAfterAPC)[1])




##Collect coordinates x and y for distances shorter then threshold
#First threshold 


My_Threshold <- 250


close_distance_Prom_ENhA <- tAfterAPC[tAfterAPC$d_Prom_EnhA < My_Threshold,c(2,3,4,5)] #22
close_distance_Prom_ENhB <- tAfterAPC[tAfterAPC$d_Prom_Enhb < My_Threshold,c(2,3,6,7)]
close_distance_EnhA_Enhb <- tAfterAPC[tAfterAPC$d_EnhA_EnhB < My_Threshold,c(4,5,6,7)]
Tricontact <- tAfterAPC[tAfterAPC$d_Prom_EnhA< My_Threshold &
                          tAfterAPC$d_EnhA_EnhB<My_Threshold &
                          tAfterAPC$d_Prom_Enhb<My_Threshold, 1]
Tricontact <- AfterAPC[AfterAPC$ID %in% Tricontact,]

pdf("R_outpout/NPC_dA/Procrustre_SuperImposition_NPC_dA_2.2CutOff.pdf", height = 6, width = 7)
a <- ggplot(AfterAPC, aes(x=Dim1, y=Dim2))+
  geom_polygon(aes(group=ID), fill=NA, colour="grey51", size=0.01)+
  geom_point(aes(colour=molecular_ID), size=1)+
  #stat_density2d(
  #  data=AfterAPC, 
  #  aes(x=Dim1,y=Dim2, fill=molecular_ID), 
  #  alpha=0.2, geom="polygon", bins=15, show.legend=F)+
  # xlim(-1500,1000)+ylim(-1500,1000)+
  xlab("Distance (nm)")+ylab("Distance (nm)")+
  scale_colour_discrete(name="Molecular\nIdentity", labels= c("Enhancer A","Enhancer B","Promoter"))+
  theme(legend.position = c(0.9,0.2))
a
a+  geom_segment(data=close_distance_Prom_ENhA,
                 aes(x=close_distance_Prom_ENhA$Xprom, 
                     xend=close_distance_Prom_ENhA$XEnhA,
                     y= close_distance_Prom_ENhA$Yprom, 
                     yend=close_distance_Prom_ENhA$YEnhA),
                 col="purple")+
  
  geom_segment(data=close_distance_Prom_ENhB,
               aes(x=close_distance_Prom_ENhB$Xprom, 
                   xend=close_distance_Prom_ENhB$XEnhb,
                   y= close_distance_Prom_ENhB$Yprom, 
                   yend=close_distance_Prom_ENhB$YEnhb),
               col="cyan")+
  
  geom_segment(data=close_distance_EnhA_Enhb,
               aes(x=close_distance_EnhA_Enhb$XEnhA, 
                   xend=close_distance_EnhA_Enhb$XEnhb,
                   y= close_distance_EnhA_Enhb$YEnhA, 
                   yend=close_distance_EnhA_Enhb$YEnhb),
               col="yellow")+
  
  geom_polygon(data=Tricontact, aes(x=Dim1, y=Dim2, group=ID),  fill=NA, colour="red")  

##Zoom in
#colors <- c("Prom - Enh. A" = "purple", "Enh A. - Enh B." = "cyan", "Prom - Enh. B" = "yellow")
ggplot(AfterAPC, aes(x=Dim1, y=Dim2))+
  geom_polygon(aes(group=ID), fill=NA, colour="grey51")+
  geom_point(aes(col=molecular_ID), size=2, show.legend = F)+
  stat_density2d(
    data=AfterAPC, 
    aes( x=Dim1,y=Dim2, fill=molecular_ID), 
    alpha=0.2, geom="polygon", bins=15, show.legend = F)+
  xlim(-500,500)+ylim(-500,500)+
  geom_segment(data=close_distance_Prom_ENhA,
               aes(x=close_distance_Prom_ENhA$Xprom, 
                   xend=close_distance_Prom_ENhA$XEnhA,
                   y= close_distance_Prom_ENhA$Yprom, 
                   yend=close_distance_Prom_ENhA$YEnhA),
               color="purple", size=0.8)+
  
  geom_segment(data=close_distance_Prom_ENhB,
               aes(x=close_distance_Prom_ENhB$Xprom, 
                   xend=close_distance_Prom_ENhB$XEnhb,
                   y= close_distance_Prom_ENhB$Yprom, 
                   yend=close_distance_Prom_ENhB$YEnhb),
               col="cyan", size=0.8)+
  
  geom_segment(data=close_distance_EnhA_Enhb,
               aes(x=close_distance_EnhA_Enhb$XEnhA, 
                   xend=close_distance_EnhA_Enhb$XEnhb,
                   y= close_distance_EnhA_Enhb$YEnhA, 
                   yend=close_distance_EnhA_Enhb$YEnhb),
               col="yellow", size=0.8)+
  geom_polygon(data=Tricontact, aes(x=Dim1, y=Dim2, group=ID),  fill=NA, colour="red", size=0.8)+
  xlab("Distance (nm)")+ylab("Distance (nm)")+
  labs(color="",title = paste0("Threshold = ",My_Threshold,"nm"))
dev.off()