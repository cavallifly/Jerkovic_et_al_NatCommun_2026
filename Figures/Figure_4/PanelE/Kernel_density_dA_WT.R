rm(list=ls())
library(ggplot2)
library(gridExtra)
library(MASS)
library(fields)
setwd("C://Users/pc/Documents/Thèse/Temp from Archive/")


#######Prepare the df#########
#Import AfterAPC for WT and dA
NPC_WT_AfterAPC <- read.csv2("R_outpout/NPC_WT/NPC_WT_AfterAPC.csv", stringsAsFactors = F)
NPC_dA_AfterAPC <- read.csv2("R_outpout/NPC_dA//NPC_dA_AfterAPC.csv", stringsAsFactors = F)
NPC_dA_AfterAPC[,2] <- NPC_dA_AfterAPC[,2]*-1
NPC_dA_AfterAPC[,3] <- NPC_dA_AfterAPC[,3]*-1
NPC_dA_AfterAPC <- NPC_dA_AfterAPC[1:1452,]


#####################################
####    With 2d Kernel Density   ####
#####################################
df_wt <- data.frame(
  x = rep(NA, 60),
  y = rep(NA, 60),
  z = rep(NA, 60),
  loci = rep(NA, 60))
df_dA <- data.frame(
  x = rep(NA, 60),
  y = rep(NA, 60),
  z = rep(NA, 60),
  loci = rep(NA, 60))

for (loci in unique(NPC_WT_AfterAPC$molecular_ID)) {
  #for WT
  x= NPC_WT_AfterAPC[NPC_WT_AfterAPC$molecular_ID==loci,2]
  y= NPC_WT_AfterAPC[NPC_WT_AfterAPC$molecular_ID==loci,3]
  f1 <- kde2d(x=x, y=y,n=20, lims = c(range(x), range(y)))
  z=interp.surface(f1,loc = matrix(c(x,y), ncol=2))
  df <- data.frame(x=x, y=y, z=z)
  df$loci <- loci
  df_wt <- rbind(df_wt, df)
  
  #for dA
  x= NPC_dA_AfterAPC[NPC_dA_AfterAPC$molecular_ID==loci,2]
  y= NPC_dA_AfterAPC[NPC_dA_AfterAPC$molecular_ID==loci,3]
  f1 <- kde2d(x=x, y=y,n=20, lims = c(range(x), range(y)))
  z=interp.surface(f1,loc = matrix(c(x,y), ncol=2))
  df <- data.frame(x=x, y=y,z=z)
  df$loci <- loci
  df_dA <- rbind(df_dA, df)
}

a1 <-
  ggplot(df_wt, aes(x=x, y=y, colour=z))+
  geom_point(size=0.2 )+
  scale_color_gradientn(colours = rainbow(20), limits=c(min(df_dA$z, na.rm = T), max(df_wt$z, na.rm = T )))+
  ylim(-1250,1000)+xlim(-1500, 1250)+
  theme(legend.position='none')+
  ggtitle("wt")

a2 <- 
  ggplot(df_dA, aes(x=x, y=y, colour=z))+
  geom_point(size=0.2)+
  scale_color_gradientn(colours = rainbow(20), limits=c(min(df_dA$z, na.rm = T), max(df_wt$z, na.rm = T )))+
  ylim(-1250,1000)+xlim(-1500, 1250)+
  ggtitle("dA")

pdf("kernel_density_WT_dA.pdf", height = 5, width = 8)
grid.arrange(a1, a2, ncol=2, widths= c(1.5,2))
dev.off()
