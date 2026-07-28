rm(list=ls())
library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(splancs)
library(MASS)
library(fields)

setwd("C://Users/pc/Documents/Thèse/Temp from Archive/")


#Import AfterAPC for WT and dA
NPC_WT_AfterAPC <- read.csv2("R_outpout/NPC_WT/NPC_WT_AfterAPC.csv", stringsAsFactors = F)
NPC_dA_AfterAPC <- read.csv2("R_outpout/NPC_dA//NPC_dA_AfterAPC.csv", stringsAsFactors = F)
NPC_dA_AfterAPC[,2] <- NPC_dA_AfterAPC[,2]*-1
NPC_dA_AfterAPC[,3] <- NPC_dA_AfterAPC[,3]*-1
NPC_dA_AfterAPC <- NPC_dA_AfterAPC[1:1452,]


#####################################
####    With 2d Kernel Density   ####
#####################################

Area_kd <- data.frame(step=rep(NA,60),
                      Area.WT=rep(NA,60),
                      Area.dA=rep(NA,60),
                      loci=rep(NA,60))
j=1

for (loci in unique(NPC_WT_AfterAPC$molecular_ID)) {
  #for WT
  x= NPC_WT_AfterAPC[NPC_WT_AfterAPC$molecular_ID==loci,2]
  y= NPC_WT_AfterAPC[NPC_WT_AfterAPC$molecular_ID==loci,3]
  f1 <- kde2d(x=x, y=y,n=20, lims = c(range(x), range(y)))
  z=interp.surface(f1,loc = matrix(c(x,y), ncol=2))
  df <- data.frame(x=x, y=y, z=z)
  df$loci <- loci
  u=j
  for (i in seq(from=0.05, to=1, by=0.05)) {
    order <- df[order(df$z, decreasing=T),]
    order <- order[1:(484*i),]
    hull <- order %>% slice(chull(x, y))
    Area_kd[u,1] <- i
    Area_kd[u,2] <- as.matrix(hull[,1:2]) %>% areapl
    Area_kd[u,4] <- loci
    u=u+1
  }

  #for dA
  x= NPC_dA_AfterAPC[NPC_dA_AfterAPC$molecular_ID==loci,2]
  y= NPC_dA_AfterAPC[NPC_dA_AfterAPC$molecular_ID==loci,3]
  f1 <- kde2d(x=x, y=y,n=20, lims = c(range(x), range(y)))
  z=interp.surface(f1,loc = matrix(c(x,y), ncol=2))
  df <- data.frame(x=x, y=y,z=z)
  df$loci <- loci
  u=j
  for (i in seq(from=0.05, to=1, by=0.05)) {
    order <- df[order(df$z, decreasing=T),]
    order <- order[1:(484*i),]
    hull <- order %>% slice(chull(x, y))
    Area_kd[u,1] <- i
    Area_kd[u,3] <- as.matrix(hull[,1:2]) %>% areapl
    u=u+1
  }
  j=u
}

tArea_kd <- melt(Area_kd, id.vars = c("step", "loci"))


#Compare Area in WT and dA for the 3 loci
for (loci in unique(NPC_WT_AfterAPC$molecular_ID)) {
  a <- ggpaired(
    tArea_kd[tArea_kd$loci==loci,], x = "variable", y = "value", color = "variable", palette = "jco",
    line.color = "gray", line.size = 0.4)+
    yscale("log10")+
    stat_compare_means(paired = T)+
    ylab("Area nm2")+
    xlab("")+
    ggtitle(loci)+
    theme(legend.position='none')
  print(a)
  pdf(paste0(loci,"_boxplot_area.pdf"), height = 4, width = 3 )
  print(a)
  dev.off()
}
