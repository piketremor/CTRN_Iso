#load packages
library(MEForLab)
library(dplyr)

#load data
setwd("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv")
trees<-read.csv("Trees2023.csv")
locs<-read.csv("Tree_locations_species.csv")
KI1_dist<-read.csv("KI_stemMap_Plot1.csv")
KI5_dist<-read.csv("KI_stemMap_Plot5.csv")

#start with KI plot 1
KItree<-filter(trees, SITEid=="KI",PLOTid=="1")
KI18<-filter(KItree,YEAR==2018)
names(KI18)[names(KI18)=="TREE"]<-"tree.x"
KIx<-KI18[,c(4,6)]
names(KIx)[names(KIx)=="DBH"]<-"DBH.x"

KI2<-left_join(KI1_dist, KIx)

##j
names(KIx)[names(KIx)=="tree.x"]<-"tree.j"
names(KIx)[names(KIx)=="DBH.x"]<-"DBH.j"

KI1<-left_join(KI2, KIx) #change name of dataframe to correspond to the correct plot

KI1["comp.heg"] <- 
  mapply(hegyi.index, DBH=KI1$DBH.x, dbhi=KI1$DBH.j, distance=KI1$distance)

##repeat for KI 5
KItree<-filter(trees, SITEid=="KI",PLOTid=="5")
KI18<-filter(KItree,YEAR==2018)
names(KI18)[names(KI18)=="TREE"]<-"tree.x"
KIx<-KI18[,c(4,6)]
names(KIx)[names(KIx)=="DBH"]<-"DBH.x"
names(KI5_dist)[names(KI5_dist)=="focal.tree"]<-"tree.x"
names(KI5_dist)[names(KI5_dist)=="competitor.tree"]<-"tree.j"

KI2<-left_join(KI5_dist, KIx)

##j
names(KIx)[names(KIx)=="tree.x"]<-"tree.j"
names(KIx)[names(KIx)=="DBH.x"]<-"DBH.j"

KI5<-left_join(KI2, KIx) #change name of dataframe to correspond to the correct plot

KI5["comp.heg"] <- 
  mapply(hegyi.index, DBH=KI5$DBH.x, dbhi=KI5$DBH.j, distance=KI5$distance)
KI5<-KI5[,2:9] #remove random x variable


