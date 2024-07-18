#load packages
library(MEForLab)
library(dplyr)
library(ggplot2)

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

names(KI1_dist)[names(KI1_dist)=="focal.tree"]<-"tree.x"
names(KI1_dist)[names(KI1_dist)=="competitor.tree"]<-"tree.j"

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

#Rump Road competition calculations
RR1_dist<-read.csv("RR_stemMap_Plot1.csv")
RR3_dist<-read.csv("RR_stemMap_Plot3.csv")


#RR1
RRtree<-filter(trees, SITEid=="RR",PLOTid=="1")
RR03<-filter(RRtree,YEAR==2003)
names(RR03)[names(RR03)=="TREE"]<-"tree.x"
RRx<-RR03[,c(4,6)]
names(RRx)[names(RRx)=="DBH"]<-"DBH.x"

names(RR1_dist)[names(RR1_dist)=="focal.tree"]<-"tree.x"
names(RR1_dist)[names(RR1_dist)=="competitor.tree"]<-"tree.j"

RR2<-left_join(RR1_dist, RRx)

##j
names(RRx)[names(RRx)=="tree.x"]<-"tree.j"
names(RRx)[names(RRx)=="DBH.x"]<-"DBH.j"

RR1<-left_join(RR2, RRx) #change name of dataframe to correspond to the correct plot

RR1["comp.heg"] <- 
  mapply(hegyi.index, DBH=RR1$DBH.x, dbhi=RR1$DBH.j, distance=RR1$distance)

##RR3
RRtree<-filter(trees, SITEid=="RR",PLOTid=="3")
RR03<-filter(RRtree,YEAR==2003)
names(RR03)[names(RR03)=="TREE"]<-"tree.x"
RRx<-RR03[,c(4,6)]
names(RRx)[names(RRx)=="DBH"]<-"DBH.x"

names(RR1_dist)[names(RR1_dist)=="focal.tree"]<-"tree.x"
names(RR1_dist)[names(RR1_dist)=="competitor.tree"]<-"tree.j"

RR2<-left_join(RR1_dist, RRx)

##j
names(RRx)[names(RRx)=="tree.x"]<-"tree.j"
names(RRx)[names(RRx)=="DBH.x"]<-"DBH.j"

RR3<-left_join(RR2, RRx) #change name of dataframe to correspond to the correct plot

RR3["comp.heg"] <- 
  mapply(hegyi.index, DBH=RR1$DBH.x, dbhi=RR1$DBH.j, distance=RR1$distance)


#KI_COF<-read.csv("~/Google Drive/My Drive/Dendrochronology/KIfull_age.csv")
#names(KI_COF)[names(KI_COF)=="Plot"]<-"PLOTid"
#names(KI_COF)[names(KI_COF)=="Tree"]<-"tree.x"
#KI_comp<-left_join(KI_COF, KI5)

#test<-KI_comp%>%
 # group_by(tree.x)

#lmKI5<-lm(sens~comp.heg,data=KI_comp)
#summary(lmKI5)
#plot(lmKI5)

#KI5 control
#ggplot(KI_comp, aes(x=comp.heg,y=sens,color=SPP))+
 # geom_point()+
  #geom_smooth(method=lm)+
  #theme_classic()
 


#KI1 thinned
KI1_comp<-left_join(KI_COF, KI1)
lmKI1<-lm(sens~comp.heg,data=KI1_comp)
summary(lmKI1)

ggplot(KI1_comp, aes(x=comp.heg,y=sens,color=SPP))+
  geom_point()+
  geom_smooth(method=lm)+
  theme_classic()
  
