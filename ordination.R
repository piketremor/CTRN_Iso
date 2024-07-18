#load packages
dev.off()
rm(list=ls())
library(dplyr)
library(mosaic)
library(forcats)
library(tidyverse)
library(nlme)
library(ggplot2)
library(leaps)
library(ggfortify)
library(ggeffects)
library(ggeasy)
library(reshape2)
library(vegan)
library(psych)
library(reshape)
library(janitor)

setwd("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv")
saplings <- read.csv("Saplings.csv")[,2:9]

saplings <- filter(saplings, SITEid == "AS" | SITEid == "DR" | SITEid == "GR" | SITEid == "HR" | SITEid == "KI" | SITEid == "LM" | SITEid == "LT" | SITEid == "PA" | SITEid == "PE" | SITEid == "RC" | SITEid == "RR" | SITEid == "SA" | SITEid == "SC" | SITEid == "SR" | SITEid == "WB") 
saplings[saplings == "SpecAld"]<-"SA"
saplings[saplings == "HM"]<-"EH"
saplings[saplings == "CH"]<-"BC"


##2018 Ordination and Data Cleaning


#Filter for just year 2018
saplings18<-filter(saplings, YEAR == 2018)

#filter for species >5%
branch<-filter(saplings18, SPP =="PB"|SPP == "OT"|SPP=="NC"|SPP=="RM"|SPP=="RS"|SPP=="WP"|SPP=="BF"|SPP=="YB"|SPP=="EH"|SPP=="QA"|SPP=="BC")

#recalculate iv values
branch <- branch%>%
  mutate(ba.half = (0.5^2*0.005454)*X1.2.inch)%>%
  mutate(ba.one = (1.0^2*0.005454)*X1.inch)%>%
  mutate(ba.two = (2.0^2)*0.005454*X2.inch)%>%
  mutate(sap.ba = (ba.half+ba.one+ba.two)*250)

branch$SAP_EXP <- 250

branch<-replace(branch, is.na(branch), 0)   


branch<-branch%>%
  mutate(total_diameter = ba.half + ba.one + ba.two)

branch<-branch%>%
  mutate(BA = ((total_diameter^2)*0.005454))

branch<-branch%>%
  mutate(SAP_BA_TPA = BA * SAP_EXP)

plot_summary_sap <- branch%>%
  group_by(SITEid, PLOTid, YEAR)%>%
  summarise(SAP_TPA_total = sum(SAP_EXP),
            SAPBA_total = sum(SAP_BA_TPA))

branch <- branch%>%
  group_by(SITEid, PLOTid, YEAR, SPP)%>%
  summarise(SAP_TPA = sum(SAP_EXP),
            SAPBA = sum(SAP_BA_TPA))

branch<-left_join(branch, plot_summary_sap)


branch <- branch%>%
  mutate(prop_tpa = (SAP_TPA/SAP_TPA_total),
         prop_ba = (SAPBA/SAPBA_total))

branch <- branch%>%
  mutate(iv = ((prop_tpa + prop_ba)/2))

#create site-plot identifier
branch$sapID<-paste0(branch$SITEid,"-",branch$PLOTid)

#get data in long format
molten <- melt(as.data.frame(branch),id=c("sapID","iv","SPP"))
sapwide <- dcast(molten,sapID~SPP,value.var = "iv",mean)
sapwide[is.na(sapwide)] <- 0
head(sapwide)
# double check they all add to 1
rowSums(sapwide[2:12])


#NMDS
sapordi<-metaMDS(sapwide[,2:12], distance = "bray")
stressplot(sapordi)
plot(sapordi)

plot(sapordi,type="n")
points(sapordi,display="sites",cex=2,pch=21,col="red", bg="yellow")
text(sapordi,display="spec",cex=1.5,col="blue")



## Constrained Ordination###
#clean the environmental data
env<-env[,c(2,3,4,27:65)]
env$sapID<-paste0(env$SITEid,"-",env$PLOTid)
env<-filter(env, YEAR==2018)
env<-env[,c(4:43)]
env<-env[,c(40, 1:39)]
env[is.na(env)] <- 0

#change from integer to numeric
env$flowdir<-as.numeric(env$flowdir)
env$Parent<-as.numeric(env$Parent)

#join together sapling and environmental data to make sure the number of observations match
bark<-left_join(env, sapwide)
bark[is.na(bark)] <- 0
names<-bark$sapID #create list of sapID to set the rownames with 

#separate the datasets back out
sapwide<-bark[c(1, 41:51)]
env<-bark[,1:40]

#add rownames
rownames(sapwide)<-names
rownames(env)<-names

#remove sapID variable
sapwide<-sapwide[,2:12]
env<-env[,2:40]

#rda
sap.rda <- rda(sapwide ~ ., data=env)
summary(sap.rda)
ordiplot(sap.rda, scaling = 2, type = "text")
step.forward <- ordiR2step(rda(sapwide ~ 1, data=env), scope=formula(sap.rda), R2scope = F, direction="forward", pstep=1000)
anova(sap.rda, by="terms", step=1000)

#reduced model
sap.rda2<-rda(sapwide ~ wd+SWI+ppt+SWC2+WD2020+Winds10+vpdmax+tri+dep+vpdmin, data=env)
summary(sap.rda2)
ordiplot(sap.rda2, scaling = 2, type = "text")











#calculate which species are below 5%
#sapwide <- dcast(saplings18,sapID~SPP, value.var = "X1.2.inch")

#sapwide<-as.data.frame(sapwide)

#sapwide[,2:21][sapwide[,2:21] > 0]<-1


#sapwide<-sapwide%>%
#  adorn_totals("row")


#RO, WS,AB, GB, SA, SM, ST, YR species <5%
#replace HM with EH


#nmds
#set.seed(1234)
#nmdssappy<-metaMDS(sapwide[,2:21], type="t")
#stressplot(nmdssappy)

#plot(nmdssappy)
#data.scores <- as.data.frame(scores(nmdssappy))  
#data.scores$site <- rownames(data.scores)  
#head(data.scores)


## 2003 saplings

#Filter for just year 2018
saplings03<-filter(saplings, YEAR == 2003)

#filter for species >5%
branch<-filter(saplings03, SPP =="PB"|SPP=="RM"|SPP=="RS"|SPP=="BF"|SPP=="EH"|SPP=="NONE")

#recalculate iv values
branch <- branch%>%
  mutate(ba.half = (0.5^2*0.005454)*X1.2.inch)%>%
  mutate(ba.one = (1.0^2*0.005454)*X1.inch)%>%
  mutate(ba.two = (2.0^2)*0.005454*X2.inch)%>%
  mutate(sap.ba = (ba.half+ba.one+ba.two)*250)

branch$SAP_EXP <- 250

branch<-replace(branch, is.na(branch), 0)   


branch<-branch%>%
  mutate(total_diameter = ba.half + ba.one + ba.two)

branch<-branch%>%
  mutate(BA = ((total_diameter^2)*0.005454))

branch<-branch%>%
  mutate(SAP_BA_TPA = BA * SAP_EXP)

plot_summary_sap <- branch%>%
  group_by(SITEid, PLOTid, YEAR)%>%
  summarise(SAP_TPA_total = sum(SAP_EXP),
            SAPBA_total = sum(SAP_BA_TPA))

branch <- branch%>%
  group_by(SITEid, PLOTid, YEAR, SPP)%>%
  summarise(SAP_TPA = sum(SAP_EXP),
            SAPBA = sum(SAP_BA_TPA))

branch<-left_join(branch, plot_summary_sap)


branch <- branch%>%
  mutate(prop_tpa = (SAP_TPA/SAP_TPA_total),
         prop_ba = (SAPBA/SAPBA_total))

branch <- branch%>%
  mutate(iv = ((prop_tpa + prop_ba)/2))

#create site-plot identifier
branch$sapID<-paste0(branch$SITEid,"-",branch$PLOTid)

#get data in long format
molten <- melt(as.data.frame(branch),id=c("sapID","iv","SPP"))
sapwide <- dcast(molten,sapID~SPP,value.var = "iv",mean)
sapwide[is.na(sapwide)] <- 0
head(sapwide)
# double check they all add to 1
rowSums(sapwide[2:7])


#NMDS
sapordi<-metaMDS(sapwide[,2:7], distance = "bray")
stressplot(sapordi)
plot(sapordi)

plot(sapordi,type="n")
points(sapordi,display="sites",cex=2,pch=21,col="red", bg="yellow")
text(sapordi,display="spec",cex=1.5,col="blue")



#calculate which species are below 5%
saplings03$sapID<-paste0(saplings03$SITEid,"-",saplings03$PLOTid)
sapwide <- dcast(saplings03,sapID~SPP, value.var = "X1.2.inch")

sapwide<-as.data.frame(sapwide)

sapwide[,2:12][sapwide[,2:12] > 0]<-1


sapwide<-sapwide%>%
  adorn_totals("row")
