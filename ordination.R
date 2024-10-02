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
env<-read.csv("CTRN_Env_vars2.csv")

saplings <- filter(saplings, SITEid == "AS" | SITEid == "DR" | SITEid == "GR" | SITEid == "HR" | SITEid == "KI" | SITEid == "LM" | SITEid == "LT" | SITEid == "PA" | SITEid == "PE" | SITEid == "RC" | SITEid == "RR" | SITEid == "SA" | SITEid == "SC" | SITEid == "SR" | SITEid == "WB") 
saplings[saplings == "SpecAld"]<-"SA"
saplings[saplings == "HM"]<-"EH"
saplings[saplings == "CH"]<-"BC"

env<-filter(env, SITEid == "AS" | SITEid == "DR" | SITEid == "GR" | SITEid == "HR" | SITEid == "KI" | SITEid == "LM" | SITEid == "LT" | SITEid == "PA" | SITEid == "PE" | SITEid == "RC" | SITEid == "RR" | SITEid == "SA" | SITEid == "SC" | SITEid == "SR" | SITEid == "WB")

#calculate overstory hil numbers
trees <- read.csv("Trees2023.csv")
locs <- read.csv("Tree_locations_species.csv")


tree_species <- locs[c(1:3,6)]
tree_species <- locs%>%
  select(SITEid, PLOTid, TREE, SPP)



over <- left_join(trees, tree_species)
over <- filter(over, SITEid == "AS" | SITEid == "DR" | SITEid == "GR" | SITEid == "HR" | SITEid == "KI" | SITEid == "LM" | SITEid == "LT" | SITEid == "PA" | SITEid == "PE" | SITEid == "RC" | SITEid == "RR" | SITEid == "SA" | SITEid == "SC" | SITEid == "SR" | SITEid == "WB") 

over.sum<-over%>%
  mutate(ef = 5)%>%
  group_by(SITEid, PLOTid, YEAR, SPP)%>%
  reframe(over.spp.total = sum(ef))

all.over<-over%>%
  mutate(ef=5)%>%
  group_by(SITEid, PLOTid, YEAR)%>%
  reframe(overstory.total = sum(ef))

head(all.over)

overstory<- left_join(over.sum, all.over)
overstory

overstory$over.spp.total[is.na(overstory$over.spp.total)]<-0
overstory$overstory.total[is.na(overstory$overstory.total)]<-0
overstory$ov.prop<-(overstory$over.spp.total/overstory$overstory.total)


overstory$ov.shann.base<-overstory$ov.prop*(log(overstory$ov.prop))
overstory

ov.dv<-overstory%>%
  group_by(SITEid,PLOTid,YEAR)%>%
  reframe(ov.Shannon = sum(ov.shann.base)*-1,
          ov.Hill = exp(ov.Shannon))
ov.dv

#overstory1<-left_join(overstory,ov.dv)

#add overstory hill numbers to env 

env<-left_join(env,ov.dv)
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
  mutate(sap.ba = (ba.half+ba.one+ba.two)*62.5)

branch$SAP_EXP <- 62.5

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
#env<-env[,c(2,3,4,27:65)]
env$sapID<-paste0(env$SITEid,"-",env$PLOTid)
env<-filter(env, YEAR==2018)
#env<-env[,c(4:43)]
#env<-env[,c(40, 1:39)]
#env[is.na(env)] <- 0

#change from integer to numeric
env$flowdir<-as.numeric(env$flowdir)
env$Parent<-as.numeric(env$Parent)

#join together sapling and environmental data to make sure the number of observations match
bark<-left_join(env, sapwide)
#bark[is.na(bark)] <- 0
names<-bark$sapID #create list of sapID to set the rownames with 

#separate the datasets back out
sapwide<-bark[,c(51:62)]
env<-bark[,4:50]

#add rownames
rownames(sapwide)<-names
rownames(env)<-names

#remove sapID variable
sapwide<-sapwide[,2:12]



#env <- subset(env, select = -sapID)

head(env)
head(sapwide)

env<-env%>% #converting all WD variables to water surplus variables and calculating WDI
  mutate(
    WDI = (ave.WD-WHC)*-1,
    WD = WD*-1, 
    cumulative.WD=cumulative.WD*-1, 
    ave.WD=ave.WD*-1,
    run.wd=run.wd*-1
    
  )


sapwide[is.na(sapwide)] <- 0
head(sapwide)


#rda with all environmental variables (0.39)
sap.rda <- rda(sapwide ~ ., data=env, na.action = na.exclude) 
summary(sap.rda)
ordiplot(sap.rda, scaling = 2, type = "text")
RsquareAdj(sap.rda)$adj.r.squared

#variable selection
step.forward <- ordiR2step(rda(sapwide ~ 1, data=env, na.action = na.exclude), scope=formula(sap.rda), R2scope = F, direction="forward", pstep=1000)
anova(sap.rda, by="terms", step=1000) 
anova(sap.rda, step=1000)

##Test individual environmental constraints

###water deficit
rda.wd<-rda(sapwide~WD,data=env,na.action=na.exclude) #0.06
summary(rda.wd)
RsquareAdj(rda.wd)$adj.r.squared

#water deficit index
rda.wdi<-rda(sapwide~WDI,data=env,na.action=na.exclude) #0.05
summary(rda.wdi)
RsquareAdj(rda.wdi)$adj.r.squared

###average water deficit ##winner
rda.avg.wd<-rda(sapwide~ave.WD,data=env,na.action=na.exclude) #0.098
summary(rda.avg.wd)
RsquareAdj(rda.avg.wd)$adj.r.squared

#cumulative water deficit
rda.cumwd<-rda(sapwide~cumulative.WD,data=env,na.action=na.exclude) #0.08
summary(rda.cumwd)
RsquareAdj(rda.cumwd)$adj.r.squared

#running water deficit
rda.runwd<-rda(sapwide~run.wd,data=env,na.action=na.exclude) #0.03
summary(rda.runwd)
RsquareAdj(rda.runwd)$adj.r.squared

#water holding capacity 
rda.whc<-rda(sapwide~WHC,data=env,na.action=na.exclude) #0.049
summary(rda.whc)
RsquareAdj(rda.whc)$adj.r.squared

#reduced model

#model with the variables that were included for the diversity model and average water deficit (0.134)
sap.rda.div <- rda(sapwide ~ dew+THIN_METH+ave.WD+actual.removed+ov.Hill, na.action = na.exclude, data=env)
summary(sap.rda.div)
ordiplot(sap.rda.div, scaling = 2, type = "text")
anova(sap.rda.div, by="terms", step=1000)
anova(sap.rda.div, step=1000)
anova(sap.rda.div, by="axis", step=1000)
RsquareAdj(sap.rda.div)$adj.r.squared


rda.test<-rda(sapwide ~ vpdmax+THIN_METH+ave.WD+actual.removed+WDI, na.action = na.exclude, data=env)
summary(rda.test)
ordiplot(rda.test, scaling = 2, type = "text")
RsquareAdj(rda.test)$adj.r.squared

#model with all the variables from the forward selection procedure (0.29)
sap.rda2<-rda(sapwide ~ ave.WD+WHC+ppt+SWC2+vpdmax+WDI+run.wd+cumulative.WD, data=env)
summary(sap.rda2)
ordiplot(sap.rda2, scaling = 2, type = "text")
anova(sap.rda2, by="terms", step=1000)
anova(sap.rda2, step=1000)
RsquareAdj(sap.rda2)$adj.r.squared


#model without ave.WD (0.25)
sap.rda3 <- rda(sapwide ~ WHC+ppt+vpdmax+WDI+run.wd+cumulative.WD, data=env)
summary(sap.rda3)
ordiplot(sap.rda3, scaling = 2, type = "text")
anova(sap.rda3, by="terms", step=1000)
anova(sap.rda3, step=1000)
RsquareAdj(sap.rda3)$adj.r.squared

#model without cumulative.WD (0.24)
sap.rda4<-rda(sapwide ~ ave.WD+WHC+ppt+vpdmax+WDI+run.wd, data=env)
summary(sap.rda4)
ordiplot(sap.rda4, scaling = 2, type = "text")
anova(sap.rda4, by="terms", step=1000)
anova(sap.rda4, step=1000)
anova(sap.rda4, by = "axis", step = 1000)
RsquareAdj(sap.rda4)$adj.r.squared

#model withouth SWC2 and with run.wd instead of cumulative wdi (0.24)
sap.rda5<-rda(sapwide ~ ave.WD+WHC+ppt+vpdmax+WDI+run.wd, data=env)
summary(sap.rda5)
ordiplot(sap.rda5, scaling = 2, type = "text")
anova(sap.rda5, by="terms", step=1000)
anova(sap.rda5, step=1000)
anova(sap.rda5, by = "axis", step = 1000)
RsquareAdj(sap.rda5)$adj.r.squared











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


