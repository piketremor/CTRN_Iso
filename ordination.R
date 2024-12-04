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


trees <- read.csv("Trees2023.csv")
locs <- read.csv("Tree_locations_species.csv")
tree_species <- locs[c(1:3,6)]
over <- left_join(trees, tree_species)
overstory <- over[c(2:9,27)]
overstory$DBH[is.na(overstory$DBH)] <- 0
overstory <- dplyr::filter(overstory,DBH>2.5)
overstory$SPP[is.na(overstory$SPP)] <- "OH"
overstory[overstory == "WB"]<-"PB"
overstory[overstory == "MA"]<-"OH"
ht.mod <- nlme(TOT_HT~4.5+exp(a+b/(DBH+1)), # nonlinear mixed model, fitting a different HT~DBH based on species
               data=overstory,fixed=a+b~1,random=a+b~1|SPP,na.action=na.pass,
               start=c(a=4.5,b=-6),control=nlmeControl(returnObject = TRUE,msMaxIter = 10000,maxIter = 5000))
ov <- overstory%>%
  mutate(ht.fit = predict(ht.mod,.),ba = DBH^2*0.005454,ef = 5,crown.width = mapply(MCW,SPP=SPP,DBH=DBH))
ht.40 <- ov%>%
  filter(.,SPP=="RS")%>%
  group_by(SITEid,PLOTid,YEAR,SPP)%>%
  top_n(8,DBH)%>%
  summarize(ht40 = mean(ht.fit))
overstory.summary <- ov%>%
  group_by(SITEid,PLOTid,YEAR)%>%
  summarize(bapa = sum(ba*ef),tpa = sum(ef),qmd = qmd(bapa,tpa),RD = bapa/sqrt(qmd),
            CCF = (sum((((crown.width/2)^2)*3.14)*ef))/43560)
sp.summary <- ov%>%
  group_by(SITEid,PLOTid,YEAR,SPP)%>%
  summarize(sp.bapa = sum(ba*ef),sp.tpa = sum(ef))
shannon <- overstory.summary%>%
  left_join(.,sp.summary)%>%
  mutate(iv = ((sp.bapa/bapa+sp.tpa/tpa)/2),Shannons = iv*log(iv))%>%
  group_by(SITEid,PLOTid,YEAR)%>%
  summarize(Shannon = sum(Shannons)*-1,over.Hill = exp(Shannon))
over.standlist <- left_join(overstory.summary,shannon)%>%
  left_join(.,ht.40)
# need to get the % of spruce, fir, and HW
spec.groups <- overstory.summary%>%
  left_join(.,sp.summary)%>%
  mutate(prop.ws = ifelse(SPP=="WS",(sp.bapa/bapa),0),
         prop.bs = ifelse(SPP=="BS",(sp.bapa/bapa),0),
         prop.rs = ifelse(SPP=="RS",(sp.bapa/bapa),0),
         prop.eh = ifelse(SPP=="EH",(sp.bapa/bapa),0),
         prop.ab = ifelse(SPP=="AB",(sp.bapa/bapa),0),
         prop.bf = ifelse(SPP=="BF",(sp.bapa/bapa),0),
         prop.hw = ifelse(SPP=="BC"|SPP=="PB"|SPP=="QA"|SPP=="RM"|SPP=="GB"|SPP=="YB"|SPP=="ST"|
                            SPP=="WILLOW"|SPP=="PC"|SPP=="BA"|SPP=="MA"|SPP=="OA"|SPP=="OH"|SPP=="WA"|SPP=="SM",sp.bapa/bapa,0))%>%
  group_by(SITEid,PLOTid,YEAR)%>%
  summarize(prop.ws.avg = sum(prop.ws),prop.bs.avg = sum(prop.bs),prop.rs.avg = sum(prop.rs),
            prop.hw.avg = sum(prop.hw),prop.ab.avg = sum(prop.ab),prop.eh.avg = sum(prop.eh),prop.bf.avg = sum(prop.bf))
over.df <- left_join(over.standlist,spec.groups)
##
site.vars <- read.csv("CTRN_SiteVariables_20240718.csv")
treat <- read.csv("CTRN_Enviro2_20240815.csv")
trt <- treat[c(1,2,4,7,8)]
sites <- left_join(over.df,site.vars)%>%
  left_join(.,trt)
wd <- read.csv("CTRN_WD_dataframe.csv")
wd <- dplyr::rename(wd,SITEid=SiteID,PLOTid=PlotID,YEAR=year)
mean.wd <- wd%>%
  group_by(SITEid,PLOTid)%>%
  summarize(mean.WD=mean(WD))
run.wd <- wd%>%
  left_join(.,trt)%>%
  mutate(tst = YEAR-TRT_YR)%>%
  filter(.,tst>0)%>%
  mutate(wd.time = cumsum(WD))
running <- run.wd[c(1,2,4,11,12)]
water.frame <- left_join(mean.wd,running)
cleaned.over <- left_join(sites,water.frame)
cleaned.over$Northing_Y[is.na(cleaned.over$Northing_Y)] <- 0
cleaned.over <- filter(cleaned.over,Northing_Y>0)
actual.rem <- read.csv("trt_list.csv")
act <- actual.rem[c(1,2,6)]
final.over <- left_join(cleaned.over,act)
final.over$wdi.time <- final.over$wd.time-final.over$SWC2
final.over<-select(final.over, -c(SPP,Lithic,Northing_Y, Easting_X,Redox,Parent,Min_depth,WD2000,WD2020))
final.over$THIN_METH<-as.factor(final.over$THIN_METH)
final.over$actual.removed[is.na(final.over$actual.removed)]<-0
names(final.over)
final.over<-left_join(saplings,final.over)
names(final.over)
final.over$sapID<-paste0(final.over$SITEid,"-",final.over$PLOTid,"-",final.over$CORNERid)
#change from integer to numeric
final.over$flowdir<-as.numeric(final.over$flowdir)
names(final.over)
final.over<-na.omit(final.over)

##2018 Ordination and Data Cleaning

#Filter for just year 2018
saplings18<-filter(saplings, YEAR == 2018)
unique(saplings18$SPP)

#filter for species >5%
branch<-filter(saplings18, SPP =="PB"|SPP == "OT"|SPP=="NC"|SPP=="RM"|SPP=="RS"|SPP=="WP"|SPP=="BF"|SPP=="YB"|SPP=="EH"|SPP=="QA"|SPP=="BC")

unique(branch$SPP)
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
  group_by(SITEid, PLOTid, CORNERid,YEAR)%>%
  summarise(SAP_TPA_total = sum(SAP_EXP),
            SAPBA_total = sum(SAP_BA_TPA))

branch <- branch%>%
  group_by(SITEid, PLOTid, CORNERid,YEAR, SPP)%>%
  summarise(SAP_TPA = sum(SAP_EXP),
            SAPBA = sum(SAP_BA_TPA))

branch<-left_join(branch, plot_summary_sap)


branch <- branch%>%
  mutate(prop_tpa = (SAP_TPA/SAP_TPA_total),
         prop_ba = (SAPBA/SAPBA_total))

branch <- branch%>%
  mutate(iv = ((prop_tpa + prop_ba)/2))


#create site-plot identifier
branch$sapID<-paste0(branch$SITEid,"-",branch$PLOTid,"-",branch$CORNERid)

#get data in long format
molten <- melt(as.data.frame(branch),id=c("sapID","iv","SPP"))
sapwide <- dcast(molten,sapID~SPP,value.var = "iv",mean)
sapwide[is.na(sapwide)] <- 0
head(sapwide)
# double check they all add to 1
rowSums(sapwide[2:12])
names<-sapwide$sapID
rownames(sapwide)<-names

#NMDS
sapordi<-metaMDS(sapwide[,2:12], distance = "bray")
stressplot(sapordi)
sapordi
plot(sapordi)

plot(sapordi,type="n")
points(sapordi,display="sites",cex=2,pch=21,col="red", bg="yellow")
text(sapordi,display="spec",cex=1.5,col="blue")

speciesscores<-as.data.frame(scores(sapordi)$species)
speciesscores$species<-rownames(speciesscores)
head(speciesscores)
sitescores<-as.data.frame(scores(sapordi)$sites)
sitescores$site<-rownames(sitescores)
head(sitescores)

#could add in a grouping variable (like shade tolerant/intolerant) to create a polygon grouping 
ggplot() + 
  geom_text(data=speciesscores,aes(x=NMDS1,y=NMDS2,label=species),size=5,vjust=0, colour="blue") +  
  #geom_point(data=sitescores,aes(x=NMDS1,y=NMDS2),size=1) + # add the point markers
  geom_text(data=sitescores,aes(x=NMDS1,y=NMDS2,label=site),size=3,vjust=0, color="red") +  
  coord_equal() +
  theme_bw()+
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),  
        axis.title.x = element_text(size=18), 
        axis.title.y = element_text(size=18), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.background = element_blank())



## Constrained Ordination###
#join together sapling and environmental data to make sure the number of observations match
final.over<-na.omit(final.over)
bark<-left_join(sapwide, final.over)
#bark[is.na(bark)] <- 0
names<-bark$sapID #create list of sapID to set the rownames with 
#separate the datasets back out
sapwide<-bark[,c(1:12)]
env<-bark[,21:74]
names(sapwide)
names(env)
#add rownames
#rownames(sapwide)<-names
#rownames(env)<-names
#remove sapID variable
sapwide<-sapwide[,2:12]
names(sapwide)
names(env)
head(env)
head(sapwide)
sapwide[is.na(sapwide)] <- 0
head(sapwide)


#rda with all environmental variables (0.39)
sap.rda <- rda(sapwide ~ ., data=env, na.action = na.exclude) 
summary(sap.rda)
ordiplot(sap.rda, scaling = 2, type = "text")
RsquareAdj(sap.rda)$adj.r.squared

#env<-na.omit(env)
#there is a group of NA values. Need to figure out how to join the datasets together to maintain the corner but not the NAs
#variable selection
step.forward <- ordiR2step(rda(sapwide ~ 1, data=env, na.action = na.exclude), scope=formula(sap.rda,na.action=na.exclude), R2scope = F, direction="forward", pstep=1000,na.action=na.exclude)
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


#model with VSURF selected variables from diversity code
sap.rda.vsurf<-rda(sapwide~TPA_TOTAL+Winds50+ppt,data=env)
summary(sap.rda.vsurf)
ordiplot(sap.rda.vsurf, scaling = 2, type = "text")
RsquareAdj(sap.rda.vsurf)$adj.r.squared









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


