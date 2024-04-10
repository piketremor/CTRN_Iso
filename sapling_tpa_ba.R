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


#setwd("G:/.shortcut-targets-by-id/1sCbm2t1PUIpbJYVOzlIrZKeVhisl4Ghv/CTRN_CFRU_Share/raw/csv")
#setwd("G:/My Drive/CTRN_CFRU_Share/raw/csv")
setwd("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv")
#setwd("~/Google Drive/My Drive/Research/CFRU/CTRN_CFRU_Share/raw/csv")
#saplings <- read.csv("~/Google Drive/My Drive/Research/CTRN_CFRU_Share/raw/csv/Saplings.csv") 
#saplings <- read.csv("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv/Saplings.csv")
saplings <- read.csv("Saplings.csv")
env<-read.csv("CTRN_Enviro.csv")

#select out for just the 15 CTRN sites
saplings <- filter(saplings, SITEid == "AS" | SITEid == "DR" | SITEid == "GR" | SITEid == "HR" | SITEid == "KI" | SITEid == "LM" | SITEid == "LT" | SITEid == "PA" | SITEid == "PE" | SITEid == "RC" | SITEid == "RR" | SITEid == "SA" | SITEid == "SC" | SITEid == "SR" | SITEid == "WB") 
saplings[saplings == "SpecAld"]<-"SA"
saplings[saplings == "HM"]<-"EH"
saplings[saplings == "CH"]<-"BC"

env<-filter(env, SITEid == "AS" | SITEid == "DR" | SITEid == "GR" | SITEid == "HR" | SITEid == "KI" | SITEid == "LM" | SITEid == "LT" | SITEid == "PA" | SITEid == "PE" | SITEid == "RC" | SITEid == "RR" | SITEid == "SA" | SITEid == "SC" | SITEid == "SR" | SITEid == "WB")





#stand metrics and importance values
saplings[is.na(saplings)] <- 0
saplings <- saplings%>%
  mutate(ba.half = (0.5^2*0.005454)*X1.2.inch)%>%
  mutate(ba.one = (1.0^2*0.005454)*X1.inch)%>%
  mutate(ba.two = (2.0^2)*0.005454*X2.inch)%>%
  mutate(sap.ba = (ba.half+ba.one+ba.two)*250)

saplings$SAP_EXP <- 250

saplings<-replace(saplings, is.na(saplings), 0)   
head(saplings)

saplings<-saplings%>%
  mutate(total_diameter = ba.half + ba.one + ba.two)

saplings<-saplings%>%
  mutate(BA = ((total_diameter^2)*0.005454))

saplings<-saplings%>%
  mutate(SAP_BA_TPA = BA * SAP_EXP)

plot_summary_sap <- saplings%>%
  group_by(SITEid, PLOTid, YEAR)%>%
  summarise(SAP_TPA_total = sum(SAP_EXP),
            SAPBA_total = sum(SAP_BA_TPA))

species_summary_sap <- saplings%>%
  group_by(SITEid, PLOTid, YEAR, SPP)%>%
  summarise(SAP_TPA = sum(SAP_EXP),
            SAPBA = sum(SAP_BA_TPA))

saplings_again<-left_join(species_summary_sap, plot_summary_sap)


saplings_again <- saplings_again%>%
  mutate(prop_tpa = (SAP_TPA/SAP_TPA_total),
         prop_ba = (SAPBA/SAPBA_total))

saplings_again <- saplings_again%>%
  mutate(iv = ((prop_tpa + prop_ba)/2))






#diversity 
saplings$spec.count <- (saplings$X1.2.inch+saplings$X1.inch+saplings$X2.inch)

sap.sum <- saplings%>%
  group_by(SITEid,PLOTid,YEAR,SPP)%>%
  summarize(sap.total = sum(spec.count))

all.in <- saplings%>%
  group_by(SITEid,PLOTid,YEAR)%>%
  summarize(sapling.total = sum(spec.count))

sappy <- left_join(sap.sum,all.in)
head(sappy)

trees <- read.csv("Trees2023.csv")
locs <- read.csv("Tree_locations_species.csv")
#trees <- read.csv("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv/Trees2023.csv")
#locs <- read.csv("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv/Tree_locations_species.csv")

tree_species <- locs[c(1:3,6)]
tree_species <- locs%>%
  select(SITEid, PLOTid, TREE, SPP)

#colnames(tree_species) <- c("SITEid", "PLOTid", "TREE", "SPP")

over <- left_join(trees, tree_species)
over$tally <- 1

over.sum<-over%>%
  group_by(SITEid, PLOTid, YEAR, SPP)%>%
  summarise(over.total = sum(tally))

all.over<-over%>%
  group_by(SITEid, PLOTid, YEAR)%>%
  summarise(overstory.total = sum(tally))

overstory<- left_join(over.sum, all.over)
head(overstory)

all_tree<-left_join(overstory, sappy, by = c("SITEid", "PLOTid", "SPP", "YEAR"))
all_tree$sap.total[is.na(all_tree$sap.total)] <- 0
all_tree$sapling.total[is.na(all_tree$sapling.total)] <- 0
all_tree$SPP[is.na(all_tree$SPP)]<-"OTHER CONIFER"

all_tree$prop <- (all_tree$over.total + all_tree$sap.total)/(all_tree$overstory.total+all_tree$sapling.total)
head(all_tree)

all_tree$shann.base <- all_tree$prop*(log(all_tree$prop))

shann.frame <- all_tree%>%
  group_by(SITEid, PLOTid, YEAR)%>%
  summarise(shannon = -1*sum(shann.base),
            hill = exp(shannon))

View(shann.frame)

sarah<-filter(shann.frame, SITEid == "SR" & YEAR == 2018)

#setwd("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv")

#write.csv(sarah,"SarahsRoadDiv.csv")

site.var<- read.csv("CTRN_Enviro.csv")
site.var <- read.csv("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv/CTRN_Enviro.csv")

shan <- left_join(shann.frame, site.var)
head(shan)
shan$ph[is.na(shan$ph)]<-0
shan <- filter(shan, ph > 0)
head(shan)
mod1<-lm(hill~YEAR+BAPA+TST+THIN_METH+delta.ba+ba.growth+elevation+wdi+RAD, data=shan)
summary(mod1)
#plot(mod1)

mdl.ac <- gls(hill~YEAR+BAPA+TST+THIN_METH+elevation+wdi+RAD+tmax+wd, data=shan, 
              correlation = corAR1(form=~1|YEAR/PLOTid/SITEid),
             na.action=na.omit)
summary(mdl.ac)
par(mfrow = c(2, 2))
plot(mdl.ac)



mod2<-lm(hill~YEAR+BAPA+TST+THIN_METH+elevation+wdi+RAD+tmax+wd, data=shan)
autoplot(mod2)
#linearity, scale location, Q-Q badish
mod3<-lm(log(hill)~YEAR+BAPA+TST+THIN_METH+elevation+wdi+RAD+tmax+wd, data=shan)
autoplot(mod3)
mod4<-glm(log(hill)~YEAR+BAPA+TST+THIN_METH+elevation+wdi+RAD+tmax+wd, data=shan,
         correlation = corARI(form=~1|YEAR/PLOTid/SITEid),
         na.action=na.omit)





sub1<-regsubsets(log(hill)~YEAR+BAPA+TST+THIN_METH+elevation+wdi+RAD+tmax+wd, data = shan)
summary(sub1)
sub2<-regsubsets(hill~YEAR+BAPA+TST+THIN_METH+elevation+wdi+RAD+tmax+wd, data=shan, 
          correlation = corAR1(form=~1|YEAR/PLOTid/SITEid),
          na.action=na.omit)





mydf2<-ggpredict(mdl.ac, terms = c("BAPA", "THIN_METH", "wd"))

#mydf2 <- ggpredict(mdl.ac,terms=c("TST","THIN_METH","wd"))


png("CTRN_Diversity_prelim.png",units='in',height=6,width=15,res=1000)
theme_set(theme_bw(16))


ggplot(mydf2,aes(x=x,y=predicted,colour=group))+
  geom_line(aes(linetype=group,color=group),size=1)+
  labs(linetype="Thinning Method")+
  labs(colour = "Thinning Method")+
  labs(x="Basal Area per Acre",y="Predicted Diversity (Hill)")+
  #ylim(0,55)+
  #xlim(4,6)+
  facet_wrap(~facet)+
  theme_bw(18) 
#theme(legend.position="none")
#scale_color_manual(values=c('gray0','gray70','gray40'))+
#scale_fill_manual(values=c('gray0','gray70','gray40'), name="fill")

dev.off()


#split the regeneration data into pre, post, and 10 yrs post treatment
plots <- read.csv("Plots.csv")
#control treatment year is unique to each installation take earliest treatment year

retreat <- left_join(saplings_again, plots)

shan <- shan%>%
  mutate(TST = YEAR - TRT_YR)

retreat<-retreat%>%
  mutate(TST = YEAR - TRT_YR)

repre<-filter(retreat, TST == -1)

repost<-filter(retreat, TST == 1)

reten <- filter(retreat, TST == 10) 
  #better?

##expected results
div.plot<-left_join(shann.frame,retreat)

ggplot(div.plot,aes(x=TST,y=hill))+
  geom_point()+
  labs(x="Time Since Thinning (years)",y="Predicted Diversity (Hill)")+
  theme_bw(18) 


ggplot(div.plot, aes(x=hill))+
  geom_histogram()+
  labs(x="Diversity (Hill numbers)")+
  ylab("")

ggplot(saplings, aes(x=total_diameter))+
  geom_histogram(bins = 10)+
  labs(x="Diameter (in)")+
  ylab("")

ggplot(shan, aes(x=hill))+
  stat_boxplot()+
  ylab("")+
  xlab("Hill Numbers")+
  theme_bw()+
  facet_wrap(vars(REMOVAL))+
  coord_flip()+
  easy_remove_x_axis()



 
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
env<-env[,c(2,3,4,27:65)]
env$sapID<-paste0(env$SITEid,"-",env$PLOTid)
env<-filter(env, YEAR==2018)
env<-env[,c(4:43)]
env<-env[,c(40, 1:39)]
env[is.na(env)] <- 0

sap2<-left_join(sapwide, env)
sap2[is.na(sap2)] <- 0

sap.rda


env<-env[,2:40]
sapwide<-sapwide[,2:12]
summary(env)

sapwide<-as.data.frame(sapwide)
env<-as.data.frame(env)

env<-as.numeric(env[,2:40])

env<-rownames(env[,1])

env$flowdir<-as.numeric(env$flowdir)
env$Parent<-as.numeric(env$Parent)

env$sapID<-as.factor(env$sapID)
sapwide$sapID<-as.factor(sapwide$sapID)

name<-sapwide$sapID

list<-env$sapID

bark<-left_join(env, sapwide)

sapwide<-bark[c(1, 41:51)]
env<-bark[,1:40]

sapwide<-sapwide[,2:12]
env<-env[,2:40]
#rda
sap.rda <- rda(sapwide ~ ., data=env)
summary(sap.rda)
ordiplot(sap.rda, scaling = 2, type = "text")
step.forward <- ordiR2step(rda(sapwide ~ 1, data=env), scope=formula(sap.rda), R2scope = F, direction="forward", pstep=1000)
anova(sap.rda, by="terms", step=1000)
#reduced model
sap.rda2<-rda(sapwide ~ wd+SWI+ppt+SWC2+WD2020+Winds10+vpdmax+tri+dep, data=env)
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


sapwide<-filter(sapwide, BF > 0 | EH > 0 | NONE > 0 | PB > 0 | RM > 0 | RS > 0)

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
