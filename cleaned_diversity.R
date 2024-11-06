############################################################################
############################################################################
#### Generate predictive models of species diversity from CTRN database ####
############################################################################
############################################################################

# 2024 08 28 - Code originated from Lila Beck and Mike Premer


#load packages
dev.off()
rm(list=ls())
library(dplyr)
library(mosaic)
library(forcats)
#library(tidyverse)
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
library(VSURF)

setwd("~/Google Drive//My Drive/CTRN_CFRU_Share/raw/csv")
setwd("~/Google Drive/My Drive/Research/CFRU/CTRN_CFRU_Share/raw/csv")

saplings <- read.csv("Saplings.csv")

#select out for just the 15 CTRN sites
saplings <- filter(saplings, SITEid == "AS" | SITEid == "DR" | SITEid == "GR" | SITEid == "HR" | SITEid == "KI" | SITEid == "LM" | SITEid == "LT" | SITEid == "PA" | SITEid == "PE" | SITEid == "RC" | SITEid == "RR" | SITEid == "SA" | SITEid == "SC" | SITEid == "SR" | SITEid == "WB") 
saplings[saplings == "SpecAld"]<-"SA"
saplings[saplings == "HM"]<-"EH"
saplings[saplings == "CH"]<-"BC"

#diversity 
saplings$spec.count <- (saplings$X1.2.inch+saplings$X1.inch+saplings$X2.inch)

sap.sum <- saplings%>%
  mutate(ef=62.5)%>%
  group_by(SITEid,PLOTid,YEAR,SPP)%>%
  reframe(sap.spp.total = sum(spec.count*ef))

all.in <- saplings%>%
  mutate(ef=62.5)%>%
  group_by(SITEid,PLOTid,YEAR)%>%
  summarize(sapling.total = sum(spec.count*ef))

sappy <- left_join(all.in,sap.sum)
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

#####stop here, this is where saplings and overstory get joined together####

#sapling diversity metrics

sappy$sapling.total[is.na(sappy$sapling.total)]<-0
sappy$sap.spp.total[is.na(sappy$sap.spp.total)]<-0
sappy$sap.prop<-(sappy$sap.spp.total/sappy$sapling.total)
sappy

sappy$sap.shann.base<-sappy$sap.prop*(log(sappy$sap.prop))
sappy

sap.dv <- sappy%>%
  group_by(SITEid,PLOTid,YEAR)%>%
  reframe(sap.Shannon = sum(sap.shann.base)*-1,
          sap.Hill = exp(sap.Shannon))
sap.dv


#overstory
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

#all_tree<-left_join(overstory, sappy, by = c("SITEid", "PLOTid", "SPP", "YEAR"))
#all_tree
#all_tree$sapling.total[is.na(all_tree$sapling.total)] <- 0
#all_tree$sap.spp.total[is.na(all_tree$sap.spp.total)] <- 0
#all_tree$SPP[is.na(all_tree$SPP)]<-"OTHER CONIFER"
#all_tree
#all_tree$total.trees <- (all_tree$overstory.total+all_tree$sapling.total)
#all_tree
#all_tree$total.species <- (all_tree$over.spp.total+all_tree$sap.spp.total)
#all_tree
#plot(all_tree$total.trees,all_tree$total.species)
#all_tree$prop <- (all_tree$total.species/all_tree$total.trees)


# brute force patch

#all_tree$total.species <- ifelse(all_tree$total.species>all_tree$total.trees,all_tree$total.trees,all_tree$total.species)
#all_tree$shann.base <- all_tree$prop*(log(all_tree$prop))

#dv <- all_tree%>%
#  group_by(SITEid,PLOTid,YEAR)%>%
#  reframe(Shannon = sum(shann.base)*-1,
  #        Hill = exp(Shannon))


#head(dv)
# cleaned. 
site.vars <- read.csv("CTRN_SiteVariables_20240718.csv")
#site.vars <- dplyr::rename(site.vars,SITEid=SiteID,PLOTid=PlotID)
#head(site.vars)
#site.dv <- left_join(dv,site.vars)


#join to each dataframe
site.sap.dv<-left_join(sap.dv,site.vars)
site.ov.dv<-left_join(ov.dv,site.vars)

wd <- read.csv("CTRN_WD_dataframe.csv")
wd <- dplyr::rename(wd,SITEid=SiteID,PLOTid=PlotID,YEAR=year)
#site.water <- left_join(site.dv,wd)
#head(site.water)

site.water.sap<-left_join(site.sap.dv,wd)
site.water.ov<-left_join(site.ov.dv,wd)
#clean.df <- distinct(site.water)

clean.sap<-distinct(site.water.sap)
clean.ov<-distinct(site.water.ov)

# finally, we need the treatment data. 

### 8/15/2024 ###
##################
##################
## HOLD THE PHONE!!! WE DON'T HAVE CONTROL PLOTS IN THE FOLLOWING FILE - NEED TO PATCH THIS - IT IS A PROBLEM. 

# 8/15/2024 - patched


treat <- read.csv("CTRN_Enviro2_20240815.csv")
head(treat)
#treat <- treat%>%
#  group_by(SITEid,PLOTid)%>%
#  top_n(.,1,wt=YEAR)
#treat.rx <- treat[c(2,3,10:14)]
#cleaned <- left_join(clean.df,treat,by=c("SITEid","PLOTid"))
cleaned.sap<-left_join(clean.sap,treat,by=c("SITEid", "PLOTid"))
cleaned.ov<-left_join(clean.ov,treat,by=c("SITEid", "PLOTid"))

cleaned.sap$tst<-cleaned.sap$YEAR-cleaned.sap$TRT_YR
cleaned.ov$tst<-cleaned.ov$YEAR-cleaned.ov$TRT_YR

#cleaned$tst <- cleaned$YEAR-cleaned$TRT_YR
cleaned.sap$tst<-cleaned.sap$YEAR-cleaned.sap$TRT_YR
cleaned.ov$tst<-cleaned.ov$YEAR-cleaned.ov$TRT_YR


#cleaned$X[is.na(cleaned$X)] <- 0
#cleaned <- dplyr::filter(cleaned,X>0)


# changing these back to keep NA, otherwise we end up with 0 inflated predictors
#cleaned$Redox[is.na(cleaned$Redox)] <-0
#cleaned$Lithic[is.na(cleaned$Lithic)] <-0
#cleaned$Densic[is.na(cleaned$Densic)] <-0
#cleaned$Profile[is.na(cleaned$Profile)] <-0
#cleaned$Planform[is.na(cleaned$Planform)] <-0
#cleaned$FERT <- 0

cleaned.ov$FERT<-0
cleaned.sap$FERT<-0
#cleaned$REMOVAL <- as.factor(cleaned$REMOVAL)
#cleaned$PCT <- as.factor(cleaned$PCT)

cleaned.ov$PCT<-as.factor(cleaned.ov$PCT)
cleaned.sap$PCT<-as.factor(cleaned.sap$PCT)

#cleaned$THIN_METH <- as.factor(cleaned$THIN_METH)

cleaned.ov$THIN_METH<-as.factor(cleaned.ov$THIN_METH)
cleaned.sap$THIN_METH<-as.factor(cleaned.sap$THIN_METH)
#head(cleaned)
head(cleaned.ov)
head(cleaned.sap)

#wd.start.frame <- cleaned%>%
#  group_by(SITEid,PLOTid,TRT_YR)%>%
#  summarize(c = mean(elevation))
#wd.start.frame <- wd.start.frame[1:3]
#wd.start.frame <- dplyr::rename(wd.start.frame,START_YR=TRT_YR)
#wd.start.frame
#cleaned <- left_join(cleaned,wd.start.frame)
#wd.calibrate <- dplyr::filter(cleaned,tst>-1)

#wd.calib <- wd.calibrate%>%
#  group_by(SITEid,PLOTid,TRT_YR)%>%
#  mutate(run.wd = cumsum(WD))

#head(wd.calib)
#wd.df <- wd.calib[c(1,2,3,65)]

#overstory
wd.start.frame1 <- cleaned.ov%>%
  group_by(SITEid,PLOTid,TRT_YR)%>%
  summarize(c = mean(elevation))
wd.start.frame1 <- wd.start.frame1[1:3]
wd.start.frame1 <- dplyr::rename(wd.start.frame1,START_YR=TRT_YR)
wd.start.frame1
cleaned.ov <- left_join(cleaned.ov,wd.start.frame1)
wd.calibrate1 <- dplyr::filter(cleaned.ov,tst>-1)

wd.calib1 <- wd.calibrate1%>%
  group_by(SITEid,PLOTid,TRT_YR)%>%
  mutate(run.wd = cumsum(WD))

cleaner.ov<-left_join(cleaned.ov,wd.calib1)
#saplings
wd.start.frame2 <- cleaned.sap%>%
  group_by(SITEid,PLOTid,TRT_YR)%>%
  summarize(c = mean(elevation))
wd.start.frame2 <- wd.start.frame2[1:3]
wd.start.frame2 <- dplyr::rename(wd.start.frame2,START_YR=TRT_YR)
wd.start.frame2
cleaned.sap <- left_join(cleaned.sap,wd.start.frame2)
wd.calibrate2 <- dplyr::filter(cleaned.sap,tst>-1)

wd.calib2 <- wd.calibrate2%>%
  group_by(SITEid,PLOTid,TRT_YR)%>%
  mutate(run.wd = cumsum(WD))

cleaner.sap<-left_join(cleaned.sap,wd.calib2)
#wd.calib[c(1,2,47,55,65,66)]



#cleaner <- left_join(cleaned,wd.df)
#cleaner <- dplyr::filter(cleaner,Shannon>0)

cleaner.sap<-dplyr::filter(cleaned.sap,sap.Shannon>0)

cleaner.ov<-dplyr::filter(cleaned.ov,ov.Shannon>0)

actual.rem <- read.csv("trt_list.csv")
act <- actual.rem[c(1,2,6)]

#cleaner <- left_join(cleaner,act)
cleaner.sap<-left_join(cleaner.sap,act)
cleaner.ov<-left_join(cleaner.ov,act)
#names(cleaner)

###add in overstory hill numbers to sapling dataset
over.hill<-select(cleaner.ov, SITEid, PLOTid, ov.Hill)
cleaner.sap<left_join(cleaner.sap, over.hill,by=c("SITEid","PLOTid"))

#model.null<- lm(Hill~elevation+tri+tpi+roughness+slope+aspect+flowdir+tmin+tmean+tmax+dew+vpdmax+
 #                  vpdmin+McNab+Bolstad+Profile+Planform+Winds10+Winds50+SWI+RAD+ppt+Parent+Densic+Lithic+
  #                 Redox+Min_depth+WD+cumulative.WD+ave.WD+WHC+ex.mg+ex.ca+ph+dep+ex.k+nit+SWC2+PCT+REMOVAL+THIN_METH+tst+run.wd+actual.removed, data= cleaner, na.action=na.omit)
#summary(model.null)
require(MASS)
#mod.step <- model.null%>%
#  stepAIC(trace=FALSE,na.action="na.omit")
#summary(mod.step)
require(performance)
#check_collinearity(mod.step)

#mod1 <- lm(Hill~REMOVAL+THIN_METH+run.wd+actual.removed,data=cleaner)
#summary(mod1)
#hist(cleaner$Hill)
#no.na.data <- dplyr::filter(no.na.data,Hill>0)
#summary(mod1)
require(leaps)
#models.ov <- regsubsets(Hill~elevation+tri+tpi+roughness+slope+aspect+flowdir+tmin+tmean+tmax+dew+vpdmax+
#                          vpdmin+McNab+Bolstad+Profile+Planform+Winds10+Winds50+SWI+RAD+ppt+Parent+Densic+Lithic+
#                          Redox+Min_depth+WD+cumulative.WD+ave.WD+WHC+ex.mg+ex.ca+ph+dep+ex.k+nit+SWC2+PCT+REMOVAL+THIN_METH+tst+run.wd+actual.removed,really.big = TRUE,
#                        data = cleaner,method="exhaustive")
names(cleaner.ov)
ov.hill <- cleaner.ov[c(1,2,3,5)]

# great, now you can either run the overstory first, or just append the overstory Hill number to the understory dataframe
## 
cleaner.sap <- left_join(cleaner.sap,ov.hill)
cleaner.sap$WDI <- cleaner.sap$ave.WD-cleaner.sap$SWC2
cleaner.sap$cumulative.WDI <- cleaner.sap$cumulative.WD-cleaner.sap$SWC2


## add overstory structural variables here. 
overstory_metrics<-read.csv("overstory_metrics.csv")
cleaner.sap<-left_join(cleaner.sap, overstory_metrics)



models.sapling<-regsubsets(sap.Hill~elevation+tri+tpi+roughness+slope+aspect+flowdir+tmin+tmean+tmax+dew+vpdmax+
                             vpdmin+McNab+Bolstad+Profile+Planform+Winds10+Winds50+SWI+RAD+ppt+Parent+WD+cumulative.WD+ave.WD+WHC+ex.mg+ex.ca+ph+dep+ex.k+nit+SWC2+PCT+REMOVAL+THIN_METH+tst+actual.removed+ov.Hill+WDI+cumulative.WDI+BAPA+TPA+pp.spruce+pp.hw+pp.bf,really.big = TRUE,
                           data = cleaner.sap,method="exhaustive")

summary(models.sapling)
models.sapling$ress
plot.sap <- lm(sap.Hill~ov.Hill,data=cleaner.sap)
summary(plot.sap)



res.sum2 <- summary(models.sapling)
which.max(res.sum2$adjr2)
data.frame(
  Adj.R2 = which.max(res.sum2$adjr2),
  CP = which.min(res.sum2$cp),
  BIC = which.min(res.sum2$bic)
)

#VSURF
names(cleaner.sap)
cleaner.sap<-select(cleaner.sap, -c("OtherPlotName", "ThinnedYet", "Notes", "Active", "SitePlotID", "PlotArea", "PlotDimensions", "X_Coordinate", "Y_Coordinate", "START_YR", "Northing_Y", "Easting_X", "Densic", "Redox", "Lithic", "Redox"))
names(cleaner.sap)
cleaner.sap<-as.data.frame(cleaner.sap)

#just a test
cleaner.sap<-cleaner.sap %>% replace(is.na(.), 0)
preds <- cleaner.sap[,c(6:63)]
obs <- cleaner.sap[,5]

sapling.vsurf<-VSURF(preds, obs) 

sapling.vsurf
sapling.vsurf$varselect.thres
sapling.vsurf$varselect.interp
sapling.vsurf$nums.varselect
sapling.vsurf$varselect.pred
vars <- preds[c(22,19,49,57,10,1,44)] 
names(vars)


v.lm.sap<-lm(sap.Hill~ppt+Winds50+TPA_TOTAL+relative.density+tmax+elevation+tst,data=cleaner.sap)
summary(v.lm.sap)
check_collinearity(v.lm.sap)

temp.lm<-lm(sap.Hill~ppt+Winds50+TPA_TOTAL+relative.density+tmax+tst,data=cleaner.sap)
summary(temp.lm) #better

elv.lm<-lm(sap.Hill~ppt+Winds50+TPA_TOTAL+relative.density+elevation+tst,data=cleaner.sap)
summary(elv.lm)

v.sap.model <- lme(sap.Hill~ppt+TPA_TOTAL+Winds50+relative.density+tmax+tst,
             data=cleaner.sap,
             correlation=corAR1(form=~YEAR|SITEid/PLOTid),
             random=~1|SITEid/PLOTid,
             na.action=na.omit,method="REML")
summary(v.sap.model)
plot(v.sap.model)


# what if we do this in chunks.... 
## treat, site, soil, meteor, climate



## for the understory, recommend adding in overstory plot summaries like.... BAPA, TPA, QMD, RD, CCF, %Spruce, %Fir



models.overstory<-regsubsets(ov.Hill~elevation+tri+tpi+roughness+slope+aspect+flowdir+tmin+tmean+tmax+dew+vpdmax+
                               vpdmin+McNab+Bolstad+Profile+Planform+Winds10+Winds50+SWI+RAD+ppt+Parent+Densic+Lithic+
                               Redox+Min_depth+WD+cumulative.WD+ave.WD+WHC+ex.mg+ex.ca+ph+dep+ex.k+nit+SWC2+PCT+REMOVAL+THIN_METH+tst+actual.removed,really.big = TRUE,
                             data = cleaner.ov,method="exhaustive")


models.sapling<-regsubsets(sap.Hill~elevation+tri+tpi+roughness+slope+aspect+flowdir+tmin+tmean+tmax+dew+vpdmax+
                             vpdmin+McNab+Bolstad+Profile+Planform+Winds10+Winds50+SWI+RAD+ppt+Parent+Densic+Lithic+
                             Redox+Min_depth+WD+cumulative.WD+ave.WD+WHC+ex.mg+ex.ca+ph+dep+ex.k+nit+SWC2+PCT+REMOVAL+THIN_METH+tst+actual.removed,really.big = TRUE,
                           data = cleaner.sap,method="exhaustive")
#summary(models.ov)
summary(models.overstory)

#overstory VSURF
names(cleaner.ov)
cleaner.ov<-select(cleaner.ov, -c("Northing_Y", "Easting_X", "OtherPlotName", "ThinnedYet", "Notes", "SitePlotID", "PlotArea", "PlotDimensions", "X_Coordinate", "Y_Coordinate", "START_YR", "Densic", "Lithic", "Redox", "Min_depth", "climate"))
names(cleaner.ov)

cleaner.ov<-as.data.frame(cleaner.ov)

cleaner.ov<-cleaner.ov %>% replace(is.na(.), 0)

ov.preds<-cleaner.ov[,6:49]
ov.obs<-cleaner.ov[,5]

ov.vsurf<-VSURF(ov.preds, ov.obs)

ov.vsurf
ov.vsurf$varselect.thres
ov.vsurf$varselect.interp
ov.vsurf$nums.varselect
ov.vsurf$nums.varselect
ov.vsurf$varselect.pred
ov.vars<-ov.preds[c(8,35)] #tmin, tmax, cumulative.WD, ave.WD
names(ov.vars)

v.lm.ov<-lm(ov.Hill~tmin+cumulative.WD, data=cleaner.ov)
summary(v.lm.ov)
check_collinearity(v.lm.ov)

v.ov.model <- lme(ov.Hill~tmin+cumulative.WD,
                   data=cleaner.ov,
                   correlation=corAR1(form=~YEAR|SITEid/PLOTid),
                   random=~1|SITEid/PLOTid,
                   na.action=na.omit,method="REML")
summary(v.ov.model)




##########################################
#one by one
avg.wd.ov<-lm(ov.Hill~tmin+tmax+ave.WD,data=cleaner.ov)
summary(avg.wd.ov) #0.2866

WD.ov<-lm(ov.Hill~tmin+tmax+cumulative.WD, data=cleaner.ov)
summary(WD.ov) #0.2873 winner

v.ov.model <- lme(ov.Hill~tmin+tmax+cumulative.WD,
                   data=cleaner.ov,
                   correlation=corAR1(form=~YEAR|SITEid/PLOTid),
                   random=~1|SITEid/PLOTid,
                   na.action=na.omit,method="REML")

summary(v.ov.model)
plot(v.ov.model)

#res.sum <- summary(models.ov)
#which.max(res.sum$adjr2)
#data.frame(
#  Adj.R2 = which.max(res.sum$adjr2),
#  CP = which.min(res.sum$cp),
#  BIC = which.min(res.sum$bic)
#)



res.sum1 <- summary(models.overstory)
which.max(res.sum1$adjr2)
data.frame(
  Adj.R2 = which.max(res.sum1$adjr2),
  CP = which.min(res.sum1$cp),
  BIC = which.min(res.sum1$bic)
)

model1<-lm(ov.Hill~McNab+Winds10+Parent+Redox+Min_depth+cumulative.WD+ave.WD+ex.ca,data=cleaner.ov)
summary(model1)
check_collinearity(model1)
plot(model1)
AIC(model1)
#mod2 <- lm(Hill~slope+dew+Winds10+Lithic+ave.WD+dep+THIN_METH,data=cleaner)
#check_collinearity(mod2)
#summary(mod2)

#mod2.ov<-lm(ov.Hill~slope+dew+Winds10+Lithic+ave.WD+dep+THIN_METH,data=cleaner.ov)
#check_collinearity(mod2.ov)
#summary(mod2.ov)



#mod3 <- lm(Hill~slope+dew+dep+THIN_METH+RAD+REMOVAL,data=cleaner)
#summary(mod3)
#check_collinearity(mod3)
#plot(mod3)


#names(cleaner)


df <- cleaner[c(5,13,19,29,42,54)]
pairs(df)

AIC(mod2,mod3)
AIC(mod1,mod2)

plot(mod2)

library(nlme)
library(regclass)

check_distribution(cleaner$Hill)
plot(density(cleaned$Hill))
cleaner <- groupedData(Hill~YEAR|SITEid/PLOTid,data=cleaner)
plot(density(cleaner$Hill))
#no.na.data$lithic.dummy <- ifelse(no.na.data$Lithic>0,1,0)
#no.na.data$lithic.depth <- ifelse(no.na.data$Lithic>0,no.na.data$Lithic,NA)
#cleaner2 <- dplyr::filter(cleaner,tst>-1)
cleaner$wdi <- cleaner$WD-cleaner$SWC2
cleaner$run.wdi <- cleaner$run.wd-cleaner$SWC2

full <- lme(Hill~dew+THIN_METH+tst, ## winner
            data=cleaner,
            correlation=corAR1(form=~YEAR|SITEid/PLOTid),
            random=~1|SITEid/PLOTid,
            na.action=na.omit,method="REML")
plot(full)
summary(full)
cleaner$run.wsi <- cleaner$run.wdi*-1
cleaner$run.ws <- cleaner$run.wd*-1

full2 <- lme(Hill~dew+THIN_METH+run.ws+actual.removed, ## winner
             data=cleaner,
             correlation=corAR1(form=~YEAR|SITEid/PLOTid),
             random=~1|SITEid/PLOTid,
             na.action=na.omit,method="REML")

summary(full2)




plot(full2)

names(cleaner)

cd <- cleaner[c(5,18,68,66,53)]
pairs(cd)


full3 <- lme(Hill~dew+THIN_METH+actual.removed+PCT+run.wdi, ## winner
             data=cleaner,
             correlation=corAR1(form=~YEAR|SITEid/PLOTid),
             random=~1|SITEid/PLOTid,
             na.action=na.omit,method="REML")

AIC(full2,full3)



plot(full2)



full4 <- lme(Hill~dew+THIN_METH+run.wdi:actual.removed, ## winner
             data=cleaner,
             correlation=corAR1(form=~YEAR|SITEid/PLOTid),
             random=~1|SITEid/PLOTid,
             na.action=na.omit,method="REML")

summary(full4)
AIC(full,full2,full3,full4)
plot(full2)

rmse(full2)
plot(ranef(full2))
plot(fixef(full2))
summary(full2)
performance(full2)

require(ggeffects)
mydf2<-ggpredict(full2, terms = c("run.ws", "THIN_METH","dew"))

ggplot(mydf2,aes(x=x,y=predicted,colour=group))+
  geom_line(aes(linetype=group,color=group),linewidth=1)+
  labs(linetype="Thinning Method")+
  labs(colour = "Thinning Method")+
  #labs(x="BA Removed (%)",y="Predicted Diversity (Hill)")+
  ylim(0.5,4)+
  #xlim(4,6)+
  facet_wrap(~facet)+
  theme_bw(18) 


summary(full2)
plot(full2)
plot(full2)
rmse(full2)
mae(full2)


cleaner$run.wdi[is.na(cleaner$run.wdi)] <- 0
cleaner <- dplyr::filter(cleaner,run.wdi<0)

cleaner$fit <- predict(full2,cleaner,allow.new.levels=TRUE)
cleaner$resid <- cleaner$Hill-cleaner$fit
rmse(full2)
mae(full2)


plot(cleaner$fit,cleaner$resid)
abline(h=0)
unique(cleaner$THIN_METH)
plot(cleaner$fit,cleaner$Hill,
     ylim=c(0,6),xlim=c(0,6),
     col=factor(cleaner$THIN_METH),
     pch=16,cex=1.5,
     ylab="Observed Hill Diversity",
     xlab="Predicted Hill Diversity")
legend(0.1,6,unique(cleaner$THIN_METH),col=1:length(cleaner$THIN_METH),pch=16)
abline(0,1,col="gray30",lty=3)

require(equivalence)

tost(cleaner$Hill,cleaner$fit,var.equal=FALSE,epsilon=1)
equivalence.xyplot(cleaner$Hill~cleaner$fit,
                  alpha=0.05,b0.ii=0.25,b1.ii=0.25,
                   xlab="Predicted Hill Numbers",
                   ylab="Observed Hill Numbers",
                   xlim=c(0,5),
                   ylim=c(0,5))


##equivalence hard code


#subtract mean prediction from prediction
mean(cleaner$fit) #1.83

cleaner<-cleaner%>%
  mutate(fit2 = cleaner$fit - 1.83)

#establish equivalence interval
#intercept (mean + or - 10%) (1.647-2.013)
#slope (1 + or - 10%) (0.9-1.1)

#fit linear regression (predicted values as independent, actual measurements as dependent)
eq1<-lm(Hill~fit2, data=cleaner)
plot(cleaner$fit2, cleaner$Hill, xlim=c(-1,1), ylim=c(0,5))
     abline(lm(Hill~fit2,data=cleaner), col="red")
     abline(1.80,1.021, lty = 2)
     abline(1.84,1.097, lty = 2)
     abline(v=0, col="grey80")
     abline(h=1.58, col="grey80")
     abline(h=2.08, col="grey80")
     abline(0,0.75, col="lightgreen")
     
confint(eq1)
summary(eq1)
#qqnorm(full)
require(MuMIn)
r.squaredGLMM(full2)

# may want a 3 plot, but for now, we are good. 
# Premer, out, 20240828

###################################################
###################################################
###################################################
no.na <- na.omit(cleaner)

grid.lines <- nrow(no.na)
run.wdi.pred <- seq(min(no.na$run.wdi), max(no.na$run.wdi), length.out = grid.lines)
dew.pred <- seq(min(no.na$dew), max(no.na$dew), length.out = grid.lines)
rem.pred <- seq(min(no.na$actual.removed,max(no.na$actual.removed)),length.out=grid.lines)

xyc.low <- expand.grid(run.wdi=run.wdi.pred,
                       dew=dew.pred,
                       actual.removed=rem.pred,
                       THIN_METH="low")

xyc.cont <- expand.grid(x=run.wdi.pred,
                        y=dew.pred,
                        c=rem.pred,
                        d="control")

xyc.dom <- expand.grid(x=run.wdi.pred,
                       y=dew.pred,
                       c=rem.pred,
                       d="dominant")

z.pred <- matrix(predict(full2, newdata=xyc.low, level=0), 
                 nrow=grid.lines, ncol=grid.lines)



require(plot3D)

persp3D(z=z.pred,theta=40,
        phi = -0,
        expand=1,
        box=T,
        border=NA,
        polygon_offset=3,
        pacakge="rgl",
        zlab="Hill Diversity",
        ylab="BA Removal (%)",
        xlab="WDI",
        cex.lab=1,
        cex.axis=1,
        ticktype="detailed")

plot(z.pred)

fitpoints <- as.vector(predict(test, data))

require(splines)
require(visreg)

install.packages("rgl")
require(rgl)
plot3d(xyc.dom,type="s")


