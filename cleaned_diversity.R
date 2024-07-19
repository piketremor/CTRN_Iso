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
  reframe(sap.total = sum(spec.count*ef))

all.in <- saplings%>%
  mutate(ef=62.5)%>%
  group_by(SITEid,PLOTid,YEAR)%>%
  summarize(sapling.total = sum(spec.count*ef))

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
over <- filter(over, SITEid == "AS" | SITEid == "DR" | SITEid == "GR" | SITEid == "HR" | SITEid == "KI" | SITEid == "LM" | SITEid == "LT" | SITEid == "PA" | SITEid == "PE" | SITEid == "RC" | SITEid == "RR" | SITEid == "SA" | SITEid == "SC" | SITEid == "SR" | SITEid == "WB") 

over.sum<-over%>%
  mutate(ef = 5)%>%
  group_by(SITEid, PLOTid, YEAR, SPP)%>%
  reframe(over.total = sum(ef))

all.over<-over%>%
  mutate(ef=5)%>%
  group_by(SITEid, PLOTid, YEAR)%>%
  reframe(overstory.total = sum(ef))

overstory<- left_join(over.sum, all.over)
all_tree<-left_join(overstory, sappy, by = c("SITEid", "PLOTid", "SPP", "YEAR"))
all_tree$sap.total[is.na(all_tree$sap.total)] <- 0
all_tree$sapling.total[is.na(all_tree$sapling.total)] <- 0
all_tree$SPP[is.na(all_tree$SPP)]<-"OTHER CONIFER"
all_tree$prop <- (all_tree$over.total+all_tree$sap.total)/(all_tree$overstory.total+all_tree$sapling.total)
all_tree$shann.base <- all_tree$prop*(log(all_tree$prop))
dv <- all_tree%>%
  group_by(SITEid,PLOTid,YEAR)%>%
  reframe(Shannon = sum(shann.base)*-1,
          Hill = exp(Shannon))
head(dv)
# cleaned. 
site.vars <- read.csv("CTRN_SiteVariables_20240718.csv")
site.vars <- dplyr::rename(site.vars,SITEid=SiteID,PLOTid=PlotID)
head(site.vars)
site.dv <- left_join(dv,site.vars)

wd <- read.csv("CTRN_WD_dataframe.csv")
wd <- dplyr::rename(wd,SITEid=SiteID,PLOTid=PlotID,YEAR=year)
site.water <- left_join(site.dv,wd)
head(site.water)
clean.df <- distinct(site.water)

# finally, we need the treatment data. 

treat <- read.csv("CTRN_Enviro.csv")
treat <- treat%>%
  group_by(SITEid,PLOTid)%>%
  top_n(.,1,wt=YEAR)
treat.rx <- treat[c(2,3,10:14)]
cleaned <- left_join(clean.df,treat.rx,by=c("SITEid","PLOTid"))

head(cleaned)
cleaned$tst <- cleaned$YEAR-cleaned$TRT_YR
cleaned$X[is.na(cleaned$X)] <- 0
cleaned <- dplyr::filter(cleaned,X>0)

cleaned$Redox[is.na(cleaned$Redox)] <-0
cleaned$Lithic[is.na(cleaned$Lithic)] <-0
cleaned$Densic[is.na(cleaned$Densic)] <-0
cleaned$Profile[is.na(cleaned$Profile)] <-0
cleaned$Planform[is.na(cleaned$Planform)] <-0
cleaned$FERT <- 0
#cleaned$REMOVAL <- as.factor(cleaned$REMOVAL)
cleaned$PCT <- as.factor(cleaned$PCT)
cleaned$THIN_METH <- as.factor(cleaned$THIN_METH)
head(cleaned)

cleaned$check.yr <- cleaned$TRT_YR-cleaned$YEAR
calibrate <- dplyr::filter(cleaned,check.yr>-1)
calib.set <- calibrate%>%
  group_by(SITEid,PLOTid,YEAR)%>%
  reframe(calib.wd = cumsum(WD))



no.na.data <- na.omit(cleaned)

model.null<- lm(Hill~elevation+tri+tpi+roughness+slope+aspect+flowdir+tmin+tmean+tmax+dew+vpdmax+
                   vpdmin+McNab+Bolstad+Profile+Planform+Winds10+Winds50+SWI+RAD+ppt+Parent+Densic+Lithic+
                   Redox+Min_depth+WD+cumulative.WD+ave.WD+WHC+ex.mg+ex.ca+ph+dep+ex.k+nit+SWC2+PCT+REMOVAL+THIN_METH+tst, data= no.na.data, na.action=na.omit)
summary(model.null)
require(MASS)
mod.step <- model.null%>%
  stepAIC(trace=FALSE)
summary(mod.step)
require(performance)
check_collinearity(mod.step)

mod1 <- lm(Hill~Planform+Lithic+WD+REMOVAL+THIN_METH+tst,data=no.na.data)
summary(mod1)
require(leaps)
models.ov <- regsubsets(Hill~elevation+tri+tpi+roughness+slope+aspect+flowdir+tmin+tmean+tmax+dew+vpdmax+
                          vpdmin+McNab+Bolstad+Profile+Planform+Winds10+Winds50+SWI+RAD+ppt+Parent+Densic+Lithic+
                          Redox+Min_depth+WD+cumulative.WD+ave.WD+WHC+ex.mg+ex.ca+ph+dep+ex.k+nit+SWC2+PCT+REMOVAL+THIN_METH+tst,really.big = TRUE,
                        data = no.na.data,method="exhaustive")
summary(models.ov)
res.sum <- summary(models.ov)
which.max(res.sum$adjr2)
data.frame(
  Adj.R2 = which.max(res.sum$adjr2),
  CP = which.min(res.sum$cp),
  BIC = which.min(res.sum$bic)
)
mod2 <- lm(Hill~tmean+Profile+ppt+Densic+Lithic+THIN_METH,data=no.na.data)
check_collinearity(mod2)
summary(mod2)
AIC(mod1,mod2)
plot(mod2)
library(nlme)
library(regclass)

check_distribution(cleaned$Hill)

cleaned <- groupedData(Hill~1|SITEid/PLOTid,data=cleaned)
plot(density(cleaned$Hill))
full <- lme(Hill~dew+Lithic+THIN_METH,
            data=cleaned,
            correlation=corAR1(form=~YEAR|SITEid/PLOTid),
            random=~1|SITEid/PLOTid,
            na.action=na.omit)
summary(full)
performance(full)
plot(full)


require(ggeffects)
mydf2<-ggpredict(full, terms = c("Lithic", "THIN_METH", "dew"))

ggplot(mydf2,aes(x=x,y=predicted,colour=group))+
  geom_line(aes(linetype=group,color=group),size=1)+
  labs(linetype="Thinning Method")+
  labs(colour = "Thinning Method")+
  #labs(x="Basal Area per Acre",y="Predicted Diversity (Hill)")+
  #ylim(0,55)+
  #xlim(4,6)+
  facet_wrap(~facet)+
  theme_bw(18) 

xyplot(Hill~WD|SITEid,data=cleaned,
       type=c("p","smooth"))

plot(full)
rmse(full)

no.na.data$fit <- predict(full,no.na.data)
mae(no.na.data$fit,no.na.data$Hill)
no.na.data$resid <- no.na.data$Hill-no.na.data$fit
plot(no.na.data$fit,no.na.data$resid)
abline(h=0)
plot(no.na.data$fit,no.na.data$Hill,
     ylim=c(0,6),xlim=c(0,6))
abline(0,1)



#qqnorm(full)
require(MuMIn)
r.squaredGLMM(full)



