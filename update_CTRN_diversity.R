############################################################################
############################################################################
#### Generate predictive models of species diversity from CTRN database ####
############################################################################
############################################################################

# 2024 08 28 - Code originated from Lila Beck and Mike Premer
# updated 2024 11 07

#load packages
dev.off()
rm(list=ls())
library(dplyr)
library(nlme)
library(ggplot2)
library(leaps)
library(VSURF)
require(MEForLab)
require(lattice)
require(ggeffects)
require(lmtest)

#setwd("~/Google Drive//My Drive/CTRN_CFRU_Share/raw/csv")
setwd("~/Google Drive/My Drive/Research/CFRU/CTRN_CFRU_Share/raw/csv")
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

# first, let's just run the overstory, then we can run the understory
##
names(final.over)
require(leaps)
models.ov <- regsubsets(over.Hill~elevation+tri+tpi+roughness+slope+aspect+flowdir+tmin+tmean+tmax+dew+vpdmax+
                          vpdmin+McNab+Bolstad+Profile+Planform+Winds10+Winds50+SWI+RAD+ppt+Parent+
                          wd.time+wdi.time+mean.WD+WHC+ex.mg+ex.ca+ph+dep+ex.k+nit+SWC2+PCT+THIN_METH+
                          tst+actual.removed,really.big = TRUE,
                        data = final.over,method="exhaustive")
res.sum1 <- summary(models.ov)
which.max(res.sum1$adjr2)
data.frame(Adj.R2 = which.max(res.sum1$adjr2),CP = which.min(res.sum1$cp),BIC = which.min(res.sum1$bic))
res.sum1

# regsubsets suggests that roughness, tmean, dew, McNab, wd.time, wdi.time, WHC, SWC2, ThinMeth, all good indicators, but mind you this is for the entire dataset
mod1 <- lm(over.Hill~tmean+dew+wdi.time*actual.removed+THIN_METH,
           data=final.over)
check_collinearity(mod1)
summary(mod1)
acf(mod1$residuals)
dwtest(mod1)
final.over$mean.wdi <- final.over$mean.WD-final.over$SWC2
mod2 <- lme(Shannon~wdi.time*actual.removed+tmean+dew+THIN_METH,
            data=final.over,
            random=~1|SITEid,
            na.action="na.omit",
            method="REML",
            #weights=varPower(0.2,form=~wdi.time),
            correlation = corAR1(form=~1|SITEid))
qqnorm(mod2, ~ resid(., type = "p") | SITEid, abline = c(0, 1))
performance(mod2)
acf(mod2$residuals)
summary(mod2)
plot(mod2)
mydf2<-ggpredict(mod2, terms = c("actual.removed", "THIN_METH","dew"))
png(filename="Overstory_diversity_CTRN.png",width=18,height=7,units="in",res=800)
ggplot(mydf2,aes(x=x,y=exp(predicted),colour=group))+
  geom_line(aes(linetype=group,color=group),linewidth=2)+
  labs(linetype="Thinning Method")+
  labs(colour = "Thinning Method")+
  labs(x="RD Removed (%)",y="Predicted Diversity (Hill)")+
  ylim(0,5)+
  #xlim(4,6)+
  facet_wrap(~facet)+
  theme_bw(22) 
dev.off()
plot((mod2))
performance(mod2)



