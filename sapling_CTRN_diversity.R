#########################################################################################
#########################################################################################
#### Generate predictive models of species diversity from CTRN database (understory) ####
#########################################################################################
#########################################################################################

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
library(performance)
require(MEForLab)
require(lattice)
require(ggeffects)
require(lmtest)
library(glmmTMB)

##################Clean data and generate overstory explainatory variables##############
setwd("~/Google Drive//My Drive/CTRN_CFRU_Share/raw/csv")
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

###############understory time###########################
saplings <- read.csv("Saplings.csv")

saplings <- filter(saplings, SITEid == "AS" | SITEid == "DR" | SITEid == "GR" | SITEid == "HR" | SITEid == "KI" | SITEid == "LM" | SITEid == "LT" | SITEid == "PA" | SITEid == "PE" | SITEid == "RC" | SITEid == "RR" | SITEid == "SA" | SITEid == "SC" | SITEid == "SR" | SITEid == "WB") 
saplings[saplings == "SpecAld"]<-"SA"
saplings[saplings == "HM"]<-"EH"
saplings[saplings == "CH"]<-"BC"
saplings[saplings == "OT (red oak)"] <- "RO"
saplings$X1.2.inch[is.na(saplings$X1.2.inch)] <- 0
saplings$X1.inch[is.na(saplings$X1.inch)] <- 0
saplings$X2.inch[is.na(saplings$X2.inch)] <- 0
saplings<-filter(saplings, SPP!="NONE")
#diversity 
saplings$spec.count <- (saplings$X1.2.inch+saplings$X1.inch+saplings$X2.inch)
sap.sum <- saplings%>%
  mutate(ef=250)%>%
  group_by(SITEid,PLOTid,CORNERid,YEAR,SPP)%>%
  reframe(sap.spp.total = sum(spec.count*ef))
sap.sum
all.in <- saplings%>%
  mutate(ef=250)%>%
  group_by(SITEid,PLOTid,CORNERid,YEAR)%>%
  summarize(sapling.total = sum(spec.count*ef))
all.in
sappy <- left_join(all.in,sap.sum)
sappy
sappy$sap.spp.total <- ifelse(sappy$sap.spp.total>sappy$sapling.total,sappy$sapling.total,sappy$sap.spp.total)
sappy$sap.prop<-(sappy$sap.spp.total/sappy$sapling.total)
sappy
unique(sappy$SPP)
sappy$sap.shann.base<-sappy$sap.prop*(log(sappy$sap.prop))
sappy
sap.dv <- sappy%>%
  group_by(SITEid,PLOTid,CORNERid,YEAR)%>%
  summarise(sap.Shannon = sum(sap.shann.base)*-1)
sap.dv
forest<-left_join(sap.dv,final.over)
head(forest)
##
names(forest)
require(leaps)

forest$shann.dummy <- ifelse(forest,sap.Shannon>0,1)


run.df <- filter(forest,sap.Shannon>0)
unique(zero.df$sap.Shannon)
models.sap <- regsubsets(sap.Shannon~elevation+tri+tpi+roughness+slope+aspect+flowdir+tmin+tmean+tmax+dew+vpdmax+
                          vpdmin+McNab+Bolstad+Profile+Planform+Winds10+Winds50+SWI+RAD+ppt+Parent+
                          wd.time+wdi.time+mean.WD+WHC+ex.mg+ex.ca+ph+dep+ex.k+nit+SWC2+PCT+THIN_METH+
                          tst+actual.removed+bapa+tpa+qmd+RD+CCF+Shannon+over.Hill+ht40+prop.ws.avg+prop.bs.avg+
                          prop.rs.avg+prop.hw.avg,really.big = TRUE,data = run.df,method="exhaustive")


# for the continuous data, it looks like McNab, SWI, ex.ca, ph, Thin meth, tst, qmd, over.Hill, prop.hw.avg

# Vsurf for the binary outcome
require(VSURF)
require(randomForest)

forest$sap.dum <- ifelse(forest$sap.Shannon>0,1,0)
forest$sap.dum[is.na(forest$sap.dum)] <- 999
bin.df <- dplyr::filter(forest,sap.dum<998)

resp <- as.factor(bin.df$sap.dum)
obs <- as.data.frame(bin.df[c(6:12,14:21,24:29,31:45,53:60,64,67)])
obs <- na.roughfix(obs)



m <- VSURF(obs,resp)

m$varselect.pred
names(obs)

#ht40 and tmean for the zero inflation. 

summary(models.sap)
model1<-lm(sap.Shannon~roughness+Winds10+wd.time+wdi.time+ph+tst+actual.removed+tpa+CCF,data=forest)
summary(model1)
check_collinearity(model1)
plot(model1)


forest$bapa[is.na(forest$bapa)]<-0
forest<-filter(forest, bapa>0)
forest$actual.removed[is.na(forest$actual.removed)]<-0

sap.mod <- glmmTMB(sap.Shannon~(1|SITEid/PLOTid/CORNERid)+THIN_METH:wdi.time,
                   data=forest,
                   ziformula=~wdi.time+THIN_METH+CCF+tst+tmean+ht40,
                   family=ziGamma(link="log"),
                   na.action="na.omit") #best model? need to include autocorrelation

require(DHARMa)
res = simulateResiduals(sap.mod)
plot(res)
res = recalculateResiduals(res, group = forest$YEAR)
testTemporalAutocorrelation(res, time = unique(forest$YEAR))

summary(sap.mod)

# don't need to adjust for temporal autocorrelation according to DW test, but the variables could be beefed up
# for the continuous data, it looks like McNab, SWI, ex.ca, ph, Thin meth, tst, qmd, over.Hill, prop.hw.avg

sap.mod2 <- glmmTMB(sap.Shannon~(1|SITEid/PLOTid/CORNERid)+THIN_METH:wdi.time+over.Hill+prop.hw.avg,
                   data=forest,
                   ziformula=~wdi.time+THIN_METH+CCF+tst+tmean+ht40,
                   family=ziGamma(link="log"),
                   na.action="na.omit") #best model? need to include autocorrelation
summary(sap.mod2)
plot(sap.mod2)
performance(sap.mod2)

MuMIn::r.squaredGLMM(sap.mod2)

mydf2 <- ggpredict(sap.mod2,terms=c("tst","THIN_METH","CCF"),type="zi.prob",allow.new.levels=TRUE)


#png("~/Desktop/SMC_Sinuosity_Model_Output.png",units='in',height=5.5,width=14,res=1000)
#theme_set(theme_bw(16))

library(ggplot2)
ggplot(mydf2,aes(x=x,y=predicted,colour=group))+
  geom_line(aes(linetype=group,color=group),size=1)+
  labs(x="Time since treatment (years)",y="Probabilty of no diversity (%)")+
  #ylim(0,55)+
  #xlim(4,6)+
  facet_wrap(~facet)+
  theme_bw(12) +
  #theme(legend.position="none")+
  scale_color_manual(values=c('gray0','gray50','gray25',"gray75"))
  #scale_fill_manual(values=c('gray0','gray70','gray40'), name="fill")

performance(sap.mod)
performance(sap.mod2)

check_collinearity(sap.mod2)

anova(sap.mod,sap.mod2)

plot(residuals(sap.mod2))
plot(acf(sap.mod))
sap.mod
resid<-residuals(sap.mod)
summary(sap.mod2)

plot(resid)





