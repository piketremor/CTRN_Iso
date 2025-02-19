############################################################################
############################################################################
#### Generate predictive models of species diversity from CTRN database (overstory) ####
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
library(performance)
require(MEForLab)
require(lattice)
require(ggeffects)
require(lmtest)

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

#converting WD variables to WS variables
final.over$wdi.time<-final.over$wdi.time * -1
final.over$mean.WD<-final.over$mean.WD * -1
final.over$wd.time<-final.over$wd.time * -1

#renaming the WD variables accordingly 

colnames(final.over)[colnames(final.over) == 'wdi.time'] <- 'wsi.time'
colnames(final.over)[colnames(final.over) == 'wd.time'] <- 'ws.time'
colnames(final.over)[colnames(final.over) == 'mean.WD'] <- 'mean.WS'

# first, let's just run the overstory, then we can run the understory
##
names(final.over)
require(leaps)
models.ov <- regsubsets(Shannon~elevation+tri+tpi+roughness+slope+aspect+flowdir+tmin+tmean+tmax+dew+vpdmax+
                          vpdmin+McNab+Bolstad+Profile+Planform+Winds10+Winds50+SWI+RAD+ppt+Parent+
                          ws.time+wsi.time+mean.WS+WHC+ex.mg+ex.ca+ph+dep+ex.k+nit+SWC2+PCT+THIN_METH+
                          tst+actual.removed+PCT,really.big = TRUE,
                        data = final.over,method="exhaustive")
summary(models.ov)
#tmean,dew,McNab,Planform,wsi.time,WHC,ex.ca,PCT,THIN_METH

#all variables from regsubsets selection
model1<-lm(Shannon~tmean+dew+McNab+Planform+wsi.time+WHC+ex.ca+THIN_METH+PCT,data=final.over)
summary(model1)
check_collinearity(model1) #tmean and dew have VIF above 10
plot(model1)
AIC(model1)

mod2<-lm(Shannon~dew+McNab+Planform+wsi.time+WHC+ex.ca+THIN_METH+PCT,data=final.over)
summary(mod2)
AIC(mod2)
check_collinearity(mod2)

mod3<-lm(Shannon~tmean+McNab+Planform+wsi.time+WHC+ex.ca+THIN_METH+PCT,data=final.over)
summary(mod3)
AIC(mod3)


names(final.over)
df <- final.over[c(9,25,30,40,64,32,51,58,66)]
pairs(df)

#look for at effects of interactions within the model
in.model1<-lm(Shannon~roughness+tmean+dew+Winds50+wdi.time:actual.removed+wd.time:actual.removed+WHC:actual.removed+
                SWC2:actual.removed,data=final.over)
summary(in.model1)
check_collinearity(in.model1)


in.model2<-lm(Shannon~roughness+tmean+dew+Winds50+wdi.time:actual.removed,data=final.over)
summary(in.model2)

in.model3<-lm(Shannon~roughness+tmean+dew+Winds50+wd.time:actual.removed,data=final.over)
summary(in.model3)

in.model4<-lm(Shannon~roughness+tmean+dew+Winds50+WHC:actual.removed,data=final.over)
summary(in.model4)

in.model5<-lm(Shannon~roughness+tmean+dew+Winds50+wdi.time:actual.removed+wdi.time:wd.time,data=final.over)
summary(in.model5)


model6 <- lme(Shannon~roughness+tmean+dew+Winds50+wdi.time:actual.removed+wd.time:actual.removed+wdi.time:wd.time,
              data=final.over,
              correlation=corAR1(form=~YEAR|SITEid/PLOTid),
              random=~1|SITEid/PLOTid,
              na.action=na.omit,method="REML")
summary(model6)
AIC(model6)
model7 <- lme(Shannon~roughness+tmean+dew+Winds50+WHC:actual.removed,
              data=final.over,
              correlation=corAR1(form=~YEAR|SITEid/PLOTid),
              random=~1|SITEid/PLOTid,
              na.action=na.omit,method="REML")
summary(model7)
AIC(model7)

model8<-lme(Shannon~roughness+tmean+dew+Winds50+wdi.time:actual.removed,
            data=final.over,
            correlation=corAR1(form=~YEAR|SITEid/PLOTid),
            random=~1|SITEid/PLOTid,
            na.action=na.omit,method="REML")
summary(model8)

model9<-lme(Shannon~tmean+dew+wd.time:actual.removed,
            data=final.over,
            correlation=corAR1(form=~YEAR|SITEid/PLOTid),
            random=~1|SITEid,
            na.action=na.omit,method="REML") ##winner winner
summary(model9)
check_collinearity(model9)
performance(model9)
plot(model9)



#####
require(ggeffects)
mydf2<-ggpredict(model9, terms = c("wd.time", "actual.removed", "dew"))
mydf3<-predict_response(model9, terms = c("wd.time", "actual.removed", "dew"))

ggplot(mydf3,aes(x=x,y=predicted,colour=group))+
  geom_line(aes(linetype=group,color=group),linewidth=1)+
  labs(linetype="Thinning Method")+
  labs(colour = "Thinning Method")+
  #labs(x="BA Removed (%)",y="Predicted Diversity (Hill)")+
  ylim(0.5,1)+
  #xlim(4,6)+
  facet_wrap(~facet)+
  theme_bw(18) 




#PCT increases AIC slightly and has higher RMSE and lower R2
model10<-lme(Shannon~tmean+dew+wd.time:actual.removed+PCT,
            data=final.over,
            correlation=corAR1(form=~YEAR|SITEid/PLOTid),
            random=~1|SITEid,
            na.action=na.omit,method="REML") 
summary(model10)
check_collinearity(model10)
performance(model10)
plot(model10)


anova(model9,model10)

# model 9 has lower AIC and is preferred through parsimony and better R2 values, RMSE is the same. 
performance(model9)
