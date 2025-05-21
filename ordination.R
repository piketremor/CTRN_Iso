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
library(gridExtra)
require(MEForLab)

setwd("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv")
setwd("~/Google Drive/My Drive/Research/CFRU/CTRN_CFRU_Share/raw/csv")

saplings <- read.csv("Saplings.csv")[,2:9]

saplings <- filter(saplings, SITEid == "AS" | SITEid == "DR" | SITEid == "GR" | SITEid == "HR" | SITEid == "KI" | SITEid == "LM" | SITEid == "LT" | SITEid == "PA" | SITEid == "PE" | SITEid == "RC" | SITEid == "RR" | SITEid == "SA" | SITEid == "SC" | SITEid == "SR" | SITEid == "WB") 
saplings[saplings == "SpecAld"]<-"SA"
saplings[saplings == "HM"]<-"EH"
saplings[saplings == "CH"]<-"BC"
saplings[saplings == "NC"]<-"WC"

detach(package:plyr)
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
require(dplyr)
#detach(package:plyr)
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

over.standlist <- left_join(overstory.summary,shannon)
over.standlist <- left_join(over.standlist,ht.40)
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

#converting WD variables to WS variables
final.over$wdi.time<-final.over$wdi.time * -1
final.over$mean.WD<-final.over$mean.WD * -1
final.over$wd.time<-final.over$wd.time * -1

#renaming the WD variables accordingly 

colnames(final.over)[colnames(final.over) == 'wdi.time'] <- 'wsi.time'
colnames(final.over)[colnames(final.over) == 'wd.time'] <- 'ws.time'
colnames(final.over)[colnames(final.over) == 'mean.WD'] <- 'mean.WS'


sunames(final.over)
final.over$sapID<-paste0(final.over$SITEid,"-",final.over$PLOTid)
#change from integer to numeric
final.over$flowdir<-as.numeric(final.over$flowdir)
names(final.over)

head(final.over)



##2018 Ordination and Data Cleaning

#Filter for just year 2018
saplings18<-filter(saplings, YEAR == 2018)
unique(saplings18$SPP)

saplings$uid <- as.data.frame(paste0(saplings$SITEid,".",saplings$PLOTid,".",saplings$CORNERid))
template <- data.frame(uid=unique(saplings$uid))
names(template)
template$uid <- (template$paste0.saplings.SITEid.......saplings.PLOTid.......saplings.CORNERid.)
names(template)
templ <- data.frame(uid=template$uid)
saplings18$uid <- as.factor(paste0(saplings18$SITEid,".",saplings18$PLOTid,".",saplings18$CORNERid))
saplings.df <- left_join(templ,saplings18)



#filter for species >5%
branch <-filter(saplings.df, SPP =="PB"|SPP == "OT"|SPP=="WC"|SPP=="RM"|SPP=="RS"|SPP=="WP"|SPP=="BF"|SPP=="YB"|SPP=="EH"|SPP=="QA"|SPP=="BC")

unique(branch$SPP)
#recalculate iv values

# because there are some corners with 0s, we cannot use the corner id, so we will have to combine regen plots



# 12/19/2024
# NEED TO TURN NA'S IN THE 1.2, 1, AND 2 INCH COUNTS TO ZEROS, OTHERWISE SAP.BA RETURNS AN NA. 
# WILL WORK ON TOMORROW (12/20)

head(branch)
branch$X1.2.inch[is.na(branch$X1.2.inch)] <- 0
branch$X1.inch[is.na(branch$X1.inch)] <- 0
branch$X2.inch[is.na(branch$X2.inch)] <- 0


head(branch)

brancher <- branch%>%
  mutate(SAP_EXP = 62.5,
         bapa.half = ((0.5^2*0.005454)*X1.2.inch)*SAP_EXP,
         bapa.one = (1.0^2*0.005454)*X1.inch*SAP_EXP,
         bapa.two = (2.0^2)*0.005454*X2.inch*SAP_EXP)

brancher$sap.bapa <- brancher$bapa.half+brancher$bapa.one+brancher$bapa.two
head(brancher)


plot_summary_sap <- brancher%>%
  group_by(SITEid, PLOTid,YEAR)%>%
  summarise(SAP_TPA_total = sum(SAP_EXP),
            SAPBA_total = sum(sap.bapa))

branched <- brancher%>%
  group_by(SITEid, PLOTid,YEAR, SPP)%>%
  summarise(SAP_TPA = sum(SAP_EXP),
            SAPBA = sum(sap.bapa))

branches <- branched%>%
  left_join(.,plot_summary_sap)%>%
  mutate(prop_tpa = (SAP_TPA/SAP_TPA_total),
         prop_ba = (SAPBA/SAPBA_total),
         iv = ((prop_tpa + prop_ba)/2))
#run linear model with IV from RM and BF
Red<-filter(branches,SPP=="RM")%>%
  left_join(.,final.over)
BF<-filter(branches,SPP=="BF")%>%
  left_join(.,final.over)

RMsub<-regsubsets(iv~prop.rs.avg+prop.hw.avg+prop.ab.avg+prop.eh.avg+prop.bf.avg+prop.ws.avg+prop.bs.avg+tpi+
          roughness+slope+aspect+flowdir+tmin+over.Hill+ht40+elevation+tri+tmean+tmax+dew+vpdmax+vpdmin+McNab+Bolstad+
          Profile+Planform+Winds10+Winds50+SWI+RAD+ppt+WHC+ex.mg+ex.ca+ph+dep+ex.k+nit+SWC2+PCT+THIN_METH+mean.WS+
          tst+ws.time+actual.removed+wsi.time,really.big=TRUE,data = Red,method="exhaustive")
summary(RMsub)
rmod1<-lm(iv~tst+Densic+actual.removed,data=Red)
summary(rmod1)
rmod2<-lm(iv~tst+actual.removed,data=Red)
summary(rmod2)
plot(rmod2)
rmod3<-lm(iv~prop.eh.avg+slope+flowdir+tri+tmax+Densic+ex.mg+tst+actual.removed,data=Red)
summary(rmod3)
rmod4<-lm(iv~tst,data=Red)
summary(rmod4)
BFsub<-regsubsets(iv~prop.rs.avg+prop.hw.avg+prop.ab.avg+prop.eh.avg+prop.bf.avg+prop.ws.avg+prop.bs.avg+tpi+
                    roughness+slope+aspect+flowdir+tmin+over.Hill+ht40+elevation+tri+tmean+tmax+dew+vpdmax+vpdmin+McNab+Bolstad+
                    Profile+Planform+Winds10+Winds50+SWI+RAD+ppt+WHC+ex.mg+ex.ca+ph+dep+ex.k+nit+SWC2+PCT+THIN_METH+mean.WS+
                    tst+ws.time+actual.removed+wsi.time,really.big=TRUE,data = BF,method="exhaustive")
summary(BFsub)
bfmod1<-lm(iv~WHC,data=BF)
summary(bfmod1)
bfmod2<-lm(iv~prop.ws.avg+roughness+ht40+vpdmin+McNab+Winds50+ppt+WHC+mean.WS,data=BF)
summary(bfmod2)
bfmod2<-lm(iv~prop.ws.avg+McNab+Winds50+ppt+WHC+mean.WS,data=BF)
summary(bfmod2)
bfmod3<-lm(iv~prop.ws.avg+McNab+Winds50+mean.WS,data=BF)
summary(bfmod3)
check_collinearity(bfmod3)
#create site-plot identifier
branches$sapID<-paste0(branches$SITEid,"-",branches$PLOTid)

#get data in long format
molten <- melt(as.data.frame(branches),id=c("sapID","iv","SPP"))
sapwide <- dcast(molten,sapID~SPP,value.var = "iv",mean)
sapwide[is.na(sapwide)] <- 0
head(sapwide)
# double check they all add to 1
rowSums(sapwide[2:12])

#NMDS
set.seed(123)
sapordi<-metaMDS(sapwide[,2:12], distance = "bray",k=3)
stressplot(sapordi)
sapordi
plot(sapordi)



plot(sapordi,type="n")
points(sapordi,display="sites",cex=2,pch=21,col="red", bg="yellow")
text(sapordi,display="spec",cex=1.5,col="red")


###Keep this section, code for pulling out species and site correlations with NMDS axes###
#ce <- cor(sapwide,method="pearson",scores(sapordi,dis="si"))
#try.x <- (ce*sqrt(96-2)/
#            sqrt(1-ce^2))
#try.x
#dt(try.x,df=95)

#cannot run categorical variables with correlations
#cx <- cor(re.env,method="pearson",scores(sapordi,dis="si"))
#try <- (cx*sqrt(96-2)/
#          sqrt(1-cx^2))
#try
#try.p<-dt(try,df=95)

########################

## Constrained Ordination###
#join together sapling and environmental data to make sure the number of observations match
final.over <- filter(final.over,YEAR=="2018")
unique(final.over$YEAR)
#final.over<-na.omit(final.over)
bark<-left_join(sapwide, final.over)
#removing observations that are missing THIN_METH data because can't replace NA w 0 due to factor 
barky<-filter(bark, THIN_METH == "control" | THIN_METH == "crown" | THIN_METH == "low" | THIN_METH == "dominant")
unique(barky$THIN_METH)
names<-barky$sapID #create list of sapID to set the rownames with 


#separate the datasets back out
sapwide<-barky[,c(2:12)]
env<-barky[c(16:20,21:52,54:61,63,65:69)]
env[is.na(env)] <- 0

#environmental constraints
ef <- envfit(sapwide, env,na.rm=TRUE)
plot(ef)
sap.rda <- rda(sapwide ~ ., data=env, na.action = na.exclude) 
summary(sap.rda)
ordiplot(sap.rda, scaling = 2, type = "text")
RsquareAdj(sap.rda)$adj.r.squared



mod0 <- rda(sapwide~1,env)
mod1 <- rda(sapwide~.,env,na.action=na.exclude)
set.seed(123)
mod <- ordistep(mod0, 
                scope=formula(mod1),
                direction="both")
#wsi.time, WHC, tst, ppt,ht40,dep,ex.k,ex.mg,ws.time,RD,vpdmax


plot(mod,add=TRUE,type="t")
anova(mod,type="t")
summary(mod)
anova(mod, by = "margin")
anova(mod, by="term")
two_env<-select(env, c(wsi.time, WHC, tst, ht40,dep,ex.k,ex.mg,ws.time,RD,vpdmax))
vif.cca(mod)
#wsi.time vs ws.time
test1<-rda(sapwide ~ wsi.time+WHC+tst+ppt+ht40+dep+ex.k+ex.mg+RD+vpdmax,env)
RsquareAdj(test1)$adj.r.squared

test2<-rda(sapwide ~ ws.time+WHC+tst+ppt+ht40+dep+ex.k+ex.mg+RD+vpdmax,env)
RsquareAdj(test2)$adj.r.squared
#wsi.time is like, slightly better

#wsi.time vs tst?
test3<-rda(sapwide ~ wsi.time+WHC+ppt+ht40+dep+ex.k+ex.mg+RD+vpdmax,env)
RsquareAdj(test3)$adj.r.squared

test4<-rda(sapwide ~ WHC+tst+ppt+ht40+dep+ex.k+ex.mg+RD+vpdmax,env)
RsquareAdj(test4)$adj.r.squared

#final variables: wsi.time, WHC, ppt, ht40, dep, ex.k, ex.mg, RD, vpdmax

two_env<-select(env, c(wsi.time, WHC, ppt,ht40,dep,ex.k,ex.mg,RD,vpdmax))
two_env$HT100<-env$ht40*0.3
mod2<-rda(sapwide~.,two_env,na.action=na.exclude)
RsquareAdj(mod2)$adj.r.squared
anova(mod2,by="term")
#vpdmax is not significant, drop?? 
mod3<-rda(sapwide ~ wsi.time+WHC+ppt+ht40+dep+ex.k+ex.mg+RD,env)
RsquareAdj(mod3)$adj.r.squared #worse dropping doesn't really improve R squared value

vif.cca(mod2)
summary(mod2)
ordiplot(mod2, scaling = 2, type = "text")
ordiplot(mod2, scaling = 2, type = "none")%>%
  points("sites", pch=16)%>%
  text("species", col="red")%>%
  text(what="biplot", col="blue")
set.seed(123)
sappy<-metaMDS(sapwide, distance = "bray",k=3)
plot(sappy,type="n")
points(sappy,display="sites",cex=2,pch=21,col="red", bg="yellow")
text(sappy,display="spec",cex=1.5,col="red")
ordisurf(sappy, two_env$wsi.time, add = TRUE, family = quasipoisson)
en<-envfit(sappy, two_env[,c(2:3,5:10)], add=TRUE)
plot(en)


text(sappy,display="spec",cex=1.5,col="red")
points(sappy,display="sites",cex=2,pch=21,col="red", bg="yellow")
en<-envfit(sappy, two_env, add=TRUE)
plot(en)

#3D plot of sapling NMDS
ordiplot3d(sappy)
text(sappy, display="spec",cex=1.5,col="red")
plot(en)

sappy

#species and constraint correlations with NMDS axes
ce <- cor(sapwide,method="pearson",scores(sappy,dis="si"))
try.x <- (ce*sqrt(96-2)/
            sqrt(1-ce^2))
try.x
try.p<-dt(try.x,df=95)

#write.csv(try.x, "~/Desktop/updatedsapling_species_correlations.csv")
#write.csv(try.p, "~/Desktop/updatedsapling_species_pval.csv")

#cannot run categorical variables with correlations
cx <- cor(two_env,method="pearson",scores(sappy,dis="si"))
try <- (cx*sqrt(96-2)/
          sqrt(1-cx^2))
try
try.p2<-dt(try,df=95)
try.p2

#write.csv(try, "~/Desktop/updatedsapling_env_cor.csv")
#write.csv(try.p2, "~/Desktop/updatedsapling_env_pval.csv")

## Premer end. 

#######MRPP######

trt<-env$THIN_METH
trt<-na.omit(trt)
dis <- vegdist(sapwide)
mod <- betadisper(dis, trt)
mod

#multivariate dispersion
test<-mrpp(dat=dis, grouping=trt, permutations=999)
test

test
centroids<-data.frame(grps=rownames(mod$centroids),data.frame(mod$centroids))
vectors<-data.frame(group=mod$group,data.frame(mod$vectors))
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")
grp1.hull<-seg.data[seg.data$group=="control",1:3][chull(seg.data[seg.data$group=="control",2:3]),]
grp2.hull<-seg.data[seg.data$group=="crown",1:3][chull(seg.data[seg.data$group=="crown",2:3]),]
grp3.hull<-seg.data[seg.data$group=="dominant",1:3][chull(seg.data[seg.data$group=="dominant",2:3]),]
grp4.hull<-seg.data[seg.data$group=="low",1:3][chull(seg.data[seg.data$group=="low",2:3]),]
all.hull<-rbind(grp1.hull,grp2.hull,grp3.hull,grp4.hull)

panel.a<-ggplot() + 
  geom_point(data=centroids[1,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=16) + 
  geom_point(data=seg.data[seg.data$group %in% "control",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=16) +
  labs(title="Control",x="",y="") +
  coord_cartesian(xlim = c(-0.5,0.7), ylim = c(-0.25,0.2)) +
  theme_classic() + 
  theme(legend.position="none")

panel.b<-ggplot() + 
  geom_point(data=centroids[2,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=17) + 
  geom_point(data=seg.data[seg.data$group %in% "crown",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=17) +
  labs(title="Crown",x="",y="") +
  coord_cartesian(xlim = c(-0.5,0.7), ylim = c(-0.25,0.2)) +
  theme_classic() + 
  theme(legend.position="none")

panel.c<-ggplot() + 
  geom_point(data=centroids[3,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=15) + 
  geom_point(data=seg.data[seg.data$group %in% "dominant",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=15) +
  labs(title="Dominant",x="",y="") +
  coord_cartesian(xlim = c(-0.1,0.7), ylim = c(-0.25,0.2)) +
  theme_classic() + 
  theme(legend.position="none")

panel.d<-ggplot() + 
  geom_point(data=centroids[3,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=15) + 
  geom_point(data=seg.data[seg.data$group %in% "low",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=15) +
  labs(title="Low",x="",y="") +
  coord_cartesian(xlim = c(-0.1,0.7), ylim = c(-0.25,0.2)) +
  theme_classic() + 
  theme(legend.position="none")

panel.e<-ggplot() + 
  geom_point(data=centroids[,1:3], aes(x=PCoA1,y=PCoA2,shape=grps),size=4,colour="red") + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2,shape=group),size=2) +
  labs(title="All",x="",y="") +
  coord_cartesian(xlim = c(-0.5,0.7), ylim = c(-0.25,0.2)) +
  theme_classic() + 
  theme(legend.position="none")

grid.arrange(panel.a,panel.b,panel.c,panel.d, panel.e, nrow=2)

#add vectors and hulls
panel.a<-ggplot() +
  geom_polygon(data=all.hull[all.hull=="control",],aes(x=v.PCoA1,y=v.PCoA2),colour="blue",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data[seg.data$group %in% "control",],aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) + 
  geom_point(data=centroids[1,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=16) + 
  geom_point(data=seg.data[seg.data$group %in% "control",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=16,colour="blue") +
  labs(title="Control",x="",y="") +
  coord_cartesian(xlim = c(-0.7,0.7), ylim = c(-0.7,0.7)) +
  theme_classic() + 
  theme(legend.position="none")

panel.b<-ggplot() +
  geom_polygon(data=all.hull[all.hull=="crown",],aes(x=v.PCoA1,y=v.PCoA2),colour="orange",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data[seg.data$group %in% "crown",],aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) + 
  geom_point(data=centroids[2,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=17) + 
  geom_point(data=seg.data[seg.data$group %in% "crown",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=17,colour="orange") +
  labs(title="Crown",x="",y="") +
  coord_cartesian(xlim = c(-0.7,0.7), ylim = c(-0.7,0.7)) +
  theme_classic() + 
  theme(legend.position="none")

panel.c<-ggplot() +   
  geom_polygon(data=all.hull[all.hull=="dominant",],aes(x=v.PCoA1,y=v.PCoA2),colour="green",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data[seg.data$group %in% "dominant",],aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) +
  geom_point(data=centroids[3,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=15) + 
  geom_point(data=seg.data[seg.data$group %in% "dominant",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=15,colour="green") + 
  labs(title="Dominant",x="",y="") +
  coord_cartesian(xlim = c(-0.5,0.7), ylim = c(-0.7,0.7)) +
  theme_classic() + 
  theme(legend.position="none")


panel.d<-ggplot() +
  geom_polygon(data=all.hull[all.hull=="low",],aes(x=v.PCoA1,y=v.PCoA2),colour="purple",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data[seg.data$group %in% "low",],aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) +
  geom_point(data=centroids[3,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=15) + 
  geom_point(data=seg.data[seg.data$group %in% "low",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=15,color="purple") +
  labs(title="Low",x="",y="") +
  coord_cartesian(xlim = c(-0.5,0.7), ylim = c(-0.7,0.7)) +
  theme_classic() + 
  theme(legend.position="none")

panel.e<-ggplot() + 
  geom_polygon(data=all.hull[all.hull=="control",],aes(x=v.PCoA1,y=v.PCoA2),colour="blue",alpha=0,linetype="dashed") +
  geom_polygon(data=all.hull[all.hull=="crown",],aes(x=v.PCoA1,y=v.PCoA2),colour="orange",alpha=0,linetype="dashed") +
  geom_polygon(data=all.hull[all.hull=="dominant",],aes(x=v.PCoA1,y=v.PCoA2),colour="green",alpha=0,linetype="dashed") +
  geom_polygon(data=all.hull[all.hull=="low",],aes(x=v.PCoA1,y=v.PCoA2),colour="purple",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data,aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) + 
  geom_point(data=centroids[,1:3], aes(x=PCoA1,y=PCoA2,shape=grps),size=4,colour="red") + 
  geom_point(data=seg.data[seg.data$group %in% "control",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=16,colour="blue") +
  geom_point(data=seg.data[seg.data$group %in% "crown",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=17,colour="orange") +
  geom_point(data=seg.data[seg.data$group %in% "dominant",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=15,colour="green") +
  geom_point(data=seg.data[seg.data$group %in% "low",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=15,color="purple") +
  labs(title="All",x="",y="") +
  coord_cartesian(xlim = c(-0.7,0.7), ylim = c(-0.7,0.7)) +
  theme_classic() + 
  theme(legend.position="none")

grid.arrange(panel.a,panel.b,panel.c,panel.d, nrow=2)




#####MRPP comparing overstory and understory########################################################

#first bring in sapling data
setwd("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv")


saplings <- read.csv("Saplings.csv")[,2:9]

saplings <- filter(saplings, SITEid == "AS" | SITEid == "DR" | SITEid == "GR" | SITEid == "HR" | SITEid == "KI" | SITEid == "LM" | SITEid == "LT" | SITEid == "PA" | SITEid == "PE" | SITEid == "RC" | SITEid == "RR" | SITEid == "SA" | SITEid == "SC" | SITEid == "SR" | SITEid == "WB") 
saplings[saplings == "SpecAld"]<-"SA"
saplings[saplings == "HM"]<-"EH"
saplings[saplings == "CH"]<-"BC"
saplings[saplings == "NC"]<-"WC"


#Filter for just year 2018
saplings18<-filter(saplings, YEAR == 2018)
unique(saplings18$SPP)

saplings$uid <- as.data.frame(paste0(saplings$SITEid,".",saplings$PLOTid,".",saplings$CORNERid))
template <- data.frame(uid=unique(saplings$uid))
names(template)
template$uid <- (template$paste0.saplings.SITEid.......saplings.PLOTid.......saplings.CORNERid.)
names(template)
templ <- data.frame(uid=template$uid)
saplings18$uid <- as.factor(paste0(saplings18$SITEid,".",saplings18$PLOTid,".",saplings18$CORNERid))
saplings.df <- left_join(templ,saplings18)



#filter for species >5%
branch <-filter(saplings.df, SPP =="PB"|SPP == "OT"|SPP=="WC"|SPP=="RM"|SPP=="RS"|SPP=="WP"|SPP=="BF"|SPP=="YB"|SPP=="EH"|SPP=="QA"|SPP=="BC")

unique(branch$SPP)

head(branch)
branch$X1.2.inch[is.na(branch$X1.2.inch)] <- 0
branch$X1.inch[is.na(branch$X1.inch)] <- 0
branch$X2.inch[is.na(branch$X2.inch)] <- 0


brancher <- branch%>%
  mutate(SAP_EXP = 62.5,
         bapa.half = ((0.5^2*0.005454)*X1.2.inch)*SAP_EXP,
         bapa.one = (1.0^2*0.005454)*X1.inch*SAP_EXP,
         bapa.two = (2.0^2)*0.005454*X2.inch*SAP_EXP)

brancher$sap.bapa <- brancher$bapa.half+brancher$bapa.one+brancher$bapa.two
head(brancher)


plot_summary_sap <- brancher%>%
  group_by(SITEid, PLOTid,YEAR)%>%
  summarise(SAP_TPA_total = sum(SAP_EXP),
            SAPBA_total = sum(sap.bapa))

branched <- brancher%>%
  group_by(SITEid, PLOTid,YEAR, SPP)%>%
  summarise(SAP_TPA = sum(SAP_EXP),
            SAPBA = sum(sap.bapa))

branches <- branched%>%
  left_join(.,plot_summary_sap)%>%
  mutate(prop_tpa = (SAP_TPA/SAP_TPA_total),
         prop_ba = (SAPBA/SAPBA_total),
         iv = ((prop_tpa + prop_ba)/2))


#create site-plot identifier
branches$sapID<-paste0(branches$SITEid,"-",branches$PLOTid)

#get data in long format
molten <- melt(as.data.frame(branches),id=c("sapID","iv","SPP"))
sapwide <- dcast(molten,sapID~SPP,value.var = "iv",mean)
sapwide[is.na(sapwide)] <- 0
head(sapwide)
# double check they all add to 1
rowSums(sapwide[2:12])
names<-sapwide$sapID
rownames(sapwide)<-names
new<-data.frame(sapID=sapwide$sapID,BC=sapwide$BC,BF=sapwide$BF,EH=sapwide$EH,WC=sapwide$WC,OT=sapwide$OT,PB=sapwide$PB,
               QA=sapwide$QA,RM=sapwide$RM,RS=sapwide$RS,WP=sapwide$WP,YB=sapwide$YB,BS="0",GB="0",PC="0",
               RP="0",ST="0",WA="0",WS="0",cover="UN")

#now grab the overstory
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
overstory$EXP<-5
overstory<-filter(overstory, YEAR==2018)
overstory<-filter(overstory, SPP!="OH")
overstory<-filter(overstory, SPP!="UNK")
overstory<-filter(overstory, SPP!="OT")
unique(overstory$SPP)
overstory$BA<-(overstory$DBH^2)*0.005454
overstory<-overstory%>%
  mutate(BA_TPA = BA * EXP)
plot_summary <- overstory%>%
  group_by(SITEid, PLOTid, YEAR)%>%
  summarise(TPA_total = sum(EXP),
            BA_total = sum(BA_TPA))
overstory <- overstory%>%
  group_by(SITEid, PLOTid, YEAR, SPP)%>%
  summarise(TPA = sum(EXP),
            OVBA = sum(BA_TPA))
branch<-left_join(overstory, plot_summary)
branch <- branch%>%
  mutate(prop_tpa = (TPA/TPA_total),
         prop_ba = (OVBA/BA_total))
branch <- branch%>%
  mutate(iv = ((prop_tpa + prop_ba)/2))
branch$sapID<-paste0(branch$SITEid,"-",branch$PLOTid)
molten <- melt(as.data.frame(branch),id=c("sapID","iv","SPP"))
overwide <- dcast(molten,sapID~SPP,value.var = "iv",mean)
overwide[is.na(overwide)] <- 0
head(overwide)
names<-overwide$sapID
rownames(overwide)<-names
rowSums(overwide[2:17])
overwide$cover<-"overstory"
overwide$BC<-"0"
overwide$OT<-"0"
#bring together? i think that worked
allin<-rbind(overwide, new)

#seperate out into sp and env groups
sp<-allin[,c(1:17,19:20)]
env<-allin[,c(1,18)]
sp$BS<-as.numeric(sp$BS)
sp$GB<-as.numeric(sp$GB)
sp$PC<-as.numeric(sp$PC)
sp$RP<-as.numeric(sp$RP)
sp$ST<-as.numeric(sp$ST)
sp$WA<-as.numeric(sp$WA)
sp$WC<-as.numeric(sp$WC)
sp$WS<-as.numeric(sp$WS)
sp$BC<-as.numeric(sp$BC)
sp$OT<-as.numeric(sp$OT)
summary(sp)
#do the test
dis <- vegdist(sp[,2:19])
mod <- betadisper(dis,env$cover)
mod
test<-mrpp(dat=dis, grouping=env$cover, permutations=999)
test

centroids<-data.frame(grps=rownames(mod$centroids),data.frame(mod$centroids))
vectors<-data.frame(group=mod$group,data.frame(mod$vectors))
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")
grp1.hull<-seg.data[seg.data$group=="overstory",1:3][chull(seg.data[seg.data$group=="overstory",2:3]),]
grp2.hull<-seg.data[seg.data$group=="UN",1:3][chull(seg.data[seg.data$group=="UN",2:3]),]
all.hull<-rbind(grp1.hull,grp2.hull)

grp1.hull<-seg.data[seg.data$group=="overstory",1:3][chull(seg.data[seg.data$group=="old",2:4]),]
grp2.hull<-seg.data[seg.data$group=="UN",1:3][chull(seg.data[seg.data$group=="new",2:4]),]
all.hull<-rbind(grp1.hull,grp2.hull)

panel.a<-ggplot() +
  geom_polygon(data=all.hull[all.hull=="overstory",],aes(x=v.PCoA1,y=v.PCoA2),colour="blue",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data[seg.data$group %in% "overstory",],aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) + 
  geom_point(data=centroids[1,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=16) + 
  geom_point(data=seg.data[seg.data$group %in% "overstory",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=16,colour="blue") +
  labs(title="Overstory",x="",y="") +
  coord_cartesian(xlim = c(-0.7,0.7), ylim = c(-0.7,0.7)) +
  theme_classic() + 
  theme(legend.position="none")

panel.b<-ggplot() +
  geom_polygon(data=all.hull[all.hull=="UN",],aes(x=v.PCoA1,y=v.PCoA2),colour="orange",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data[seg.data$group %in% "UN",],aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) + 
  geom_point(data=centroids[2,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=17) + 
  geom_point(data=seg.data[seg.data$group %in% "UN",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=17,colour="orange") +
  labs(title="Understory",x="",y="") +
  coord_cartesian(xlim = c(-0.7,0.7), ylim = c(-0.7,0.7)) +
  theme_classic() + 
  theme(legend.position="none")
grid.arrange(panel.a,panel.b)

panel.e<-ggplot() + 
  geom_polygon(data=all.hull[all.hull=="overstory",],aes(x=v.PCoA1,y=v.PCoA2),colour="blue",alpha=0,linetype="dashed") +
  geom_polygon(data=all.hull[all.hull=="UN",],aes(x=v.PCoA1,y=v.PCoA2),colour="orange",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data,aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) + 
  geom_point(data=centroids[,1:3], aes(x=PCoA1,y=PCoA2,shape=grps),size=4,colour="red") + 
  geom_point(data=seg.data[seg.data$group %in% "overstory",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=16,colour="blue") +
  geom_point(data=seg.data[seg.data$group %in% "UN",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=17,colour="orange") +
  labs(x="",y="") +
  coord_cartesian(xlim = c(-0.7,0.7), ylim = c(-0.7,0.7)) +
  theme_classic() + 
  theme(legend.position="none")
panel.e

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


##Example graphics
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




#########Extra#######
# Updates using Mike's code from overstory ordination 
head(sapwide)
#sap.rda <- rda(sapwide ~ tmean+dew+actual.removed+wd.time, data=env, na.action = na.exclude)
sap.dist <- vegdist(sapwide)
sap.ord<-metaMDS(sapwide, distance = "bray", k=3,na.action=na.exclude)
sap.ord
sap.rdaall<-rda(sapwide~.,data=env,na.action="na.omit")
#env[is.na(env)] <- 0
efsap <- envfit(sapwide, env,na.rm=TRUE)
efsap
sapmod0 <- rda(sapwide~1,env,na.rm=TRUE)
sapmod1 <- rda(sapwide~.,env,na.rm=TRUE)
bstick(sapmod0)
screeplot(sapmod0,bstick=TRUE,type="lines")
summary(eigenvals(sapmod0))

set.seed(123)
mod <- ordistep(sapmod0, scope=formula(sapmod1),direction="both",na.action=na.exclude)
anova(mod,type="t")
anova(mod, by = "margin")
vif.cca(mod)
re.env <- env[c(11,27,36,24,29,4,26,38)]
ordisurf(vare.ord ~ actual.removed, env, bubble = 1, display="species")
text(vare.ord,display="spec",cex=1.5,col="blue")
fit <- ordisurf(vare.ord~wdi.time,env,family=quasipoisson)
calibrate(fit)


names<-sapwide$sapID
rownames(sapwide)<-names
final.over<-filter(final.over,bapa>0)
#join together sapling and environmental data to make sure the number of observations match
final.over<-na.omit(final.over)
bark<-left_join(sapwide, final.over)
#bark[is.na(bark)] <- 0
names<-bark$sapID #create list of sapID to set the rownames with 
#separate the datasets back out
sapwide<-bark[,c(1:12)]
env<-bark[,16:69]
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

