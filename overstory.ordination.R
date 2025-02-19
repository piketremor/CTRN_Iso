###################################################################
###################################################################
###################################################################
#######################OVERSTORY ORDINATION########################
###################################################################
###################################################################

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
library(MEForLab)
library(gridExtra)

setwd("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv")
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
set.seed(123)
over.ord<-metaMDS(overwide[,2:17], distance = "bray")
stressplot(over.ord)
over.ord
plot(over.ord)

plot(over.ord)
plot(over.ord,type="n")
points(over.ord,display="sites",cex=2,pch=21,col="red", bg="yellow")
text(over.ord,display="spec",cex=1.5,col="red")
over.ord
summary(over.ord)
over.ord$stress
#over.ord<-as.data.frame(scores(over.ord$species))
#species.scores <-
#speciesscores$species<-rownames(speciesscores)
#head(speciesscores)
#sitescores<-as.data.frame(scores(over.ord)$sites)
#sitescores$site<-rownames(sitescores)
#head(sitescores)

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

################RDA (add in environmental constraints)#########################
trees <- read.csv("Trees2023.csv")
locs <- read.csv("Tree_locations_species.csv")
tree_species <- locs[c(1:3,6)]
over <- left_join(trees, tree_species)
overstory <- over[c(2:9,27)]
overstory<-filter(overstory, YEAR == 2018)
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
final.over<-select(final.over, -c(SPP,Lithic,Northing_Y, Easting_X,bapa, tpa, qmd, RD, CCF, Shannon, over.Hill, ht40, prop.ws.avg,prop.bs.avg,prop.rs.avg,prop.hw.avg,prop.ab.avg,prop.eh.avg,prop.bf.avg,Redox,Parent,Min_depth,WD2000,WD2020))
final.over$THIN_METH<-as.factor(final.over$THIN_METH)
final.over$actual.removed[is.na(final.over$actual.removed)]<-0
names(final.over)
final.over$sapID<-paste0(final.over$SITEid,"-",final.over$PLOTid)
#change from integer to numeric
final.over$flowdir<-as.numeric(final.over$flowdir)
names(final.over)

#converting WD variables to WS variables
final.over$wdi.time<-final.over$wdi.time * -1
final.over$mean.WD<-final.over$mean.WD * -1
final.over$wd.time<-final.over$wd.time * -1

#renaming the WD variables accordingly 

colnames(final.over)[colnames(final.over) == 'wdi.time'] <- 'wsi.time'
colnames(final.over)[colnames(final.over) == 'wd.time'] <- 'ws.time'
colnames(final.over)[colnames(final.over) == 'mean.WD'] <- 'mean.WS'


#final.over<-na.omit(final.over)
#join together sapling and environmental data to make sure the number of observations match
bark<-left_join(final.over, overwide)
names<-bark$sapID #create list of sapID to set the rownames with 
#separate the datasets back out
overwide<-bark[,c(43:59)]
names(overwide)
env<-bark[,4:43]
names(env)
#add rownames
overwide<-as.data.frame(overwide)
env<-as.data.frame(env)
rownames(overwide)<-names
rownames(env)<-names
#remove sapID variable
overwide<-overwide[,2:17]
names(overwide)
names(env)
head(env)
head(overwide)
env<-select(env, -sapID)
names(env)

#write.csv(final.over,"CTRN_EnvMatrix.csv")


# Premer pickup here. 
head(overwide)
set.seed(123)
#sap.rda <- rda(sapwide ~ tmean+dew+actual.removed+wd.time, data=env, na.action = na.exclude)
vare.dist <- vegdist(overwide)
vare.ord<-metaMDS(overwide, distance = "bray",
                  trymax=250,autotransform =TRUE,k=3)
over.rdaall<-rda(overwide~.,data=env,na.action="na.omit")
env[is.na(env)] <- 0
ef <- envfit(overwide, env,na.action="na.pass")
ef
plot(ef)
mod0 <- rda(overwide~1,env)
mod1 <- rda(overwide~.,env)


bstick(mod0)
screeplot(mod0,bstick=TRUE,type="lines")
summary(eigenvals(mod0))


#sig<-select(env, c(WHC,dew,roughness))
#siggy<-envfit(overwide,sig,na.action="na.pass")
#plot(siggy)

set.seed(123)
mod <- ordistep(mod0, scope=formula(mod1),direction="both")
summary(mod)
#PCT,dew,ph,tst,WHC,ex.k,roughness,ex.ca,THIN_METH, actual.removed
anova(mod,type="t")
anova(mod, by = "margin")

#significance of constraints
anova(mod, by="term")
dis<-vegdist(overwide)
test1<-adonis2(mod,data=env,by="terms")
vif.cca(mod)
re.env <- select(env, c(dew,ph,tst,WHC,ex.k,roughness,ex.ca,THIN_METH,actual.removed))
ordisurf(vare.ord ~ actual.removed, re.env, bubble = 1)
text(vare.ord,display="spec",cex=1,col="blue")
env$wdi.time<-env$wdi.time*-1
fit <- ordisurf(vare.ord~WHC,re.env,family=quasipoisson)
calibrate(fit)

RsquareAdj(mod)$adj.r.squared
summary(mod)
final.mod<-rda(overwide~.,re.env)
anova(final.mod, by="term")
#species correlations with axes and associated p values
ce <- cor(overwide,method="pearson",scores(vare.ord,dis="si"))
try.x <- (ce*sqrt(96-2)/
            sqrt(1-ce^2))
try.x
dt(try.x,df=95)


#cannot run categorical variables with correlations
cx <- cor(re.env,method="pearson",scores(vare.ord,dis="si"))
try <- (cx*sqrt(96-2)/
          sqrt(1-cx^2))
try
try.p<-dt(try,df=95)


#write.csv(cx, "~/Desktop/overstory.environmental.csv")
#write.csv(try.p,"~/Desktop/overstory.pval.csv")

## Premer end. 12/6/2024

#head(sapwide)
#sap.rda <- rda(sapwide ~ tmean+dew+actual.removed+wd.time, data=env, na.action = na.exclude)
#sap.rdaall<-rda(overwide~.,data=env,na.action="na.omit")

#summary(sap.rdaall)
#summary(sap.rdaall)
#ordiplot(sap.rdaall, scaling = 2, type = "text")
#step.forward <- ordiR2step(rda(overwide ~ 1, data=env), scope=formula(sap.rdaall), R2scope = F, direction="forward", pstep=1000)
#anova(sap.rdaall, by="terms", step=1000) 
#anova(sap.rdaall, step=1000)


#names(env)


env<-select(env, c(PCT, dew, ph,tst,WHC,ex.k, THIN_METH))
sapordi<-metaMDS(overwide, distance = "bray")
en<-envfit(overwide, env, na.action="na.pass")
summary(sapordi)
sapordi
plot(sapordi)
stressplot(sapordi)
anova(en,by="term")

#plotting with the ellipses grouped by thinning method
gof<-goodness(sapordi)
plot(en)
trt<-env$THIN_METH
data.scores<-as.data.frame(scores(mod)$species) 
data.scores$site <-as.data.frame(scores(mod)$sites)
data.scores$THIN_METH = trt
ggplot(data=data.scores) + 
  stat_ellipse(aes(x=RDA1,y=RDA2,colour=THIN_METH),level = 0.50) +
  geom_point(aes(x=RDA1,y=RDA2,colour=THIN_METH),size=4) + 
  theme_classic()+
  labs(color="Thinning Method")
  
adon.results<-adonis2(sapwide ~ trt, method="bray",perm=999)
print(adon.results)


slice <- left_join(over.ID,slice,by="uid")
slice$tst <- as.numeric(slice$tst)
slice$age <- as.numeric(slice$age)
slice.ff <- slice[order(slice$uid),]
ff.all <- slice.ff[,c(2:17)]
ff.env <- slice.ff[,c(2:15,18,19)]













#mulitvariate dispersion grouped by thinning method
dis <- vegdist(overwide)
mod <- betadisper(dis, trt)
mod

test<-mrpp(dat=dis, grouping=trt, permutations=999)

#adonis (PERMANOVA) test

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
  coord_cartesian(xlim = c(-0.5,0.7), ylim = c(-0.25,0.2)) +
  theme_classic() + 
  theme(legend.position="none")

panel.b<-ggplot() +
  geom_polygon(data=all.hull[all.hull=="crown",],aes(x=v.PCoA1,y=v.PCoA2),colour="orange",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data[seg.data$group %in% "crown",],aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) + 
  geom_point(data=centroids[2,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=17) + 
  geom_point(data=seg.data[seg.data$group %in% "crown",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=17,colour="orange") +
  labs(title="Crown",x="",y="") +
  coord_cartesian(xlim = c(-0.5,0.7), ylim = c(-0.25,0.2)) +
  theme_classic() + 
  theme(legend.position="none")

panel.c<-ggplot() +   
  geom_polygon(data=all.hull[all.hull=="dominant",],aes(x=v.PCoA1,y=v.PCoA2),colour="green",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data[seg.data$group %in% "dominant",],aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) +
  geom_point(data=centroids[3,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=15) + 
  geom_point(data=seg.data[seg.data$group %in% "dominant",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=15,colour="green") + 
  labs(title="Dominant",x="",y="") +
  coord_cartesian(xlim = c(-0.5,0.7), ylim = c(-0.5,0.4)) +
  theme_classic() + 
  theme(legend.position="none")


panel.d<-ggplot() +
  geom_polygon(data=all.hull[all.hull=="low",],aes(x=v.PCoA1,y=v.PCoA2),colour="purple",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data[seg.data$group %in% "low",],aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) +
  geom_point(data=centroids[3,1:3], aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=15) + 
  geom_point(data=seg.data[seg.data$group %in% "low",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=15,colour="purple") +
  labs(title="Low",x="",y="") +
  coord_cartesian(xlim = c(-0.5,0.7), ylim = c(-0.25,0.2)) +
  theme_classic() + 
  theme(legend.position="none")


panel.e<-ggplot() + 
  geom_polygon(data=all.hull[all.hull=="low",],aes(x=v.PCoA1,y=v.PCoA2),colour="purple",alpha=0,linetype="dashed") +
  geom_polygon(data=all.hull[all.hull=="dominant",],aes(x=v.PCoA1,y=v.PCoA2),colour="green",alpha=0,linetype="dashed") +
  geom_polygon(data=all.hull[all.hull=="crown",],aes(x=v.PCoA1,y=v.PCoA2),colour="orange",alpha=0,linetype="dashed") +
  geom_polygon(data=all.hull[all.hull=="control",],aes(x=v.PCoA1,y=v.PCoA2),colour="blue",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data,aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) + 
  geom_point(data=centroids[,1:3], aes(x=PCoA1,y=PCoA2,shape=grps),size=4,colour="red") + 
  geom_point(data=seg.data[seg.data$group %in% "low",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=15,colour="purple") +
  geom_point(data=seg.data[seg.data$group %in% "dominant",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=15,colour="green") +
  geom_point(data=seg.data[seg.data$group %in% "crown",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=17,colour="orange") +
  geom_point(data=seg.data[seg.data$group %in% "control",], aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=16,colour="blue") +
  labs(title="All",x="",y="") +
  coord_cartesian(xlim = c(-0.5,0.7), ylim = c(-0.5,0.4)) +
  theme_classic() + 
  theme(legend.position="none")

grid.arrange(panel.a,panel.b,panel.c,panel.d, panel.e, nrow=2)
#more graphing
data.scores = as.data.frame(scores(mod)$sites)
data.scores$THIN_METH = env$THIN_METH

Type<-c("Softwood", "Softwood","Softwood")

en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)
speciesscores<-as.data.frame(scores(mod, display = "species"))

gg = ggplot(data = data.scores, aes(x = RDA1, y = RDA2)) + 
  geom_point(data = data.scores, aes(colour = THIN_METH), size = 3, alpha = 0.5) + 
  scale_colour_manual(values = c("orange", "steelblue", "lightblue2", "gold1"))  + 
  geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  geom_point(data = en_coord_cat, aes(x = RDA1, y = RDA2), 
             shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  geom_text(data = en_coord_cat, aes(x = RDA1, y = RDA2+0.04), 
            label = row.names(en_coord_cat), colour = "navy", fontface = "bold") + 
  geom_text(data = en_coord_cont, aes(x = RDA1, y = RDA2), colour = "grey30", 
            fontface = "bold", label = row.names(en_coord_cont)) + 
  geom_text(data=speciesscores,aes(x=RDA1,y=RDA2, label=row.names(speciesscores)),size=5,vjust=0, colour="navy")+
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Thinning Method")

gg


ordiplot(final.mod, scaling = 2, type = "none")%>%
  points("sites", pch=16)%>%
  text("species", col="red")%>%
  text(what="biplot", col="blue")


#graphing grouped by softwood or hardwood

grp<-c("Softwood", "Softwood","Softwood", "Hardwood", "Hardwood", "Hardwood","Hardwood","Hardwood","Softwood","Softwood","Hardwood","Hardwood","Softwood","Softwood","Softwood","Hardwood")
speciesscores$grp<-grp
speciesscores<-na.omit(speciesscores)

grp.a <- speciesscores[speciesscores$grp == "Softwood", ][chull(speciesscores[speciesscores$grp == 
                                                                   "Softwood", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.b <- speciesscores[speciesscores$grp == "Hardwood", ][chull(speciesscores[speciesscores$grp == 
                                                                   "Hardwood", c("NMDS1", "NMDS2")]), ]  # hull values for grp B

hull.data <- rbind(grp.a, grp.b)  #combine grp.a and grp.b
hull.data

ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=grp,group=grp),alpha=0.30) + # add the convex hulls
  geom_text(data=speciesscores,aes(x=NMDS1,y=NMDS2,label=rownames(speciesscores)),alpha=1.5) +  # add the species labels
  #geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=4) + # add the point markers
  #scale_colour_manual(values=c("Softwood" = "blue", "Hardwood" = "red")) +
  coord_equal() +
  theme_bw() + 
  labs(fill="Type")+ 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
