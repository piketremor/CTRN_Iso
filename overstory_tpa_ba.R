#edits made by Lila Beck on 11/3/23 10:10am
#load packages
library(dplyr)
library(mosaic)
library(forcats)
#library(tidyverse)

#load in data
rm(list=ls())
setwd("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv")
setwd("~/Google Drive/My Drive/Research/CFRU/CTRN_CFRU_Share/raw/csv")

#setwd("G:/.shortcut-targets-by-id/1sCbm2t1PUIpbJYVOzlIrZKeVhisl4Ghv/CTRN_CFRU_Share/raw/csv")
trees <- read.csv("Trees2023.csv")
tree_locations <- read.csv("Tree_locations_species.csv")
#trees <- read.csv("C:/users/lila.beck/Desktop/CTRN-Thesis data/raw/Trees.csv")
#tree_locations <- read.csv("C:/Users/lila.beck/Desktop/CTRN-Thesis data/raw/Tree_locations_species.csv")

#######################Basal Area and TPA#######################################
#select for SITEid, PLOTid, and TREE to be in tree locations dataset
tree_species <- tree_locations%>%
  select(SITEid, PLOTid, TREE, SPP)

#join the datasets
tree_join <- left_join(trees, tree_species)
tree_join <- filter(tree_join, SITEid == "AS" | SITEid == "DR" | SITEid == "GR" | SITEid == "HR" | SITEid == "KI" | SITEid == "LM" | SITEid == "LT" | SITEid == "PA" | SITEid == "PE" | SITEid == "RC" | SITEid == "RR" | SITEid == "SA" | SITEid == "SC" | SITEid == "SR" | SITEid == "WB")


#filter for the living trees
tree_join_alive <- filter(tree_join, STATUS == 1 | STATUS == 2 | STATUS == 4 | STATUS == 10 | STATUS == 5 | STATUS == 14)

#add in expansion factor
tree_join_alive <- tree_join_alive%>%
  mutate(TREES_EXP = 5,
         BA = (DBH^2)*0.005454,
         TREE_BAPA = BA * TREES_EXP)


#remove NA values in basal area


tree_join_alive$BA[is.na(tree_join_alive$BA)]<-0

#basal area and tpa calculations for plot level per year
plot_summary <- tree_join_alive%>%
  group_by(SITEid, PLOTid, YEAR)%>%
  summarise(TPA_TOTAL = sum(TREES_EXP),
            BAPA = sum(TREE_BAPA))

#basal area and tpa calculations by species
species_summary <- tree_join_alive%>%
  group_by(SITEid, PLOTid, YEAR, SPP)%>%
  summarise(TPA = sum(TREES_EXP),
            TreeBA = sum(TREE_BAPA))

#join species and plot tables together
trees_again <- left_join(species_summary, plot_summary)

#############QMD############################################################
trees_again<-trees_again%>%
  mutate(QMD = sqrt((BAPA/TPA_TOTAL)/0.005454))

ov.metrics<-trees_again%>%
  group_by(SITEid, PLOTid, YEAR)%>%
  summarise(TPA_TOTAL = mean(TPA_TOTAL),
            BAPA = mean(BAPA),
            QMD = mean(QMD))


#full<-left_join(ov.metrics, pp.ov)
##########################Percent Spruce/Fir/HW################################
percent<-trees_again%>%
  mutate(percent.spp = (TPA/TPA_TOTAL)*100)

trees_again.rs<-percent%>%
  filter(SPP=="RS")%>%
  mutate(pp.rs = percent.spp)

trees_again.rs<-trees_again.rs[,c(1:4,11)]

trees_again.ws<-percent%>%
  filter(SPP=="WS")%>%
  mutate(pp.ws = percent.spp)

trees_again.ws<-trees_again.ws[,c(1:4,11)]

trees_again.bs<-percent%>%
  filter(SPP=="BS")%>%
  mutate(pp.bs = percent.spp)

trees_again.bs<-trees_again.bs[,c(1:4,11)]

trees_again.bf<-percent%>%
  filter(SPP=="BF")%>%
  mutate(pp.bf = percent.spp)

trees_again.bf<-trees_again.bf[,c(1:4,11)]

trees_again.hw<-percent%>%
  filter(SPP=="RM" | SPP=="WA"|SPP=="YP"|SPP=="SM"|SPP=="MM"|SPP=="QA"|SPP=="BC"|SPP=="AB"|SPP=="GB")%>%
  mutate(pp.hw = percent.spp)

trees_again.hw<-trees_again.hw[,c(1:4,11)]

trees_again<-left_join(trees_again,trees_again.rs, by=join_by(SITEid,PLOTid,YEAR,SPP))

trees_again<-left_join(trees_again,trees_again.ws, by=join_by(SITEid,PLOTid,YEAR,SPP))

trees_again<-left_join(trees_again,trees_again.bs, by=join_by(SITEid,PLOTid,YEAR,SPP))

trees_again<-left_join(trees_again, trees_again.bf, by=join_by(SITEid,PLOTid,YEAR,SPP))

trees_again<-left_join(trees_again, trees_again.hw, by=join_by(SITEid,PLOTid,YEAR,SPP))

#trees_again<-select(trees_again, -c("TreeBA.x","TPA_TOTAL.x","BAPA.x", "QMD.x", "TPA.x","TreeBA.y","TPA_TOTAL.y","BAPA.y","QMD.y","TPA.y","percent.spp.x","percent.spp.y"))


trees_again$pp.rs[is.na(trees_again$pp.rs)]<-0
trees_again$pp.ws[is.na(trees_again$pp.ws)]<-0
trees_again$pp.bs[is.na(trees_again$pp.bs)]<-0
trees_again$pp.bf[is.na(trees_again$pp.bf)]<-0
trees_again$pp.hw[is.na(trees_again$pp.hw)]<-0

pp.summary<-trees_again%>%
  group_by(SITEid, PLOTid, YEAR)%>%
  summarise(pp.rs.avg = mean(pp.rs),
            pp.ws.avg = mean(pp.ws),
            pp.bs.avg = mean(pp.bs),
            pp.bf.avg = mean(pp.bf),
            pp.hw.avg = mean(pp.hw))

#write.csv(pp.summary, "~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv/percent_summary.csv")


#trees_again<-select(trees_again, "SITEid", "PLOTid", "YEAR", "SPP", "pp.spruce", "pp.bf", "pp.hw")

#metrics<-read.csv("overstory_metrics.csv")

#full<-left_join(metrics, pp.summary)



#write.csv(full, "~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv/overstory_metrics.csv")


##############Relative Density###################################################
ov<-read.csv("overstory_metrics.csv")

ov<-ov%>%
  mutate(
    relative.density = BAPA/sqrt(QMD)
  )


#write.csv(ov, "~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv/overstory_metrics.csv")
###############CCF#############################################################
library(MEForLab)

#calculating maximum crown width
tree_join_alive["MCW"]<-
  mapply(MCW, SPP = tree_join_alive$SPP, DBH = tree_join_alive$DBH)

#calculating crown area
tree_join_alive["crown.area"]<-pi*((tree_join_alive$MCW/2)^2)

#multiply by expansion factor
tree_join_alive<-tree_join_alive%>%
  mutate(crown.area = crown.area*5)

#calculating crown competition factor
tree_join_alive<-tree_join_alive%>%
  group_by(SITEid, PLOTid, YEAR)%>%
  summarise(plot.ca = sum(crown.area))

tree_join_alive["CCF"]<-tree_join_alive$plot.ca/43560
#tree_join_alive<-tree_join_alive[,c(1:3,5)]

#ov<-read.csv("overstory_metrics.csv")

check <- left_join(trees_again,tree_join_alive)

plot(check$BAPA,check$CCF)

# looks good here !! 

#trees<-left_join(ov,tree_join_alive)

#write.csv(trees, "~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv/overstory_metrics.csv")

################################################################################
#importance values
trees_again <- trees_again%>%
  mutate(prop_tpa = (TPA/TPA_TOTAL),
         prop_ba = (TreeBA/BAPA))

trees_again <- trees_again%>%
  mutate(iv = ((prop_tpa + prop_ba)/2))

library(lattice)
xyplot(BAPA~YEAR|SITEid,data=plot_summary,type="p")

spruce <- plot_summary%>%
  group_by(SITEid,PLOTid)%>%
  arrange(YEAR,.by_group = TRUE)

picea <- spruce %>%
  group_by(SITEid,PLOTid)%>%
  mutate(delta.ba = BAPA-lag(BAPA))
picea[is.na(picea)] <- 0
picea$ba.growth <- ifelse(picea$delta.ba>0,picea$delta.ba,picea$delta.ba)

xyplot(ba.growth~YEAR|SITEid,data=picea)

plots <- read.csv("Plots.csv")
ctrn <- left_join(picea,plots)
ctrn <- ctrn%>%
  #filter(.,THIN_METH!="control")%>%
  mutate(TST = YEAR-TRT_YR)%>%
  filter(.,TST>0)

ctrn2 <- ctrn%>%
  group_by(SITEid,PLOTid)%>%
  mutate(comp.ba = cumsum(ba.growth),
         mai = comp.ba/TST,
         prev.ba = lag(BAPA))

sites <- read.csv("CTRN_Site_Variables.csv")
ctrn3 <- left_join(ctrn2,sites)

xyplot(comp.ba~TST|THIN_METH,data=ctrn3,)
ctrn3$wdi <- ((ctrn3$WD2000+ctrn3$WD2020)/2)-ctrn3$SWC2
ctrn3$wd <- (ctrn3$WD2000+ctrn3$WD2020)/2
ctrn3 <- dplyr::filter(ctrn3,SITEid!="AP")
#ctrn3$mai[ctrn3$mai]

#latest <- ctrn3%>%
#  group_by(SITEid)%>%
#  top_n(n=1,wt=YEAR)
ctrn3$REMOVAL[is.na(ctrn3$REMOVAL)] <- 0
ctrn3$GE <- ctrn3$ba.growth/ctrn3$prev.ba
ctrn3 <- dplyr::filter(ctrn3,)

trt.sum <- ctrn3%>%
  group_by(SITEid,PLOTid)%>%
  summarize(TRT = mean(TRT_YR))

#ctrn3 <- dplyr::filter(ctrn3,mai<8)


easy <- lm(BAPA~+TST+THIN_METH+wd+REMOVAL+PCT
           ,data=ctrn3)
summary(easy)

ctrn3$comp.ba[is.na(ctrn3$comp.ba)] <- 0
#plot(density(log(ctrn3$comp.ba)))

#write.csv(ctrn3,"~/Desktop/CTRN_Enviro.csv")

#ctrn3 <- dplyr::filter(ctrn3,SITEid!="AP")
library(nlme)
#ctrn3$REMOVAL <- as.factor(ctrn3$REMOVAL)
#ctrn3$ba.growth[ctrn3$ba.growth>25] <- 25
#ctrn3$qmd <- sqrt(((ctrn3$BAPA/ctrn3$TPA)/0.005454))
ctrn3$GE[is.na(ctrn3$GE)] <- 999
#ctrn4 <- dplyr::filter(ctrn3,mai<10)
mdl.ac <- gls(BAPA~wd+(REMOVAL:TST)+PCT+REMOVAL+vpdmax, data=ctrn3, 
              correlation = corAR1(form=~1|YEAR/PLOTid/SITEid),
              na.action=na.omit)

summary(mdl.ac)
plot(mdl.ac)
AIC(mdl.ac)

car::vif(easy)
library(ggeffects)

AIC(easy,mdl.ac)


mydf2 <- ggpredict(mdl.ac,terms=c("TST","REMOVAL","wd"))


#png("CTRN_BA_Prelim3.png",units='in',height=6,width=15,res=1000)
theme_set(theme_bw(16))

library(ggplot2)
theme_set(theme_bw(18))
ggplot(mydf2,aes(x=x,y=predicted,colour=group))+
  geom_line(aes(linetype=group,color=group),size=1)+
  labs(x="Time since Treatment (Years)",y="BAPA post-treatment")+
  labs(linetype="RD Reduction (%)")+
  labs(colour = "RD Reduction (%)")+
  #xlim(5,20)+
  #ylim(0,25)+ƒ
  facet_wrap(~facet)+
  theme_bw(18) +
  #theme(legend.position="none")
  scale_color_manual(values=c('gray0','gray70','gray40'))+
  scale_fill_manual(values=c('gray0','gray70','gray40'), name="fill")

dev.off()

# i think we have enough to work with. 

goose <- picea%>%
  group_by(SITEid,PLOTid)


#QMD
trees_again<-trees_again%>%
  mutate(QMD = sqrt((BAPA/TPA_TOTAL)/0.005454))

# nice, you've got the right idea, but try it with a function

qmd <- function(bapa,tpa){
  qmd = sqrt((bapa/tpa)/0.005454)
  return(qmd)
}

trees_again$qmd <- qmd(trees_again$BAPA,trees_again$TPA_TOTAL)
honk<-filter(trees_again, SITEid == "SR" & YEAR == 2018)

# Premer updates 10/17/2023
ht.fit <- tree_join_alive%>%
  filter(.,SITEid=="SR"&TOT_HT>0&YEAR==2018)
mod1 <- lm(TOT_HT~log(DBH+0.01),data=ht.fit)


#everything in functions
sp.fit <- function(dbh){
  fit.ht = (mod1$coefficients[1]+mod1$coefficients[2]*log(dbh+0.01))
  return(fit.ht)
}

tree_join_alive$TOT_HT[is.na(tree_join_alive$TOT_HT)] <- 0
sr.df <- tree_join_alive%>%
  left_join(plots)%>%
  filter(.,SITEid=="SR"&YEAR==2018)%>%
  mutate(fit.ht = sp.fit(DBH+0.01))%>% # use the function here. 
  group_by(PLOTid,REMOVAL,THIN_METH,TRT_YR)%>%
  summarise(bapa = sum(TREE_BAPA),
            tpa = sum(TREES_EXP),
            qmd = sqrt((bapa/tpa)/0.005454))

ht.pop <- tree_join_alive %>%
  left_join(plots)%>%
  filter(.,SITEid=="SR"&YEAR==2018)%>%
  mutate(fit.ht = sp.fit(DBH))%>%
  group_by(PLOTid)%>%
  slice_max(.,DBH,n=8)

ducks <- ht.pop%>%
  group_by(PLOTid)%>%
  summarize(ht.40 = mean(fit.ht))%>%
  left_join(sr.df,.)
 
xyplot(TOT_HT~DBH|SPP,data=ht.fit)

sr <- tree_join_alive%>%
  filter(.,SITEid=="SR"&YEAR==2018)%>%
  mutate(fit.ht = sp.fit(DBH+0.01))

xyplot(fit.ht~DBH|SPP,data=sr,type="l")

rd <- function(ba,qmd){
  curtis = ba/sqrt(qmd)
  return(curtis)
} 

ducks$rd <- rd(ducks$bapa,ducks$qmd)
print(ducks)


write.csv(ducks,"SarahsRoad.csv")



#quack<-honk%>%
#  left_join(.,plots)%>%
#  group_by(PLOTid,REMOVAL,PCT,THIN_METH,TRT_YR)%>%
#  summarise(QMD = mean(QMD),
#            BAPA = mean(TreeBA),
#            TPA = mean(TPA))


#overstory split into pre-treat, post-treat, 10 yrs post-treat

plots <- read.csv("Plots.csv")

treat <- left_join(trees_again, plots)

treat<-treat%>%
  mutate(TST = YEAR - TRT_YR)


pre<-filter(treat, TST == 1)

post<-filter(treat, TST == -1)

ten <-filter(treat, TST == 10)

#i think this works?