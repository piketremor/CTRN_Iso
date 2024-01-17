# edits by Lila Beck on 10/7/23
dev.off
rm(list=ls())

#load packages
library(dplyr)
library(mosaic)
library(forcats)
library(tidyverse)
library(openxlsx)
library(forestmangr)
library(vegan)

#saplings <- read.csv("C:/users/lila.beck/Desktop/CTRN-Thesis data/raw/Saplings.csv")
saplings <- read.csv("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv/Saplings.csv")
plot <- read.csv("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv/Plots.csv")
#plot <- read.csv("C:/Users/lila.beck/Desktop/CTRN-Thesis data/raw/Plots.csv")

#saplings<-replace(saplings, is.na(saplings), 0)

saplings <- saplings%>%
  mutate(ba.half = (0.5^2*0.005454)*X1.2.inch,
         ba.one = (1.0^2*0.005454)*X1.inch,
         ba.two = (2.0^2*0.005454)*X2.inch,
         SAP_EXP = 250)
         
  

saplings<-saplings%>%
  mutate(bapa = (ba.half+ba.one+ba.two)*SAP_EXP)


plot_summary_sap <- saplings%>%
  group_by(SITEid, PLOTid, YEAR, CORNERid)%>%
  summarise(total.tpa = sum(SAP_EXP),
            total.bapa = sum(bapa))

species_summary_sap <- saplings%>%
  group_by(SITEid, PLOTid, YEAR, CORNERid, SPP)%>%
  summarise(spp.tpa = sum(SAP_EXP),
            spp.bapa = sum(bapa))

saplings_again<-left_join(species_summary_sap, plot_summary_sap)

saplings_again <- saplings_again%>%
  mutate(prop_tpa = (spp.tpa/total.tpa),
         prop_ba = (spp.bapa/total.bapa))

saplings_again <- saplings_again%>%
  mutate(iv = ((prop_tpa + prop_ba)/2))

#join the saplings_again data with plot data
plot_treatments <- plot%>%
  select(SITEid, PLOTid, REMOVAL, THIN_METH, TRT_YR)
saplings_treatment <- left_join(plot_treatments, saplings_again)

#saplings_treatment<-replace(saplings_treatment, is.na(saplings_treatment), 0)

saplings_treatment<-saplings_treatment%>%
  mutate(TST = YEAR - TRT_YR) #TST = time since thinning

target <- dplyr::filter(saplings_treatment,SPP=="RS"|SPP=="RM"|SPP=="BF")
library(lattice)



xyplot(iv~TST|SPP*THIN_METH,data=target,type="p")
target[is.na(target)] <- 0

sites <- read.csv("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv/CTRN_Site_Variables.csv")
#sites <- read.csv("C:/Users/lila.beck/Desktop/CTRN-Thesis data/raw/CTRN_Site_Variables.csv")
summary(sites)
cal <- left_join(target,sites)
rs <- dplyr::filter(cal,SPP=="RS")

lm1 <- lm(iv~THIN_METH+YEAR,data=rs)
summary(lm1)
# suggests there is a larger effect across sites than treatments

# look at them fall off at ~2012 - i wonder if there was a big recruitment event. 

# now, bring in the overstory data and join it so that you have overstory ba and species proportion, 
trees <- read.csv("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv/Trees2023.csv")
#trees <- read.csv("C:/Users/lila.beck/Desktop/CTRN-Thesis data/raw/Trees2023.csv")
tree_locations <- read.csv("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv/Tree_locations_species.csv")
#tree_locations <- read.csv("C:/Users/lila.beck/Desktop/CTRN-Thesis data/raw/Tree_locations_species.csv")
summary(trees)


tree_species <- tree_locations%>%
  select(SITEid, PLOTid, TREE, SPP)

colnames(tree_species) <- c("SITEid", "PLOTid", "TREE", "TREE_SPP")

tree_join <- left_join(trees, tree_species)


tree_join_alive <- filter(tree_join, STATUS == 1 | STATUS == 2 | STATUS == 4 | STATUS == 10 | STATUS == 5 | STATUS == 14)

tree_join_alive <- tree_join_alive%>%
  mutate(TREES_EXP = 5,
         BA = (DBH^2)*0.005454,
         TREE_BAPA = BA * TREES_EXP)

plot_summary <- tree_join_alive%>%
  group_by(SITEid, PLOTid, YEAR)%>%
  summarise(TPA = sum(TREES_EXP),
            BAPA = sum(TREE_BAPA))

species_summary <- tree_join_alive%>%
  group_by(SITEid, PLOTid, YEAR, TREE_SPP)%>%
  summarise(spp_TPA = sum(TREES_EXP),
            TreeBA = sum(TREE_BAPA))


trees_again <- left_join(species_summary, plot_summary)


trees_again <- trees_again%>%
  mutate(treeprop_tpa = (spp_TPA/TPA),
         treeprop_ba = (TreeBA/BAPA))

trees_again <- trees_again%>%
  mutate(tree_iv = ((treeprop_tpa + treeprop_ba)/2))
summary(trees_again)

rs2 <- filter(trees_again, TREE_SPP == "RS")


cal2 <- left_join(rs, rs2)



# then, re-run lm1 with the addition of the RS over story basal area before and after thinning, and calculate the change in ba
lm2 <- lm(iv~THIN_METH+YEAR+total.bapa,data=rs)
summary(lm2) 


