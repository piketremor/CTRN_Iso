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


setwd("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv")
saplings <- read.csv("Saplings.csv")[,2:9]
sap_tpa_ba <- read.csv("sapling_tpa_ba.csv")[,2:12]

saplings <- filter(saplings, SITEid == "AS" | SITEid == "DR" | SITEid == "GR" | SITEid == "HR" | SITEid == "KI" | SITEid == "LM" | SITEid == "LT" | SITEid == "PA" | SITEid == "PE" | SITEid == "RC" | SITEid == "RR" | SITEid == "SA" | SITEid == "SC" | SITEid == "SR" | SITEid == "WB") 
saplings[saplings == "SpecAld"]<-"SA"
saplings[saplings == "HM"]<-"EH"
saplings[saplings == "CH"]<-"BC"


#diversity 
saplings<-saplings%>%
  mutate(X1.2.inch = replace_na(X1.2.inch, 0),
         X1.inch = replace_na(X1.inch, 0),
         X2.inch = replace_na(X2.inch, 0))

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


tree_species <- locs[c(1:3,6)]
tree_species <- locs%>%
  select(SITEid, PLOTid, TREE, SPP)

over <- left_join(trees, tree_species)
over$tally <- 1

over.sum<-over%>%
  group_by(SITEid, PLOTid, YEAR, SPP)%>%
  summarize(over.total = sum(tally))

all.over<-over%>%
  group_by(SITEid, PLOTid, YEAR)%>%
  summarize(overstory.total = sum(tally))

overstory<- left_join(over.sum, all.over)
head(overstory)

all_tree<-left_join(overstory, sappy, by = c("SITEid", "PLOTid", "SPP", "YEAR"))
all_tree$sap.total[is.na(all_tree$sap.total)] <- 0
all_tree$sapling.total[is.na(all_tree$sapling.total)] <- 0
all_tree$SPP[is.na(all_tree$SPP)]<-"OTHER CONIFER"

all_tree$prop <- (all_tree$over.total + all_tree$sap.total)/(all_tree$overstory.total+all_tree$sapling.total)
head(all_tree)

all_tree$shann.base <- all_tree$prop*(log(all_tree$prop))
all_tree$shann.base <- -1*(all_tree$shann.base)

shann.frame <- all_tree%>%
  group_by(SITEid, PLOTid, YEAR)%>%
  mutate(shannon = sum(shann.base),
         hill = exp(shannon))

shann.sum <- shann.frame%>%
  group_by(SITEid, PLOTid, YEAR)%>%
  summarise(shannon = mean(shannon))

View(shann.frame)




