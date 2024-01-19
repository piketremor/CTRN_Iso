#edits made by Lila Beck on Dec. 21, 2023
#load packages
library(dplyr)
library(tidyverse)
library(lme4)
library(MEForLab)

#load data
setwd("~/Google Drive/My Drive/Research/CFRU/CTRN_CFRU_Share/raw/csv")
#setwd("G:/My Drive/Research/CFRU/CTRN_CFRU_Share/raw/csv")
overstory<-read.csv("Trees2023.csv")
tree_locations <- read.csv("tree_locations_species.csv")


#select for SITEid, PLOTid, and TREE to be in tree locations dataset
tree_species <- tree_locations%>%
  select(SITEid, PLOTid, TREE, SPP)

#join the datasets
tree_join <- left_join(overstory, tree_species)

#filter for the living trees
tree_join <- filter(tree_join, STATUS == 1 | STATUS == 2 | STATUS == 4 | STATUS == 10 | STATUS == 5 | STATUS == 14)

#filter for heights above 0
height.frame <- filter(tree_join, TOT_HT > 0)

#create model
ht.lm <- lmer(TOT_HT ~ log(DBH) + (1 | SPP/YEAR/SITEid), data = height.frame)
summary(ht.lm)

#predict heights
tree_predict <- tree_join %>% 
  mutate(PRD_HT = predict(ht.lm, tree_join, re.form = NULL, allow.new.levels = TRUE))

tree_predict$TOT_HT[is.na(tree_predict$TOT_HT)] <- 0
tree_predict$PRD_HT[is.na(tree_predict$PRD_HT)] <- 0
tree_predict$fin.ht <- ifelse(tree_predict$TOT_HT<1,tree_predict$PRD_HT,tree_predict$TOT_HT)
tree_predict$fin.ht <- ifelse(tree_predict$fin.ht<0,0,tree_predict$fin.ht)

tree_predict$SPP[is.na(tree_predict$SPP)] <- "OT"



#here's where i'm getting stuck...
#apply the function to the dataframe and create new volume variable


tree_predict["vol"] <- 
  mapply(vol_calc,SPP=tree_predict$SPP,DBH=tree_predict$DBH,HT=tree_predict$fin.ht)


library(lattice)

xyplot(vol~DBH|SPP,data=tree_predict)

xyplot(fin.ht~DBH|SPP,data=tree_predict)

xplot(fin.ht~DBH,data=tree_predict)



##scratch##  
#gdf <- groupedData(TOT_HT ~ log(DBH)|SPP/YEAR/SITEid, data=tree_join)

#plotlm = lmList(TOT_HT~log(DBH), data=gdf, na.action = na.pass)

#coef(summary(plotlm))

#plotlm.frame <- bind_rows(as.data.frame(coef(summary(plotlm))), as.data.frame(coef(summary(plotlm))))

#plotlm.frame <- plotlm.frame %>% 
  #separate(c("SPP", "YEAR", "SITEid"))
#how separate column with no name

#install.packages("~/Google Drive/My Drive/MEForLab",
                 #repos = NULL,
                 #type = "source")


