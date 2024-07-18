rm(list=ls())
library(dplyr)
# generate inter-stem distances from CTRN dataset
# Premer, M.I. 3/4/2024
setwd("~/Google Drive/My drive/CTRN_CFRU_Share/raw/csv")
trees <- read.csv("Trees2023.csv")
locs <- read.csv("Tree_locations_species.csv")
names(locs)
ex <- locs
ex <- dplyr::filter(locs,SITEid=="GR",PLOTid=="4")
ex[is.na(ex)] <- 0
ex<- dplyr::filter(ex,DIST>0)
ki <- ex
ki <- ki%>%
  mutate(k=1)
ki <- ki%>% 
  full_join(ki, by = "k") %>% 
  #filter(TREE.x != TREE.y) %>%
  mutate(dist = sqrt((X.x - X.y)^2 + (Y.x - Y.y)^2)) %>%
  select(-k)
new.df <- data.frame(SITEid=ki$SITEid.x,PLOTid=ki$PLOTid.x,focal.tree=ki$TREE.x,competitor.tree=ki$TREE.y,distance=ki$dist)
new.df <- arrange(new.df,SITEid,PLOTid,focal.tree,competitor.tree)
new.df <- dplyr::filter(new.df,distance>0)
plot(density(new.df$distance))
#View(new.df)
write.csv(new.df,"~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv/GR_stemMap_Plot4.csv")
dev.off()
rm(list=ls())
