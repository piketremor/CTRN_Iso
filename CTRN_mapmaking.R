##Script to make map of all 15 CTRN sites
dev.off()
rm(list=ls())

library(readxl)
library(oce)
library(terra)
library(dplyr)
library(lattice)
library(nlme)
#library(forfun)
library(sf)
library(sp)
library(raster)
library(maptools)
library(rgeos)
library(ggplot2)
library(gstat)
library(snow)
library(rpart)
library(Cubist)
library(fBasics)
library(nortest)
library(devtools)
library(randomForest)
#library(ithir)
library(USAboundaries)
library(equivalence)
library(pdp)
library(ggmap)
library(tidyr)
library(nlme)
library(forfun)
require(MEForLab)

locs<-read.csv("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv/CTRN_Plot_Location.csv")
test<-utm2lonlat(locs$Easting_X, locs$Northing_Y, zone = 19, hemisphere = "N", km = FALSE)
testy<-as.data.frame(test)                            
locs2<-cbind(locs,testy) 
flacs<-locs2[,c(1:2,5:6)]

detach(package:plyr)
flacs2 <- flacs%>%
  group_by(SiteID)%>%
  summarize(mean.lat = mean(latitude),
            mean.long = mean(longitude))


state_names <- "Maine"
states <- us_states(resolution="high",states=state_names)
plot(states$geometry)


ggm1 <- ggplot()+
  geom_sf(data=states,fill="white")+
  theme_classic()+
  geom_point(data=flacs2,aes(x=mean.long,y=mean.lat,group=SiteID,colour=SiteID,stroke=1.5),
             size=2,alpha=1)+
  theme(axis.line = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+
  theme(axis.text = element_blank())+
  guides(colour=guide_legend(title="CTRN Sites"))
ggm1


####Google maps#####

register_google(key="AIzaSyBRj7TixooXrAm05JQwPfoPhjuc-3aIsnk")

ggmap(get_googlemap(center=c(lon=-69,lat=45.5),
                    zoom=7,scale=4,color='bw'))+
  geom_point(data=flacs2,aes(x=mean.long,y=mean.lat,group=SiteID,colour=SiteID,stroke=1.5),
             size=2,alpha=1)
