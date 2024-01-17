library(terra)
map <- read.csv("C:/Users/michael.premer/Desktop/steMapTrial.csv")
map[is.na(map)] <- 0
head(map)
plot(map$X,map$Y)

library(tidyverse)
map <- dplyr::filter(map,X>0)

mapper <- map%>%
  mutate(k=1)

mapping <- map%>%
  full_join(mapper,by="k")%>%
  mutate(dist=sqrt((X.x-X.y)^2+(Y.x-Y.y)^2))%>%
  select(-k)

neighbors <- dplyr::filter(mapping,dist<30)

hood <- neighbors%>%
  group_by(TREE.x)%>%
  summarize(s)

mapping <- map%>%
  full_join(mapper, by = "k") %>% 
  filter(ID.x != ID.y) %>%
  mutate(dist = sqrt((X.x - X.y)^2 + (Y.x - Y.y)^2)) %>%
  select(-k)

head(mapping)


tree <- map$TREE
X <- map$X
Y <- map$Y

hose <- cbind(X,Y)


pts <- data.frame(tree=map$TREE,X=map$X,Y=map$Y)
head(pts)
pts[is.na(pts)] <- 0

dis <- dist(pts)
dis
