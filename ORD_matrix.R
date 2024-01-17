# 2023 November 6
# Example of generating a species matrix from a data frame listing Species IV by Plot and Installation (In this example, plot and stand)

ov.s <- read.csv("~/Desktop/Matrix_example.csv")
ov.s$uid <- paste(ov.s$Stand_ID,".",ov.s$Plot_ID)
ov.s <- ov.s[3:5]
library(reshape2)
molten <- melt(as.data.frame(ov.s),id=c("uid","IV"))
arf <- dcast(molten,uid~value,value.var = "IV")
arf[is.na(arf)] <- 0
arf$srata <- "overstory"
write.csv(arf,"~/Desktop/Overstory_matrix.csv")
dev.off()
rm(list=ls())
