#load packages 
library(dplR)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(gdata)

setwd("~/Google Drive/My drive/CTRN_CFRU_Share/raw/raw")

## KI Spruce ##
KI<-read.tucson("KI_PIRU2.raw", header = NULL, long = FALSE, encoding = getOption("encoding"), edge.zeros = TRUE)
KI_test<-read.csv("~/Desktop/KI_test.csv")
KI_test<-KI_test[,2:39]
#data checking
dim(KI) 
colnames(KI)
class(KI)
rwl.report(KI)
summary(KI)
sens1(KI)


#prelim plotting
plot(KI)
#spag
plot(KI, plot.type="spag", col='blue')

#average acrosss cores in a single tree
ids<-read.ids(KI_test, stc=c(2, 4, 1))
KImeans <- treeMean(KI, ids, na.rm=TRUE)


#detrend the series
chrono_i=i.detrend(KImeans, nyrs = NULL, f=0.5, pos.slope = FALSE)

chrono_i$YR=row.names(chrono_i)
head(chrono_i)
colnames(chrono_i)
chrono_i

#combine
CHRONO.crn <- chron(chrono_i, prefix = "Plymouth")
plot(CHRONO.crn, add.spline=TRUE, nyrs=20, xlim=c(1900, 2025), ylim=c(0.5, 1.5))

#export
CHRONO.crn <- cbind(rownames(CHRONO.crn), data.frame(CHRONO.crn, row.names=NULL))
CHRONO.crn <- rename.vars(CHRONO.crn, c("rownames(CHRONO.crn)"), c("YEAR"))
write.table(format(CHRONO.crn, digits=2), file="KImeanscor.TXT", append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)

##KI Balsam fir##

KIbf<-read.tucson("KI_ABBA.raw", header = NULL, long = FALSE, encoding = getOption("encoding"), edge.zeros = TRUE)
#KIbf<-write.csv(KIbf, "KIbf_fix.csv")
KIbf_fix<-read.csv("~/Desktop/KIbf_fix.csv")
KIbf_fix<-KIbf_fix[,2:38]

ids<-read.ids(KIbf_fix, stc=c(2, 4, 1))
KImeansbf <- treeMean(KIbf, ids, na.rm=TRUE)

chrono_bf=i.detrend(KImeansbf, nyrs = NULL, f=0.5, pos.slope = FALSE)

CHRONO.bf <- chron(chrono_bf, prefix = "Plymouth")
plot(CHRONO.bf, add.spline=TRUE, nyrs=20, xlim=c(1900, 2025), ylim=c(0.5, 1.5))

CHRONO.bf <- cbind(rownames(CHRONO.bf), data.frame(CHRONO.bf, row.names=NULL))
CHRONO.bf <- rename.vars(CHRONO.bf, c("rownames(CHRONO.bf)"), c("YEAR"))
write.table(format(CHRONO.crn, digits=2), file="KIbfcor.TXT", append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)


