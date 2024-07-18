#load packages 
library(dplR)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(gdata)

## KI Spruce ##
KI<-read.tucson("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/raw/KI_PIRU2.raw", header = NULL, long = FALSE, encoding = getOption("encoding"), edge.zeros = TRUE)
KI_test<-read.csv("~/Google Drive/My Drive/Dendrochronology/KI_test.csv")
year<-KI_test$YEAR
rownames(KI_test)<-year
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

#detrend all at once
KI.rw<-detrend(rwl=KImeans, method="Friedman")
KIrs.crn<-chron(x=KI.rw, prefix = "KI RS", biweight = TRUE, prewhiten = TRUE)
write.csv(KIrs.crn, file="~/Google Drive/My Drive/Dendrochronology/KIRS_chron.csv")
plot.crn(x= KIrs.crn[,-1], add.spline = TRUE, nyrs = 20, f = 0.5, crn.line.col = "grey",
         spline.line.col = "red", samp.depth.col = "grey90", samp.depth.border.color = "grey80",
         crn.lwd = 1, spline.lwd=2.0, abline.pos=1, abline.col = "black",abline.lty=1,
         xlab="Time", ylab="RWI")
plot.crn(x=KIrs.crn, add.spline=TRUE, nyrs=NULL)

#combine
CHRONO.crn <- chron(chrono_i, prefix = "Plymouth")
plot(CHRONO.crn, add.spline=TRUE, nyrs=20, xlim=c(1900, 2025), ylim=c(0.5, 1.5))

#export
CHRONO.crn <- cbind(rownames(CHRONO.crn), data.frame(CHRONO.crn, row.names=NULL))
CHRONO.crn <- rename.vars(CHRONO.crn, c("rownames(CHRONO.crn)"), c("YEAR"))
write.table(format(CHRONO.crn, digits=2), file="KImeanscor.TXT", append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)

##KI Balsam fir##

KIbf<-read.tucson("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/raw/KI_ABBA.raw", header = NULL, long = FALSE, encoding = getOption("encoding"), edge.zeros = TRUE)
#KIbf<-write.csv(KIbf, "KIbf_fix.csv")
KIbf_fix<-read.csv("~/Google Drive/My Drive/Dendrochronology/KIbf_fix.csv")
year<-KIbf_fix$x
rownames(KI_test)<-year
KIbf_fix<-KIbf_fix[,2:38]

ids<-read.ids(KIbf_fix, stc=c(2, 4, 1))
KImeansbf <- treeMean(KIbf, ids, na.rm=TRUE)

KIbf.rw<-detrend(rwl=KImeansbf, method="Friedman")
KIbf.crn<-chron(x=KIbf.rw, prefix = "KI BF", biweight = TRUE, prewhiten = TRUE)
write.csv(KIrs.crn, file="~/Google Drive/My Drive/Dendrochronology/KIBF_chron.csv")
plot.crn(x= KIbf.crn[,-1], add.spline = TRUE, nyrs = 20, f = 0.5, crn.line.col = "grey",
         spline.line.col = "red", samp.depth.col = "grey90", samp.depth.border.color = "grey80",
         crn.lwd = 1, spline.lwd=2.0, abline.pos=1, abline.col = "black",abline.lty=1,
         xlab="Time", ylab="RWI")
plot.crn(x=KIrs.crn, add.spline=TRUE, nyrs=NULL)
chrono_bf=i.detrend(KImeansbf, nyrs = NULL, f=0.5, pos.slope = FALSE)


CHRONO.bf <- chron(chrono_bf, prefix = "Plymouth")
plot(CHRONO.bf, add.spline=TRUE, nyrs=20, xlim=c(1900, 2025), ylim=c(0.5, 1.5))

CHRONO.bf <- cbind(rownames(CHRONO.bf), data.frame(CHRONO.bf, row.names=NULL))
CHRONO.bf <- rename.vars(CHRONO.bf, c("rownames(CHRONO.bf)"), c("YEAR"))
write.table(format(CHRONO.crn, digits=2), file="KIbfcor.TXT", append = FALSE, quote = FALSE, sep = ",", row.names = FALSE, col.names = TRUE)


