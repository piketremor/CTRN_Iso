#load packages 
library(dplR)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(gdata)

setwd("~/Google Drive/My drive/CTRN_CFRU_Share/raw/raw")

## KI Spruce
KI<-read.tucson("KI_PIRU2.raw", header = NULL, long = FALSE, encoding = getOption("encoding"), edge.zeros = TRUE)
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

#detrend the series
chrono_i=i.detrend(KI, nyrs = NULL, f=0.5, pos.slope = FALSE)
chrono_i$YR=row.names(chrono_i)
head(chrono_i)
colnames(chrono_i)
chrono_i
