library(tidyverse)
library(MEForLab)
#first a demo
demo <- read.csv("~/Desktop/buck_demo.csv")
merch.demo <- Kozak_merchR(demo)


#now try on a subset of CTRN data
trees<-read.csv("~/Google Drive/My Drive/CTRN_CFRU_Share/raw/csv/CTRNoverstorytrees.csv")[,2:10]
#select out just stand, plot, dbh, and height
trees<-filter(trees, YEAR==2018)
barky<-trees[,c(1:3,5:6,9)]
#just try WS for now
barkws<-filter(barky, SPP=="WS")
colnames(barkws)[colnames(barkws) == 'TOT_HT'] <- 'tht'
colnames(barkws)[colnames(barkws) == 'DBH'] <- 'dbh'
head(barkws)
barkna<-na.omit(barkws)
merch.ws<-KozakmerchWS(barkna)
merch.ws<-Kozak_merchR(barkna) #original





#Kozak function with parameters for white spruce

KozakmerchWS<-function (x, log.breaks = c(2, 5, 7, 8, 12, 51), log.grades = c("junk", 
                                                                               "pulp", "cns", "ssl", "lsl", "psl"), display.stems = FALSE) 
{
  if (!is.data.frame(x)) {
    stop("x needs to be a data.frame")
  }
  grade.vols <- as.data.frame(matrix(rep(0, nrow(x)), nrow(x), 
                                     length(log.grades), byrow = TRUE))
  names(grade.vols) <- log.grades
  bcmof.diDBH <- function(dbh, tht, cr, hi) {
    dbh <- dbh
    tht <- tht
    hi <- hi
    if (hi >= tht) {
      retval <- 0
    }
    else {
      a0 <- 0.732
      a1 <- 0.958
      a2 <- 0.159
      b1 <- 0.264
      b2 <- -0.4246
      b3 <- 0.551
      b4 <- -0.1269
      b5 <- 0.115
      b6 <- -0.6249
      XI <- (1 - ((hi/tht))^(1/3))/(1 - (4.5/tht)^(1/3))
      QI <- 1 - ((hi/tht)^(1/3))
      ZI <- (hi/tht)
      partA <- (a0 * (dbh^a1) * (tht^a2))
      partB.1 <- b1 * (ZI^4)
      partB.2 <- b2 * 1/(exp(dbh/tht))
      partB.3 <- b3 * ((1 - (ZI))/(1 - (4.5/tht)))^0.1
      partB.4 <- b4 * (1/dbh)
      partB.5 <- b5 * (tht^QI)
      partB.6 <- b6 * XI
      PartB <- partB.1 + partB.2 + partB.3 + partB.4 + 
        partB.5 + partB.6
      retval <- partA * (XI^PartB)
    }
    retval
  }
  merch.height.func <- function(hi, dbh, tht, cr, md) {
    dib <- bcmof.diDBH(dbh, tht, cr, hi)
    diff <- dib - md
    diff
  }
  scribner <- function(d, l) {
    scribner <- (0.79 * d^2 - 2 * d - 4) * l/16
    scribner
  }
  logs <- NULL
  for (s in 1:nrow(x)) {
    if (x[s, ]$dbh <= 5) 
      next
    mh.bks <- rep(0, length(log.breaks))
    for (i in 1:length(log.breaks)) {
      dbh <- x[s, ]$dbh
      tht <- x[s, ]$tht
      if (log.breaks[i] <= bcmof.diDBH(dbh, tht, cr, 0)) {
        mh <- uniroot(merch.height.func, c(0, tht), 
                      dbh = dbh, tht = tht, cr = 0.6, md = log.breaks[i])
        mh.bks[i] <- mh$root
      }
      else {
        mh.bks[i] <- 0
      }
    }
    vol.bks <- data.frame(log.grades, log.breaks, mh.bks)
    vol.bks.temp <- vol.bks[mh.bks > 0.3, ]
    dib <- bcmof.diDBH(dbh, tht, 0.6, 0.3)
    last.grade <- log.grades[nrow(vol.bks.temp) + 1]
    last.grade.df <- data.frame(log.grades = last.grade, 
                                log.breaks = dib, mh.bks = 0.3)
    vol.bks <- rbind(vol.bks.temp, last.grade.df)
    vol.bks.null <- data.frame(log.grades, log.breaks, mh.bks)
    n.sorts <- nrow(vol.bks)
    diff.rows <- nrow(vol.bks.null) - nrow(vol.bks)
    vol.bks <- rbind(vol.bks, tail(vol.bks.null, diff.rows))
    sed <- c(0, vol.bks[1:(n.sorts - 1), ]$log.breaks, rep(0, 
                                                           diff.rows))
    led <- c(vol.bks[1:n.sorts, ]$log.breaks, rep(0, diff.rows))
    vol.bks <- cbind(vol.bks, sed, led)
    log.lens <- c(tht - vol.bks[1, ]$mh.bks, abs(diff(vol.bks$mh.bks))[1:n.sorts - 
                                                                         1], rep(0, diff.rows))
    vol.bks <- cbind(vol.bks, log.lens)
    vol.bks$sm.vol <- scribner(vol.bks$sed, vol.bks$log.lens)
    rownames(vol.bks) <- 1:nrow(vol.bks)
    if (display.stems) {
      hi <- 0:tht
      dbh <- rep(dbh, length(hi))
      tht <- rep(tht, length(hi))
      cr <- rep(0.6, length(hi))
      stem.tpr <- data.frame(cbind(dbh, tht, cr, hi))
      stem.tpr$dib <- NA
      for (i in 1:nrow(stem.tpr)) {
        stem.tpr[i, ]$dib <- bcmof.diDBH(stem.tpr[i, 
        ]$dbh, stem.tpr[i, ]$tht, stem.tpr[i, ]$cr, 
        stem.tpr[i, ]$hi)
      }
      plot(stem.tpr$dib ~ stem.tpr$hi, type = "l")
      abline(v = vol.bks$mh.bks, lty = 2, lwd = 2)
      abline(h = vol.bks$log.breaks, lty = 3)
      abline(h = 0)
      text(x = vol.bks$mh.bks + 1.5, y = vol.bks$log.breaks + 
             1, labels = paste(vol.bks$log.breaks, "@", round(vol.bks$mh.bks, 
                                                              2)))
    }
    grade.vols[s, ] <- vol.bks$sm.vol
  }
  vol <- rowSums(grade.vols)
  x <- cbind(x, vol, grade.vols)
  x
}
