#! /usr/bin/Rscript

###
## Plot dCV/dPV vs Time.
##
## Time-stamp: <2018-10-16 10:24:46 gepr>
###
argv <- commandArgs(TRUE)

dCV <- "dCV∈[0,8)"
dPV <- "dPV∈[5,10)"
types <- c("exposure", "hsolute")

is.infinite.data.frame <- function(obj){
    sapply(obj,FUN = function(x) is.infinite(x))
}

exp <- argv[1]

for (type in types) {
  if (type == "exposure") {
    numeratorsuffix <- paste("_entries-avg-pHPC-pMC-", dCV, "-exposure", sep="")
    denominatorsuffix <- paste("_entries-avg-pHPC-pMC-", dPV, "-exposure", sep="")
    outfileroot <- paste("_entries-", type, "-", sep="")
  } else if (type == "hsolute") {
    numeratorsuffix <- paste("_hsolute-", dCV, sep="")
    denominatorsuffix <- paste("_hsolute-", dPV, sep="")
    outfileroot <- "_hsolute-"
  }
  numerator <- read.csv(paste(exp, "-reduced/", exp, numeratorsuffix, ".csv",sep=""))
  denominator <- read.csv(paste(exp, "-reduced/", exp, denominatorsuffix, ".csv",sep=""))

  ratio <- numerator[2:ncol(numerator)]/denominator[2:ncol(denominator)]
  ratio <- as.data.frame(cbind(numerator[,1],ratio))


  colnames(ratio) <- c("Time",colnames(ratio)[2:ncol(numerator)])
  ratio[is.na(ratio)] <- 0
  #print(head(ratio,n=20))
  ratio[is.infinite.data.frame(ratio)] <- NA
  #print(head(ratio,n=20))

  write.csv(ratio, paste(exp, "-reduced/", exp, outfileroot, dCV, "-to-", dPV, ".csv", sep=""), row.names=F)

}
