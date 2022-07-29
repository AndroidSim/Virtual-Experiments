#! /usr/bin/Rscript

###
## Plot dCV/dPV vs Time.
##
## Time-stamp: <2019-05-14 11:29:38 gepr>
###
argv <- commandArgs(TRUE)

usage <- function() {
  print("Usage: dcvdpv.r pvMin pvMax cvMin cvMax <exp directories>")
  print("  directories should contain files like mobileObject-dPV-[0-9]+.csv.gz")
  quit()
}
is.infinite.data.frame <- function(obj){
    sapply(obj,FUN = function(x) is.infinite(x))
}

if (length(argv) < 5) usage()

pvMin <- argv[1]
pvMax <- argv[2]
dPV <- paste("dPV∈[", pvMin, ",", pvMax, ")", sep="")
cvMin <- argv[3]
cvMax <- argv[4]
dCV <- paste("dCV∈[", cvMin, ",", cvMax, ")", sep="")

exps <- argv[-(1:4)]

for (exp in exps) {
  exp <- paste(dirname(exp), basename(exp), sep="/") # removes trailing slash from the argument\

  types <- c("exposure", "mobileObject")
  for (type in types) {
    if (type == "exposure") {
      #numeratorsuffix <- paste("_entries-avg-pHPC-pMC-", dCV, "-exposure", sep="")
      #denominatorsuffix <- paste("_entries-avg-pHPC-pMC-", dPV, "-exposure", sep="")
      numeratorsuffix <- paste("_exposure-entries-", dCV, sep="")
      denominatorsuffix <- paste("_exposure-entries-", dPV, sep="")
      outfileroot <- paste("_entries-", type, "-", sep="")
    } else if (type == "mobileObject") {
      numeratorsuffix <- paste("_mobileObject-", dCV, sep="")
      denominatorsuffix <- paste("_mobileObject-", dPV, sep="")
      outfileroot <- "_mobileObject-"
    }
    
    numerfile <- paste(exp, "-reduced/", basename(exp), numeratorsuffix, ".csv",sep="")
    denomfile <- paste(exp, "-reduced/", basename(exp), denominatorsuffix, ".csv",sep="")
    nonfile <- F
    nodfile <- F
    if (!file.exists(numerfile)) {
      print(paste(numerfile,"doesn't exist."))
      nonfile <- T
    }
    if (!file.exists(denomfile)) {
      print(paste(denomfile,"doesn't exist."))
      nodfile <- T
    }
    if (nonfile || nodfile) {
      next
    } else {
      numerator <- read.csv(numerfile)
      denominator <- read.csv(denomfile)
	}
    ratio <- numerator[2:ncol(numerator)]/denominator[2:ncol(denominator)]
    ratio <- as.data.frame(cbind(numerator[,1],ratio))


    colnames(ratio) <- c("Time",colnames(ratio)[2:ncol(numerator)])
    ratio[is.na(ratio)] <- 0
    ratio[is.infinite.data.frame(ratio)] <- NA

    write.csv(ratio, paste(exp, "-reduced/", exp, outfileroot, dCV, "-to-", dPV, ".csv", sep=""), row.names=F)

  } # end for (type in types) {
} # end for (e in exps) {
