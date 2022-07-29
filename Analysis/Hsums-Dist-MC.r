#! /usr/bin/Rscript

###
## Calculates the sum of Hepatocytes for each dCV and dPV bands over all MC trials
## using the raw data of the number of Hepatoctyes at each distance in the hcount files.
##
## Time-stamp: <2019-11-11 18:27:16 gepr>
###

argv <- commandArgs(T)

require(stats) # for statistics

usage <- function() {
  print("Usage: Measure_per_HPC-inband.r <dCV|dPV> <dMin> <dMax> <exp directories>")
  print("  directories should contain files like mobileObject-d[CP]V-[0-9]+.csv.gz, ")
  print (" celladj-d[CP]V-[0-9]+.csv.gz, and hcount-d[CP]V-[0-9]+.csv")
  quit()
}

if (length(argv) < 4) usage()

tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

direction <- argv[1]
dMin <- as.numeric(argv[2])
dMax <- as.numeric(argv[3])
exps <- argv[-(1:3)] ## all remaining args

Hfb <- "hcount"

for (x in exps) {
  x <- paste(dirname(x), basename(x), sep="/") # removes trailing slash from the argument
  print(paste("Working on", x))

  if (!file.exists(x)) {
    print(paste(x,"doesn't exist."))
    next
  }

  hfile <- paste(Hfb,"-",direction,"-[0-9]+.csv",sep="")
  allhfiles <- list.files(path=x, pattern=hfile, recursive=T,  full.names=T)
  nHfiles <- length(allhfiles)

  if (nHfiles <= 0) {
    print(paste("no hepatocyte counts for exp ", x))
    next
  }

  nMCtrials <- nHfiles
  ## progress bar
  print(paste("distance = ",direction))
  pb <- txtProgressBar(min=0,max=nMCtrials,style=3)
  setTxtProgressBar(pb,0)

  ## loop over each MC trial
  for (i in 1:nMCtrials) {
    ## read hepatocyte count data from file
    ## first row (distances) becomes column names when read
    ##cat("\n")
    Hdat <- read.csv(allhfiles[i], check.names=F)
    ##print(Hdat)
    ##cat("\n")
    ##Hdat <- read.csv(allhfiles[i], colClasses="numeric")
    ##print(Hdat)

    ## two sets of data for each measure and H count
    ## one = totals over whole lobule, other = totals within distance band
    ## select measure data within band
    HbandMC <- Hdat[,dMin<=as.numeric(colnames(Hdat)) & as.numeric(colnames(Hdat))<dMax]

    ## if there is no data for this band, place a column of zeros
    if (ncol(as.matrix(HbandMC)) < 1) HbandMC <- as.data.frame(t(seq(0,ncol(Hdat))))
    ## if there is only 1 column for this band, convert it to a 1 column DF
    if (ncol(as.matrix(HbandMC)) == 1) HbandMC <- as.data.frame(t(HbandMC))

    ## calculate the total number of Hepatocytes and within the band
    if (exists("HsumbandMC")) {
       if (ncol(HbandMC) > 1) {
         HsumbandMC <- cbind(HsumbandMC, rowSums(HbandMC))
       } else {
         HsumbandMC <- cbind(HsumbandMC, HbandMC)
      }
    } else {
      if (length(HbandMC) > 1) {
        HsumbandMC <- rowSums(HbandMC)
      } else {
        HsumbandMC <- HbandMC
      }
    }
    if (exists("HsumallMC")) {
      HsumallMC <- cbind(HsumallMC,rowSums(Hdat))
    } else {
      HsumallMC <- rowSums(Hdat)
    }

    ## set data from 1st file as totals, then sum with subsequent data
    if (i == 1) {
      Htotals <- Hdat
    } else {
      ## combine H counts by distance from each trial
      Htotals <- pad1stColumns(Htotals, Hdat)
      Hpad <- pad1stColumns(Hdat, Htotals)
      Htotals <- Htotals[,order(names(Htotals))]
      Hpad <- Hpad[,order(names(Hpad))]
      Htotals <- Htotals + Hpad
    }
    setTxtProgressBar(pb,i); ## progress bar
  } ## loop over MC trials
  setTxtProgressBar(pb,i); ## progress bar
  close(pb) ## progress bar

  ## two sets of data for H count
  ## one = totals over whole lobule, other = totals within distance band
  Hinband <- Htotals[,dMin<=as.numeric(colnames(Htotals)) & as.numeric(colnames(Htotals))<dMax]
  HsumallMC <- as.data.frame(as.matrix(HsumallMC))
  colnames(HsumallMC) <- paste(1:nMCtrials)
  colnames(HsumbandMC) <- paste(1:nMCtrials)

  ## write hepatocyte files, pD = per Distance, pT = per MC trial
  Hband_pT <- HsumbandMC
  Htotal_pT <- HsumallMC
  Hband_pD <- Hinband
  Htotal_pD <- Htotals

  fileprefix <- paste(x, "_", "Hsums-pMC-", direction, sep="")
  write.csv(cbind(Htotal_pT, total = rowSums(Htotal_pT)), file=paste(fileprefix, ".csv", sep=""), row.names=F)

  fileprefix <- paste(x, "_", "Hsums-pDist-", direction, sep="")
  write.csv(cbind(Htotal_pD, total = rowSums(Htotal_pD)), file=paste(fileprefix, ".csv", sep=""), row.names=F)

  fileprefix <- paste(x, "_", "Hsums-pMC-", direction, sep="")
  write.csv(cbind(Hband_pT, total = rowSums(Hband_pT)), file=paste(fileprefix, "∈[", dMin, ",", dMax, ").csv", sep=""), row.names=F)

  fileprefix <- paste(x, "_", "Hsums-pDist-", direction, sep="")
  write.csv(cbind(Hband_pD, total = ifelse(ncol(Hband_pD) < 1, Hband_pD, rowSums(Hband_pD))), file=paste(fileprefix, "∈[",dMin, ",", dMax, ").csv", sep=""), row.names=F)

  ## remove variables for next measure distance pair
  remove(HsumbandMC)
  remove(HsumallMC)
  remove(Hinband)
  remove(Htotals)

} # loop over experiments
