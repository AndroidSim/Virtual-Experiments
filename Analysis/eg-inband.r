#! /usr/bin/Rscript
###
## Sum the columns with the 1st element of the column name is within the band
## specified by [dMin,dMax).  Column elements are separated with ":".
##
## Time-stamp: <2019-11-11 18:12:09 gepr>
###

argv <- commandArgs(TRUE)

usage <- function() {
    print("Usage: eg-inband.r <dPV|dCV> <dMin> <dMax> <exp directories>")
    print("  exp directories should contain an <exp>_enzymes.csv file that looks like:")
    print("    Time, 0:Phase1, 0:Phase2, ...")
    print("  Should work for either dPV or dCV data produced by the ")
    print("  reduction scripts.")
    quit()
}

if (length(argv) < 4) usage()
direction <- argv[1]
dMin <- as.numeric(argv[2])
dMax <- as.numeric(argv[3])
exps <- argv[-(1:3)]

tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

for (exp in exps) {
  exp <- paste(dirname(exp), basename(exp), sep="/")
  inFile <- paste(exp, "_enzymes-", direction, ".csv", sep="")
  dat <- read.csv(inFile, check.names=F)

  sums <- sumBandByLastTag(dat, c(dMin,dMax))
  if (length(sums) >= 1)
    write.csv(sums,paste(exp,"_enzymes-",direction,"∈[",dMin,",",dMax,").csv",sep=""),row.names=F)
  else
    cat("No EGs in",direction,"∈[",dMin,",",dMax,")\n")
}

