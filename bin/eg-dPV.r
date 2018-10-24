#! /usr/bin/Rscript
#! /bin/bash setR
argv <- commandArgs(TRUE)
#dev.off()

###
## This script averages Enzyme Group capacities per dPV, per EG
##
## gepr 2017-03-24
###


usage <- function() {
    print("Usage: eg-dPV.r <exp directories>")
    print("  directories should contain files named enzymes-dPV-[0-9]+.csv.gz")
    quit()
}

if (length(argv) < 1) usage()

source("~/R/misc.r")

exps <- argv


filebase <- "enzymes"
# for each experiment
for (d in exps) {
  d <- unlist(strsplit(d,'/')) ## remove trailing slash if present
  print(paste("Working on", d))

  if (!file.exists(d)) {
    print(paste(d,"doesn't exist."))
    quit("no")
  }

  allfiles <- list.files(path = d, pattern = "enzymes-dPV-[0-9]+.csv.gz", recursive=T, full.names=T)

  totals <- avgByColumn(allfiles)

  write.csv(totals,
    file=paste(d,"_",filebase,".csv",sep=""),
    row.names=F)

} ## end for (d in exps)

#q()
