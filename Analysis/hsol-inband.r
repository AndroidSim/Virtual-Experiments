#! /usr/bin/Rscript
##! /bin/bash setR
argv <- commandArgs(T)
##dev.off()

###
## Averages Hepatocyte Solute counts per dPV and dCV, per Solute.
##
## gepr 2017-06-15
###

usage <- function() {
  print("Usage: hsol-inband.r dMin dMax <exp directories>")
  print("  directories should contain files like hsolute-dPV-[0-9]+.csv.gz")
  quit()
}

if (length(argv) < 3) usage()

tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

dMin <- as.numeric(argv[1])
dMax <- as.numeric(argv[2])
exps <- argv[-(1:2)] ## all remaining args

filebase <- "hsolute"
for (d in exps) {
  print(paste("Working on", d))

  if (!file.exists(d)) {
    print(paste(d,"doesn't exist."))
    quit("no")
  }

  ## for both directions
  for (dir in c("dPV","dCV")) {
    pattern <- paste(filebase,"-",dir,"-[0-9]+.csv.gz",sep="")
    allfiles <- list.files(path=d, pattern=pattern, recursive=T,  full.names=T)
    totals <- avgByColumn(allfiles)
    fileprefix <- paste(d,"_",filebase,"-",dir,sep="")
    write.csv(totals,
              file=paste(fileprefix,".csv",sep=""),
              row.names=F)

    inband <- sumBandByLastTag(totals,c(dMin,dMax))
    write.csv(inband,
              file=paste(fileprefix,"∈[",dMin,",",dMax,").csv",sep=""),
              row.names=F)
  }

}
