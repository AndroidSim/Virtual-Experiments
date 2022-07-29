#! /usr/bin/Rscript

###
## Averages intraHepatocyte Solute counts per dPV and dCV, per Solute
## and extraCellular Solute, as represented in the celladj files.
##
## Time-stamp: <2019-07-11 09:47:23 gepr>
###

argv <- commandArgs(T)

usage <- function() {
  print("Usage: inextra-inband.r <dPV|dCV> dMin dMax <exp directories>")
  print("  directories should contain files like mobileObject-dPV-[0-9]+.csv.gz")
  print("  and celladj-d[CP]V-[0-9]+.csv.gz")
  quit()
}

if (length(argv) < 4) usage()

tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

direction <- argv[1]
dMin <- as.numeric(argv[2])
dMax <- as.numeric(argv[3])
exps <- argv[-(1:3)] ## all remaining args

filebases <- c("celladj", "mobileObject", "boundMObject", "entries", "exits", "rejects", "traps")

for (d in exps) {
  d <- paste(dirname(d), basename(d), sep="/") # removes trailing slash from the argument
  print(paste("Working on", d))

  if (!file.exists(d)) {
    print(paste(d,"doesn't exist."))
    next
  }

  ## for both measures
  for (measure in filebases) {

    pattern <- paste(measure,"-",direction,"-[0-9]+.csv.gz",sep="")
    allfiles <- list.files(path=d, pattern=pattern, recursive=T,  full.names=T)
    if (length(allfiles) <= 0) {
      print(paste("skipping ", measure, "-", direction, sep=""))
      next
    }
    totals <- avgByColumn(allfiles)
    fileprefix <- paste(d,"_",measure,"-",direction,sep="")
    write.csv(totals,
              file=paste(fileprefix,".csv",sep=""),
              row.names=F)
    distances <- unlist(strsplit(colnames(totals),":"))
    distances <- unique(distances[c(F,T)])
    Maxdist <- max(as.numeric(distances))
    if (dMin > Maxdist) {
      print(paste("dMin, ",dMin,", is > max distance, ",Maxdist,", for direction, ",direction))
      next
    }
    inband <- sumBandByLastTag(totals,c(dMin,dMax))
    write.csv(inband,
              file=paste(fileprefix,"∈[",dMin,",",dMax,").csv",sep=""),
              row.names=F)

  } ## end   for (measure in filebases) {
} ## end for (d in exps) {
