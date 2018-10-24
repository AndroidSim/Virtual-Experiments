#! /usr/bin/Rscript

###
## Averages intraHepatocyte Solute counts per dPV and dCV, per Solute
## and extraCellular Solute, as represented in the celladj files.
##
## Time-stamp: <2018-02-14 14:56:24 gepr>
###

argv <- commandArgs(T)

usage <- function() {
  print("Usage: inextra-inband.r dMin dMax <exp directories>")
  print("  directories should contain files like hsolute-dPV-[0-9]+.csv.gz")
  print("  and celladj-d[CP]V-[0-9]+.csv.gz")
  quit()
}

if (length(argv) < 3) usage()

source("~/R/misc.r")

dMin <- as.numeric(argv[1])
dMax <- as.numeric(argv[2])
exps <- argv[-(1:2)] ## all remaining args

filebases <- c("celladj", "hsolute", "entries", "exits", "rejects", "traps")

for (d in exps) {
  print(paste("Working on", d))

  if (!file.exists(d)) {
    print(paste(d,"doesn't exist."))
    next
  }

  ## for both measures
  for (measure in filebases) {
    ## for both directions
    for (dir in c("dPV","dCV")) {
		pattern <- paste(measure,"-",dir,"-[0-9]+.csv.gz",sep="")
		allfiles <- list.files(path=d, pattern=pattern, recursive=T,  full.names=T)
		if (length(allfiles) <= 0) {
			print(paste("skipping ", measure, "-", dir, sep=""))
			next
		}
		totals <- avgByColumn(allfiles)
		fileprefix <- paste(d,"_",measure,"-",dir,sep="")
		write.csv(totals,
                file=paste(fileprefix,".csv",sep=""),
                row.names=F)
		distances <- unlist(strsplit(colnames(totals),":"))
		distances <- unique(distances[c(F,T)])
		Maxdist <- max(as.numeric(distances))
		if (dMin > Maxdist) {
			print(paste("dMin, ",dMin,", is > max distance, ",Maxdist,", for direction, ",dir))
			next
		}
		inband <- sumBandByLastTag(totals,c(dMin,dMax))
		write.csv(inband,
                file=paste(fileprefix,"∈[",dMin,",",dMax,").csv",sep=""),
                row.names=F)
    }
  }
}
