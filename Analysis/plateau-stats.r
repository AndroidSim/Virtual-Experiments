#! /usr/bin/Rscript
##
# Read multiple *.csv files and print averages and standard deviations.
# Script has been designed and tested on extRatio .csv files.
# Parsing of files utilizes code from cmp-plot.r.
# Script does not explicitly check for a valid range.
#
# Time-stamp: <2021-04-22 13:14:00 rck>
#


argv <- commandArgs(TRUE)

if (length(argv) < 3) {
    print("Usage: plateau-stats.r range_start range_end <analysis .csv file> <analysis .csv file>")
    print("  e.g. plateau-stats.r 3000 4000 vCul057_extRatio.csv vCul058_extRatio.csv")
    print("Note that this assumes data is in columns 2:4 and that a valid range is provided.")
    quit()
}

range.start <- as.numeric(argv[1])+1
range.end <- as.numeric(argv[2])+1
files <- argv[-1:-2]

data <- vector("list")
## get component name from basename of 1st argument
files.basename <- basename(files[1])
compname <- substr(files.basename,tail(gregexpr('_',files.basename)[[1]], n=1)+1,nchar(files.basename))
fileName.base <- substr(compname, 0, regexpr('(_|.csv)', compname)-1)
print(compname)
cat(paste(fileName.base,'\n'))
filenum <- 1
for (f in files) {
  nxtName <- substr(basename(f),0,tail(gregexpr('_',basename(f))[[1]], n=1)-1)
  raw <- read.csv(f)
  raw <- raw[c(range.start:range.end),]
  rawm <- colMeans(raw)
  cat("\n", nxtName, ": [" , range.start-1, ":", range.end-1, "]\n", sep='')
  cat("Averages \n")
  print(rawm[2:4])
  stdevs <- apply(raw,2,sd)
  cat("Standard deviations\n")
  print(stdevs[2:4])
  filenum <- filenum+1
}



