#! /usr/bin/Rscript
###
## Calculate peri-central to mid-zonal ratios.
##
## Time-stamp: <2019-02-25 11:21:13 gepr>
###
bandRadius <- 2

argv <- commandArgs(T)

usage <- function() {
  print("Usage: pc-to-mz.r pcMin pcMax <xpt_file>")
  print("  xpt_file should look like: exp_xpt.csv")
  print("  The script will derive the MZ interval from the xpt file.")
  print("  If either the PC file or the the MZ file don't exist, it will ")
  print("  execute the reduction script using those intervals.")
  print("  It then calls ratio.r for the PC/MZ calculation.")
  quit()
}

## get_Rscript_filename() stolen from getopt package:
## https://github.com/trevorld/r-getopt/blob/master/R/utils.R
## Copyright (c) 2011-2013 Trevor L. Davis <trevor.l.davis@gmail.com>
get_Rscript_filename <- function() {
    prog <- sub("--file=", "", grep("--file=", commandArgs(), value=TRUE)[1])
    if( .Platform$OS.type == "windows") {
        prog <- gsub("\\\\", "\\\\\\\\", prog)
    }
    prog
}

if (length(argv) != 3) usage()
error <- 0

pcMin <- argv[1]
pcMax <- argv[2]
xptFile <- argv[3]

scriptDir <- dirname(get_Rscript_filename())

xptData <- read.csv(xptFile, colClasses="numeric")
mzMin <- round(xptData[["dCV"]] - bandRadius)
mzMax <- round(xptData[["dCV"]] + bandRadius)

exp <- substr(basename(xptFile), 0, regexpr('_', basename(xptFile))-1)

for (type in c("celladj", "mobileObject")) {
  pcFile <- paste(exp, "-reduced/", exp, "_", type, "-dCV∈[", pcMin, ",", pcMax, ").csv", sep="")
  mzFile <- paste(exp, "-reduced/", exp, "_", type, "-dCV∈[", mzMin, ",", mzMax, ").csv", sep="")

  runPC <- !file.exists(pcFile)
  runMZ <- !file.exists(mzFile)
  if (runPC || runMZ) {
    ifelse(runPC, print(paste("Warning! ", pcFile, "missing.  Running reduction script.")),"")
    ifelse(runMZ, print(paste("Warning! ", mzFile, "missing.  Running reduction script.")),"")
    command <- paste(scriptDir, "/reduce-bands.sh", sep="")
    command <- paste(command, pcMin, pcMax, mzMin, mzMax, exp)
    error <- system(command)
  }

  runPC <- !file.exists(pcFile)
  runMZ <- !file.exists(mzFile)
  if (error != 0 || runPC || runMZ) {
    print("Error running reduction script.")
    quit()
  }

  ## Execute ratio.r
  command <- paste(scriptDir, "/ratio.r", " \"", pcFile, "\" \"", mzFile, "\"", sep="")
  error <- system(command)
  if (error != 0) {
    print("Error running ratio.r.")
    quit()
  }
  error <- system(paste("mv ", exp, "_*-ratio*.csv ", exp, "-reduced/", sep=""))
  if (error != 0) {
    print("Error moving ratio files.")
    quit()
  }

}
