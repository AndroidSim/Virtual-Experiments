#! /usr/bin/Rscript
###
## This script averages Enzyme Group capacities per dPV, per EG
##
## Time-stamp: <2019-07-11 09:46:28 gepr>
###

argv <- commandArgs(TRUE)

usage <- function() {
    print("Usage: eg-preproc.r <exp directories>")
    print("  directories should contain files named enzymes-dPV-[0-9]+.csv.gz")
    quit()
}

if (length(argv) < 1) usage()

tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

exps <- argv

filebase <- "enzymes"
## for each experiment
for (d in exps) {
  d <- paste(dirname(d), basename(d), sep="/") # removes trailing slash from the argument

  if (!file.exists(d)) {
    print(paste(d,"doesn't exist."))
    next
  }

  for (direction in c("dPV", "dCV")) {
    pattern <- paste(filebase, "-", direction, "-[0-9]+.csv.gz", sep="")
    cat(paste("Working on", pattern, "\n"))
    allfiles <- list.files(path = d, pattern = pattern, recursive=T, full.names=T)

    totals <- avgByColumn(allfiles)

    write.csv(totals,
              file=paste(d, "_", filebase, "-", direction, ".csv", sep=""),
              row.names=F)
  } ## end for (direction in paste("dPV", "dCV")) {
} ## end for (d in exps)


