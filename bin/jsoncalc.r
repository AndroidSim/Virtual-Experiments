#! /usr/bin/Rscript
###
## Time-stamp: <2018-05-22 11:29:08 gepr>
##
## Gathers summary stats on hepinit-*.json files.  Print them for each trial
## and a final table with experiments totals.
##
## Note that the hepinit-????.json output contains the (potentially false)
## IVars (i.e. instance variables) inside the Hepatocyte objects.
## So if hepInitRead=true, then these
## stats will not be representative of the actual dPV/dCV values.
##
###

argv <- commandArgs(T)

usage <- function() {
  print("Usage: jsoncalc.r <exp directories>")
  print("  directories should contain files like hepinit-????.json")
  quit()
}

if (length(argv) < 1) usage()

exps <- argv

filebase <- "hepinit"

directions <- c("dPV","dCV")
parameters <- c("min","max","μ","η","σ")
columns <- vector()
for (dir in directions)
  for (param in parameters)
    columns <- c(columns,paste(dir,".",param,sep=""))

whole <- vector()
for (d in exps) {
  print(paste("Processing",d))
  if (!file.exists(d)) {
    print(paste(d,"doesn't exist."))
    quit("no")
  }
  pattern <- paste(filebase,"-[0-9]+.json",sep="")
  files <- list.files(path=d, pattern=pattern, recursive=F, full.names=T)
  pb <- txtProgressBar(min=0,max=length(files),style=3)
  setTxtProgressBar(pb,0)

  allTrials <- vector()

  summ <- vector()
  trials <- vector()
  for (f in files) {
    hi <- jsonlite::fromJSON(f)
    hi.hr <- unlist(hi$hepatocyte_records)
    allTrials <- c(allTrials, hi.hr)
    row <- vector()
    for (dir in directions) {
      row <- c(row, min(hi.hr[grep(dir,names(hi.hr))]),
               max(hi.hr[grep(dir,names(hi.hr))]),
               mean(hi.hr[grep(dir,names(hi.hr))]),
               median(hi.hr[grep(dir,names(hi.hr))]),
               sd(hi.hr[grep(dir,names(hi.hr))]))
    }
    trial <- substr(f,regexpr("[0-9]+.json",f),regexpr(".json",f)-1)
    trials <- c(trials,trial)
    summ <- rbind(summ,row)

    setTxtProgressBar(pb,getTxtProgressBar(pb)+1) ## progress bar

  }
  close(pb)
  summ <- as.data.frame(summ, row.names=trials)
  colnames(summ) <- columns
  print(summ,digits=4)

  row <- vector()
  for (dir in directions) {
    row <- c(row, min(allTrials[grep(dir,names(allTrials))]),
             max(allTrials[grep(dir,names(allTrials))]),
             mean(allTrials[grep(dir,names(allTrials))]),
             median(allTrials[grep(dir,names(allTrials))]),
             sd(allTrials[grep(dir,names(allTrials))]))
  }
  whole <- rbind(whole, row)
}
rownames(whole) <- exps
colnames(whole) <- columns
print(whole, digits=4)
