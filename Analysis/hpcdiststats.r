#! /usr/bin/Rscript
###
## Time-stamp: <2019-08-13 15:09:35 gepr>
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
  print("Usage: hpcdiststats.r <exp directories>")
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
    sum.dist <- 0
    for (dir in directions) {
      row <- c(row, min(hi.hr[grep(dir,names(hi.hr))]),
               max(hi.hr[grep(dir,names(hi.hr))]),
               mean(hi.hr[grep(dir,names(hi.hr))]),
               median(hi.hr[grep(dir,names(hi.hr))]),
               sd(hi.hr[grep(dir,names(hi.hr))]))
      sum.dist <- sum.dist + hi.hr[grep(dir,names(hi.hr))]
    }
    row <- c(row,min(sum.dist),max(sum.dist),mean(sum.dist),median(sum.dist),sd(sum.dist))
    trial <- substr(f,regexpr("[0-9]+.json",f),regexpr(".json",f)-1)
    trials <- c(trials,trial)
    summ <- rbind(summ,row)

    setTxtProgressBar(pb,getTxtProgressBar(pb)+1) ## progress bar

  }
  close(pb)
  summ <- as.data.frame(summ, row.names=trials)
  sum.names <- paste("dPV+dCV.",parameters,sep="")
  colnames(summ) <- c(columns,sum.names)
  print(summ,digits=4)

  row <- vector()
  sum.dist <- 0
  for (dir in directions) {
    row <- c(row, min(allTrials[grep(dir,names(allTrials))]),
             max(allTrials[grep(dir,names(allTrials))]),
             mean(allTrials[grep(dir,names(allTrials))]),
             median(allTrials[grep(dir,names(allTrials))]),
             sd(allTrials[grep(dir,names(allTrials))]))
    sum.dist <- sum.dist + allTrials[grep(dir,names(allTrials))]
  }
  row <- c(row,min(sum.dist),max(sum.dist),mean(sum.dist),median(sum.dist),sd(sum.dist))
  whole <- rbind(whole, row)
}
rownames(whole) <- exps
sum.names <- paste("dPV+dCV.",parameters,sep="")
colnames(whole) <- c(columns,sum.names)
print(whole, digits=4)
