#! /usr/bin/Rscript
###
## Time-stamp: <2019-08-14 10:30:25 gepr>
##
## Gathers summary stats on RxnProbMap elements in hepinit-*.json
## files.  Write them for each trial and a final table with experiment
## totals.
##
## Note that the hepinit-????.json output contains the (potentially
## false) d[CP]V IVars (i.e. instance variables) inside the Hepatocyte
## objects.  If hepInitRead=true, then the d[CP]V values will not
## reflect the true values given the HepStruct used. However, these
## RxnProb stats will *still* be correct, commensurate with the d[CP]V
## values found in these hepinit-*.json files, because they are
## calculated based on what the Hepatocyte *thinks* its d[CP]V are.
##
###
argv <- commandArgs(T)

usage <- function() {
  print("Usage: hpcrpstats.r <exp directories>")
  print("  directories should contain files like hepinit-????.json")
  quit()
}

if (length(argv) < 1) usage()

library(matrixStats)

exps <- argv

filebase <- "hepinit"

parameters <- c("min","max","μ","η","σ")
columns <- vector()

whole <- vector()
for (d in exps) {
  cat(paste("Processing",d,"\n"))
  if (!file.exists(d)) {
    print(paste(d,"doesn't exist."))
    quit("no")
  }
  pattern <- paste(filebase,"-[0-9]+.json",sep="")
  files <- list.files(path=d, pattern=pattern, recursive=F, full.names=T)
  pb <- txtProgressBar(min=0,max=length(files)*5+1,style=3)
  setTxtProgressBar(pb,0)

  all <- numeric()
  summ <- vector()
  for (f in files) {
    trial <- numeric()
    trial.name <- substr(f,regexpr("[0-9]+.json",f),regexpr(".json",f)-1)
    hi <- jsonlite::fromJSON(f)
    hi.hr <- hi$hepatocyte_records
    setTxtProgressBar(pb,getTxtProgressBar(pb)+1) ## progress: read the file
    ## assume the enzyme groups are all the same for all vHPCs!
    if (!exists("egs")) egs <- names(hi.hr[[1]]$RxnProbMap)

    vHPC.ids <- names(hi.hr[1:length(hi.hr)])
    for (id in vHPC.ids) {
      prs <- hi.hr[id][[1]]$RxnProbMap
      row <- as.numeric(prs)
      trial <- rbind(trial, row)
      all <- rbind(all, row)
    }
    mins <- colMins(trial)
    setTxtProgressBar(pb,getTxtProgressBar(pb)+1) ## progress: mins
    maxs <- colMaxs(trial)
    setTxtProgressBar(pb,getTxtProgressBar(pb)+1) ## progress: maxs
    means <- colMeans(trial)
    setTxtProgressBar(pb,getTxtProgressBar(pb)+1) ## progress: means
    medians <- colMedians(trial)
    setTxtProgressBar(pb,getTxtProgressBar(pb)+1) ## progress: medians
    sds <- colSds(trial)
    setTxtProgressBar(pb,getTxtProgressBar(pb)+1) ## progress: sds
    trial.stats <- rbind(mins, maxs, means, medians, sds)
    trial.stats <- as.data.frame(trial.stats)
    rownames(trial.stats) <- parameters
    colnames(trial.stats) <- egs

    write.csv(trial.stats, paste(d, "_rxnprob-", trial.name, ".csv", sep=""), quote=F)

  }
  all.stats <- rbind(colMins(all), colMaxs(all), colMeans(all), colMedians(all), colSds(all))
  setTxtProgressBar(pb,getTxtProgressBar(pb)+1) ## progress bar

  all.stats <- as.data.frame(all.stats)
  rownames(all.stats) <- parameters
  colnames(all.stats) <- egs
  write.csv(all.stats, paste(d, "_rxnprob-alltrials.csv", sep=""), quote=F)

  close(pb)
}
