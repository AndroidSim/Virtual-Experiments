#! /usr/bin/Rscript
###
## Time-stamp: <2019-08-14 12:01:04 gepr>
##
## Gathers summary stats on RxnProbMap elements in hepinit-*.json
## files.  Write them for each dPV and dCV.
##
## Note that the hepinit-????.json output contains the (potentially
## false) d[CP]V IVars (i.e. instance variables) inside the Hepatocyte
## objects.  If hepInitRead=true, then the d[CP]V values will not
## reflect the true values given the HepStruct used. However, the
## RxnProb stats will be commensurate with these d[CP]V because they
## are calculated based on what the Hepatocyte *thinks* its d[CP]V
## are.
##
###
argv <- commandArgs(T)

usage <- function() {
  print("Usage: hpcrppdstats.r <exp directories>")
  print("  directories should contain files like hepinit-????.json")
  quit()
}

if (length(argv) < 1) usage()

library(matrixStats)

exps <- argv

filebase <- "hepinit"

directions <- c("dPV", "dCV")
parameters <- c("min","max","μ","η","σ")
columns <- vector()

for (d in exps) {
  cat(paste("Processing",d,"\n"))
  if (!file.exists(d)) {
    print(paste(d,"doesn't exist."))
    quit("no")
  }
  pattern <- paste(filebase,"-[0-9]+.json",sep="")
  files <- list.files(path=d, pattern=pattern, recursive=F, full.names=T)

  ## load all the trials into memory, extract max dists and EG names
  all <- vector("list")
  maxd <- vector("list")
  for (f in files) {
    trial <- numeric()
    trial.name <- substr(f,regexpr("[0-9]+.json",f),regexpr(".json",f)-1)
    cat(paste("  Loading", f,"\n"))
    all[[trial.name]] <- jsonlite::fromJSON(f)
    hi <- all[[trial.name]]
    hi.hr <- hi$hepatocyte_records

    vHPC.ids <- names(hi.hr[1:length(hi.hr)])
    ## first get the maximum d[CP]V values
    hi.hr.unlisted <- unlist(hi.hr)
    maxdPV <- max(hi.hr.unlisted[grep("dPV",names(hi.hr.unlisted))])
    maxd[[directions[1]]] <- max(maxd[[directions[1]]], maxdPV)
    maxdCV <- max(hi.hr.unlisted[grep("dCV",names(hi.hr.unlisted))])
    maxd[[directions[2]]] <- max(maxd[[directions[2]]], maxdCV)

    ## assume the enzyme groups are all the same for all vHPCs!
    if (!exists("egs")) egs <- names(hi.hr[[1]]$RxnProbMap)

  } ## done reading all the trial files

  ## write 2 (d[CP]V) files for each EG
  for (eg in egs) {
    cat(paste(" ", eg))
    for (direction in directions) {
      outfilename <- paste(d, "_rxnprob-",direction,"-", eg, ".csv", sep="")
      if (exists("egstats")) rm("egstats")
      cols <- seq(0,maxd[[direction]])
      for (col in cols) {
        egvals <- vector()
        for (trial in all) {
          trial.hr <- trial$hepatocyte_records
          for (hpcnum in 1:length(trial.hr)) {
            ## if this hpc is at this dPV, include it
            if (trial.hr[[hpcnum]][direction] == col) {
              egvals <- c(egvals, trial.hr[[hpcnum]]$RxnProbMap[eg])
            } ## at this dPV
          } ## ∀ hpcs
        } ## ∀ trials
        egvals <- as.numeric(egvals)
        colstats <- c(min(egvals), max(egvals), mean(egvals), median(egvals), sd(egvals))
        if (exists("egstats")) egstats <- cbind(egstats, colstats)
        else egstats <- as.data.frame(colstats)
      } ## ∀ distance
      colnames(egstats) <- cols
      rownames(egstats) <- parameters
      write.csv(egstats, outfilename, quote=F)
    } ## ∀ direction
  }
  cat("\n")

}
