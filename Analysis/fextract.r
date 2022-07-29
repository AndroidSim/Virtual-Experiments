#! /usr/bin/Rscript
###
## Time-stamp: <2019-07-11 09:51:48 gepr>
##
## Average rawInput and rawOutput.  Then calculate:
##
##   Output Fraction: output/input
##   Extraction Ratio: (input - output)/input
##
###

PREWINDOW <- 31

usage <- function() {
    print(paste("Usage: fextract.r <exp directories>"))
    print("  Requires there be files below the PWD with a common name and")
    print("  common headers for columns. e.g.:")
    print("  ./2012-08-22-1645/rawInput-00[0-9][0-9].csv")
    quit()
}

argv <- commandArgs(TRUE)
if (length(argv) <= 0) usage()

exps <- argv

tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

profileError <- function (innames,outnames) {
  print("Column names don't match.")
  print(paste("colnames(ins) =",innames))
  print(paste("colnames(outs) =",outnames))
  quit()
}

outFract <- function(ins, outs) {
  of <- outs/ins
  of[1] <- ins[1]
  of[is.nan(of)] <- 0.0
  of[is.infinite(of)] <- 0.0
  of[is.na(of)] <- 0.0
  if (!all(colnames(ins) == colnames(outs))) profileError(colnames(ins),colnames(outs))
  colnames(of) <- colnames(ins)
  return(of)
}
extRatio <- function(ins, outs) {
  er <- (ins-outs)/ins
  er[1] <- ins[1]
  er[is.nan(er)] <- 0.0
  er[is.infinite(er)] <- 0.0
  er[is.na(er)] <- 0.0
  if (!all(colnames(ins) == colnames(outs))) profileError(colnames(ins),colnames(outs))
  colnames(er) <- colnames(ins)
  return(er)
}


# for each experiment
for (exp in exps) {
  exp <- paste(dirname(exp), basename(exp), sep="/") # removes trailing slash from the argument

  print (paste("exp = ",exp))

  for (pm in c("rawInput", "rawOutput")) {
    timeisset <- FALSE
    run <- vector()
    ## get all the enzymes files in that experiment
    pattern = paste(pm,"-[0-9]+.csv",sep="")
    files <- list.files(path = exp, pattern = pattern, recursive = TRUE)

    # for each trial file (run)
    whole <- vector("list") # list mode vector, list of matrices
    mindfsize <- Inf
    trial <- 1
    for (f in files) {

      odata <- read.csv(file = paste(exp, f, sep="/"), colClasses = "numeric")
      dims <- dim(odata)
      ##print(paste("str odata = ",str(odata)))
      ## time column
      if (timeisset == FALSE) {
          run <- odata[1]
          timeisset <- TRUE
      }

      ## compensate for component ordering (e.g. Body,Harness vs. Harness,Body) with a moving average over PREWINDOW
      PREWINDOW <- ma.check(odata[,2:ncol(odata)], PREWINDOW)
      odata.prewindow <- apply(odata[,2:ncol(odata)], 2, ma.cent, n=PREWINDOW)
      whole[[trial]] <- odata.prewindow
      mindfsize <- min(mindfsize,dim(whole[[trial]])[1])
      trial <- trial + 1

    }

    ## truncate outsized data frames
    trial <- 1
    for (i in 2:length(whole)-1) {
        if (dim(whole[[trial]])[1] > mindfsize) {
            whole[[trial]] <- whole[[trial]][1:mindfsize,]
            run <- run[1:mindfsize,]
        }
    }

    ##print(paste("length of whole = ",length(whole)))
    if (length(whole) == 1) {
      expmean <- whole[[1]]
    } else {
      ## get the averages for each node, for each time, over all trials
      trialsum <- whole[[1]]
      for (i in 2:length(whole)-1) {

        ##print(paste("adding in", whole[[i+1]]))
        try(trialsum+whole[[i+1]])
        if (geterrmessage() != "") q()
        trialsum <- trialsum + whole[[i+1]]
      }
      expmean <- trialsum/length(whole)
    }

    ## table of outflow per node per time
    run <- cbind(run, expmean)
    result <- run
    ## set the column names
    colnames(result) <- colnames(odata)[1:length(colnames(odata))]

    if (pm == "rawInput") {
      outFile <- "avgInput"
      ins <- result
    } else {
      outFile <- "avgOutput"
      outs <- result
    }

    ## write the data to a file
    write.csv(x=result, file=paste(exp, "_", outFile, ".csv", sep=""), row.names=F)
  }

  of <- outFract(ins,outs)
  er <- extRatio(ins,outs)
  write.csv(of,file=paste(exp,"_","outFract",".csv",sep=""),row.names=F)
  write.csv(er,file=paste(exp,"_","extRatio",".csv",sep=""),row.names=F)

}

#q()
