#! /usr/bin/Rscript

###
## Plot ratio of two data sets vs Time.
##
## Time-stamp: <2019-07-11 09:47:46 gepr>
###

argv <- commandArgs(TRUE)

ma.window <- 181

require(stats) # for statistics
tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

usage <- function() {
  print("Usage: ratio.r <numerator CSV file1> <denominator CSV file2> ...")
  print(" 			  the arguments are in numerator denominator pairs")
  print("e.g. ratio.r exp_celladj-dCV-avg-pHPC-pMC∈[0,8).csv... ")
  print("				...exp_mobileObject-dPV-avg-pHPC-pMC∈[5,10).csv")
  quit()
}
if (length(argv) < 2) usage()

is.infinite.data.frame <- function(obj) {
  sapply(obj,FUN = function(x) is.infinite(x))
}

parse.file.name <- function(f) {
  ## parse file name
  seps <- gregexpr('/',f)[[1]] # get all the '/' locations
  aftersep <- substr(f,seps[length(seps)]+1,nchar(f)) # get everything after the last '/'
  expname <- substr(aftersep,0,regexpr('_',aftersep)-1)
  compname <- substr(f,regexpr('_',f)+1,nchar(f))
  fileName.base <- paste(expname,substr(compname, 0, regexpr('.csv', compname)-1),sep='_')
  return(fileName.base)
}

calc.ratio <- function(df1,df2) {
  ## calculate ratio of rate over concentration, df1/df2
  numerator <- df1
  denominator <- df2
  ratio <- numerator/denominator
  ##ratio <- numerator[2:ncol(numerator)]/denominator[2:ncol(denominator)]
  ##ratio <- as.data.frame(cbind(numerator[,1],ratio))
  ##colnames(ratio) <- c("Time",colnames(ratio)[2:ncol(numerator)])
  ratio[is.na(ratio)] <- 0
  ratio[is.infinite.data.frame(ratio)] <- NA
  return(ratio)
}

nArgs <- length(argv)
if (nArgs%%2 != 0) {
	print("the number or arguments must be even as numerator denominator pairs.")
    q("no")
 }

##exps <- argv[-(1:2)] ## all remaining args
datafiles <- argv
##for (f in datafiles) {
##while (!is.empty(datafiles)) {
while (length(datafiles) != 0) {
  numf <- datafiles[1]
  denf <- datafiles[2]
  print(paste("Working on numerator", numf,"and denominator", denf))

  if (!file.exists(numf) || !file.exists(denf)) {
    if (!file.exists(numf)) {
      print(paste("numerator",numf,"doesn't exist."))
    }
    if (!file.exists(denf)) {
      print(paste("denominator",denf,"doesn't exist."))
    }
    print(paste("skipping", numf,"to", denf,"ratio pair"))
    datafiles <- datafiles[-(1:2)]
    next
  }

  ## parse numerator and denominator file names
  numfb <- parse.file.name(numf)
  denfb <- parse.file.name(denf)

  ## read data and process
  numdat <- read.csv(numf)
  dendat <- read.csv(denf)
  numdat.time <- numdat[,1]
  dendat.time <- dendat[,1]
  if (length(numdat.time) != length(dendat.time)) {
    print(paste("number of time points don't match for",numf,"and",denf))
    print(paste("skipping", numf,"to", denf,"ratio pair"))
    datafiles <- datafiles[-(1:2)]
    next
  }
  dat.time <- numdat.time

  ## find common column names between numerator and denominator
  numcoln <- colnames(numdat)[2:ncol(numdat)]
  dencoln <- colnames(dendat)[2:ncol(dendat)]
  vCcoln <- intersect(numcoln,dencoln)

  ## calculate ratio
  for (vCom in vCcoln) {
    tmp <- calc.ratio(numdat[,vCom],dendat[,vCom])
    if (exists("ratio")) {
      ratio <- cbind(ratio,tmp)
    } else {
      ratio <- tmp
    }
  }
  ratio <- as.data.frame(cbind(dat.time,ratio))
  colnames(ratio) <- c("Time",vCcoln)

  ## apply moving average
  ratio[is.na(ratio)] <- 0 # replace NAs with zeros
  ma.window <- ma.check(ratio[,2:ncol(ratio)], ma.window)
  ratio.ma <- apply(ratio[,2:ncol(ratio)], 2, ma.cent, n=ma.window)
  ratio.ma <- cbind(dat.time, ratio.ma)
  ratio.ma <- as.data.frame(ratio.ma)
  colnames(ratio.ma) <- colnames(ratio)

  write.csv(ratio,paste(numfb,"-to-",denfb,"-ratio.csv",sep=""),row.names=F)
  write.csv(ratio.ma,paste(numfb,"-to-",denfb,"-ratio-ma.csv",sep=""),row.names=F)

  remove("ratio")
  datafiles <- datafiles[-(1:2)]
}
quit()
