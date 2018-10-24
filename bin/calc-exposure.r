#! /usr/bin/Rscript

###
## Read multiple *.csv files and calculate exposure for each column.
##
###

argv <- commandArgs(T)

require(stats) # for statistics
source("~/R/misc.r")

usage <- function() {
  print("Usage: calc-exposure.r dTime Hcount <CSV file1> <CSV file2> ...")
  print("       dTime = time interval to sample data for exposure")
  print("               default = 50")
  print("       Hcount = number of aHPCs to scale data")
  print("               default = 1000")
  print("e.g. calc-exposure.r 50 1000 exp_celladj-dCV-avg-pHPC-pMC∈[0,8).csv")
  quit()
}

if (length(argv) < 3) usage()

dTime <- as.numeric(argv[1])
Hcount <- as.numeric(argv[2])
datafiles <- argv[-(1:2)] ## all remaining args

if (!file.exists("graphics")) dir.create("graphics")

edata <- vector("list")

for (f in datafiles) {
   print(paste("Working on", f))

   if (!file.exists(f)) {
      print(paste(f,"doesn't exist."))
      next
   }

   ## parse file name
   seps <- gregexpr('/',f)[[1]] # get all the '/' locations
   aftersep <- substr(f,seps[length(seps)]+1,nchar(f)) # get everything after the last '/'
   expname <- substr(aftersep,0,regexpr('_',aftersep)-1)
   compname <- substr(f,regexpr('_',f)+1,nchar(f))
   fileName.base <- paste(expname,substr(compname, 0, regexpr('.csv', compname)-1),sep='_')

   ## read data and process
   dat <- read.csv(f)
   dat.time <- dat[,1]
   dat.tmp <- dat

   if (grepl("entries", fileName.base)) {
     ## take the derivative
     dat.dxdt <- diff(as.matrix(dat[,2:ncol(dat)]))
     dat.tmp <- as.data.frame(dat.dxdt)
     dat.tmp <- cbind(dat.time[-length(dat.time)],dat.dxdt)
     colnames(dat.tmp) <- colnames(dat)
   }

   dat.tmp[,2:ncol(dat.tmp)] <- Hcount*dat.tmp[,2:ncol(dat.tmp)]
   dat.tmp[is.na(dat.tmp)] <- 0 # replace NAs with zeros
   dat.dxdt <- as.data.frame(dat.tmp)
   dat.ma <- apply(dat.dxdt[,2:ncol(dat.dxdt)], 2, ma.cent, n=dTime)
   dat.ma <- cbind(dat.dxdt[,1], dat.ma)
   dat.ma <- as.data.frame(dat.ma)
   colnames(dat.ma) <- colnames(dat)

   
   write.csv(dat.dxdt,paste(fileName.base,"-exposure.csv",sep=""),row.names=F)
   write.csv(dat.ma,paste(fileName.base,"-exposure-ma.csv",sep=""),row.names=F)
}

quit()
