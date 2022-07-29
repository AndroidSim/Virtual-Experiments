#! /usr/bin/Rscript

###
## Read multiple *.csv files and calculate derivative wrt time for each column.
##
## Time-stamp: <2019-07-11 09:45:58 gepr>
###

argv <- commandArgs(T)

ma.window <- 181

require(stats) # for statistics
tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

usage <- function() {
  print("Usage: derivative.r <CSV file1> <CSV file2> ...")
  print("e.g. derivative.r exp_celladj-dCV-avg-pHPC-pMC∈[0,8).csv")
  quit()
}
if (length(argv) < 1) usage()

datafiles <- argv

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
	dat <- read.csv(f, check.names=FALSE, colClasses="numeric")
	dat.time <- dat[,1]
	maxtime <- max(dat.time)
	
	## take the derivative
	dat.dxdt <- diff(as.matrix(dat[,2:ncol(dat)]))
	dat.dxdt <- cbind(dat.time[-length(dat.time)],dat.dxdt)
	# add row at end to get original number of time data points
	#dat.dxdt <- rbind(dat.dxdt,c(maxtime,matrix(0, nrow=1, ncol=ncol(dat.dxdt)-1)))
	dat.dxdt <- rbind(dat.dxdt,c(maxtime,dat.dxdt[nrow(dat.dxdt),2:ncol(dat.dxdt)]))
	colnames(dat.dxdt) <- colnames(dat)

	# apply moving average
	dat.dxdt[is.na(dat.dxdt)] <- 0 # replace NAs with zeros
	dat.dxdt <- as.data.frame(dat.dxdt)
    ma.window <- ma.check(dat.dxdt[,2:ncol(dat.dxdt)], ma.window)
	dat.ma <- apply(dat.dxdt[,2:ncol(dat.dxdt)], 2, ma.cent, n=ma.window)
	dat.ma <- cbind(dat.dxdt[,1], dat.ma)
	dat.ma <- as.data.frame(dat.ma)
	colnames(dat.ma) <- colnames(dat)

	write.csv(dat.dxdt,paste(fileName.base,"-dxdt.csv",sep=""),row.names=F)
	write.csv(dat.ma,paste(fileName.base,"-dxdt-ma.csv",sep=""),row.names=F)
}
quit()
