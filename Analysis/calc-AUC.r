#! /usr/bin/Rscript

###
## Read multiple *.csv files and calculate the area-under-the-curve (AUC) for each column.
##
## Time-stamp: <2020-12-07 aks>
###

argv <- commandArgs(T)

require(stats) # for statistics
tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

usage <- function() {
  print("Usage: calc-AUC.r start.time end.time <CSV file1> <CSV file2> ...")
  print("e.g. ./bin/calc-AUC.r 100 1000 exp_celladj-dCV-avg-pHPC-pMC∈[0,8).csv")
  quit()
}
if (length(argv) < 3) usage()

start.time <- as.numeric(argv[1])
end.time <- as.numeric(argv[2])
datafiles <- argv[-(1:2)] ## all remaining args
nfiles <- length(datafiles)

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

	## read data
	dat <- read.csv(f, check.names=FALSE, colClasses="numeric")
	dat.time <- dat[,1]
	
	## calculate AUC by quadrature
	columns <- colnames(dat[,2:ncol(dat)])
    AUC <- as.data.frame(matrix(data = NA, nrow = 1, ncol = length(columns)))
    colnames(AUC) <- columns
    istart <- which(dat.time == start.time)
	iend <- which(dat.time == end.time)
	i <- (istart+1):iend
    for (cn in columns) {
		AUC[,cn] <- ((dat[i,cn] + dat[i-1,cn])/2) %*% (dat.time[i] - dat.time[i-1])
    }	
    print(AUC)
	write.csv(AUC,paste(fileName.base,"-AUC.csv",sep=""),row.names=F)
}

quit()
