#! /usr/bin/Rscript

###
## Read multiple *.csv files and calculate exposure for each column.
##
## Time-stamp: <2019-07-11 09:44:48 gepr>
###

argv <- commandArgs(T)

require(stats) # for statistics
tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

usage <- function() {
  print("Usage: calc-exposure.r <dCV|dPV> <dMin> <dMax> <exp directories> ...")
  print("e.g. calc-exposure.r dCV 0 100 Mouse")
  quit()
}
if (length(argv) < 4) usage()

calc.derivative <- function(dat.time, dat) {
	dat.dxdt <- diff(as.matrix(dat[,2:ncol(dat)]))
	dat.tmp <- as.data.frame(dat.dxdt)
	dat.tmp <- cbind(dat.time[-length(dat.time)],dat.dxdt)
	colnames(dat.tmp) <- colnames(dat)
	return(dat.tmp)
}

read.data.file <- function(x, xdir, measure, drange) {
	if (any(c(grepl("body-avg", measure),grepl("extra-avg", measure)))) {
		suffix <- paste("_",measure,".csv",sep="")
	} else {
		suffix <- paste("_",measure,"-",drange,".csv",sep="")
	}
	infile <- paste(xdir,basename(x),suffix,sep="")
	if (!file.exists(infile)) {
		print(paste(infile,"doesn't exist for clearance. exiting."))
		quit()
	} else {
		print(paste("Reading file",infile))
		dat <- read.csv(infile)
		#dat <- read.csv(file = filename, colClasses = "numeric")
	}
	return(dat)
}

direction <- argv[1]
dMin <- as.numeric(argv[2])
dMax <- as.numeric(argv[3])
exps <- argv[-(1:3)] ## all remaining args
nxs <- length(exps)
drange <- paste(direction,"∈[",dMin,",",dMax,")",sep="")

measures <- c("entries","celladj")

## read in data from files in <exp>-reduced directory
for (x in exps) {
	x <- paste(dirname(x), basename(x), sep="/") 
	xdir <- paste(x,"-reduced/",sep="")
	if (!file.exists(xdir)) {
		print(paste("directory",xdir,"doesn't exist. skipping."))
		next
	} else {
		print(paste("Working on exposure for experiment", xdir))
	}	
	
	## read amount and other statistics of vHPCs 
    infile <- paste(xdir,basename(x),"_Hcounts-all-",direction,".csv",sep="")
	if (!file.exists(infile)) {
		print(paste(infile,"doesn't exist. skipping."))
		next
	} else {
		print(paste("Reading file",infile))
	}
	# read all available data: each trial value, total, avg, stddev
	allHdat <- read.csv(infile, row.names=1, check.names=F)
	
	# keep the trial vs distance matrix with total columns appended
	Hdat.avg <- allHdat[1:(nrow(allHdat)-1),1:(ncol(allHdat)-2)]
	Hdat.total <- allHdat[1:(nrow(allHdat)-2),1:(ncol(allHdat)-2)]
	
	colnms <- colnames(Hdat.total)
	colnms <- colnms[1:length(colnms)-1]
	distances <- as.numeric(colnms)
	ndpts <- length(distances)
	totalH <- Hdat.total["total","total"]
	dtotalH <- Hdat.total["total", ]
	davgH <- Hdat.avg["avg", ]
	dtotalH <- dtotalH[ ,1:(ncol(dtotalH)-1)]
	davgH <- davgH[ ,1:(ncol(davgH)-1)]
	rownames(dtotalH) <- c()
	rownames(davgH) <- c()
	ntrials <- nrow(Hdat.total)-1
	
	maxdist <- max(distances)
	if (dMin > maxdist) {
		print(paste("dMin",dMin,"is > max distance",maxdist,"for direction",direction,"in experiment",x))
		next
	}
	
	totalH.drange <- sum(dtotalH[,dMin <= distances & distances < dMax])
	avgH.drange <- sum(davgH[,dMin <= distances & distances < dMax])
	temp <- totalH.drange/ntrials
	
	## read cumulative entries
	entries <- read.data.file(x, xdir, "entries", drange)
	celladj <- read.data.file(x, xdir, "celladj", drange)
		
	## calculate derivative to get rate
	entries.time <- entries[,1]
	entries.dxdt <- calc.derivative(entries.time, entries)
		
	if (any(is.na(entries.dxdt))) {
		entries.dxdt[is.na(entries.dxdt)] <- 0 # replace NAs with zeros
	} 
	
	# calculate exposure per vHPC
	# the average number of vHPCs is a sum of averages within a distance range
	exposure.pH <- cbind(entries.dxdt[,1],entries.dxdt[,2:ncol(entries.dxdt)]/avgH.drange)
	colnames(exposure.pH) <- colnames(entries.dxdt)
	
	write.csv(entries.dxdt,paste(xdir,basename(x),"_exposure-entries-",drange,".csv",sep=""),row.names=F)
	write.csv(exposure.pH,paste(xdir,basename(x),"_exposure-entries-perH-",drange,".csv",sep=""),row.names=F)
}

quit()
