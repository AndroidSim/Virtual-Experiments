#! /usr/bin/Rscript

###
## Read multiple *.csv files and calculate clearance for each column.
##
## Time-stamp: <2019-11-11 19:12:52 gepr>
###

argv <- commandArgs(T)

require(stats) # for statistics
tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

usage <- function() {
  print("Usage: calc-clearance.r <dPV|dCV> dMin dMax <exp directories> ...")
  print("       dMin = min dCV or dPV distance")
  print("       dMax = max dCV or dPV distance")
  print("e.g. calc-clearance.r dPV 0 100 data-time-format")
  print("  directories should contain files like: ")
  print("  <exp>_body-avg.csv and <exp>_entries-dPV∈[0,100).csv")
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

## read in data from files in <exp>-reduced directory
for (x in exps) {
	x <- paste(dirname(x), basename(x), sep="/") 
	xdir <- paste(x,"-reduced/",sep="")
	if (!file.exists(xdir)) {
		print(paste("directory",xdir,"doesn't exist. skipping."))
		next
	} else {
		print(paste("Working on clearance for experiment", xdir))
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
	
	## read cumulative entries and exits
	entries <- read.data.file(x, xdir, "entries", drange)
	exits <- read.data.file(x, xdir, "exits", drange)
	
    ## read amounts in Body
    bodyavg <- read.data.file(x, xdir, "body-avg", drange)
    bodyavg.data <- bodyavg[,2:ncol(bodyavg)]
    
    ## clearance = d(entries - exits)/dt 
    inminusout <- (entries[,2:ncol(entries)]-exits[,2:ncol(exits)])
    inminusout <- cbind(entries[,1],inminusout)
    colnames(inminusout) <- colnames(entries)
    ## calculate derivative to get rate
    inminusout.time <- inminusout[,1]
    inminusout.dxdt <- calc.derivative(inminusout.time, inminusout)
    
    # calculate clearance per vHPC
	# the average number of vHPCs is a sum of averages within a distance range
	inminusout.pH <- cbind(inminusout.dxdt[,1],inminusout.dxdt[,2:ncol(inminusout.dxdt)]/avgH.drange)
	colnames(inminusout.pH) <- colnames(inminusout.dxdt)
    
    ## Calculate clearance (CLint) = dose/AUC using linear trapezoidal method
    ## code from ishc:
    #normy <- dat[,compound]/dat[1,compound]    	
    #i <- 2:new_final_time
    #AUC <- as.double( (x[i] - x[i-1]) %*% (normy[i] + normy[i-1])) / 2
    #doseToAUC[trial] <- (1-normy[new_final_time])/(D*AUC)
    
    ## specify what data to use for AUC calculation
    #AUC.data <- totalamt.data
    #AUC.time <- totalamt.time
    AUC.data <- bodyavg.data
    AUC.time <- bodyavg[,1]
    
    ## calculate AUC and dose/AUC
    columns <- colnames(AUC.data)
    AUC <- as.data.frame(matrix(data = NA, nrow = 1, ncol = length(columns)))
    colnames(AUC) <- columns
    doseToAUC <- as.data.frame(matrix(data = NA, nrow = 1, ncol = length(columns)))
    colnames(doseToAUC) <- columns
    rate.removed <- as.data.frame(matrix(data = NA, nrow = 1, ncol = length(columns)))
    colnames(rate.removed) <- columns
    for (cn in columns) {
		maxdata <- max(AUC.data[,cn])
		istart <- which.max(AUC.data[,cn])
		maxtime <- AUC.time[istart]
		iend <- length(AUC.time)
		i <- (istart+1):iend
		#i <- (maxtime+1):AUC.time[length(AUC.time)]
		AUC[,cn] <- ((AUC.data[i,cn] + AUC.data[i-1,cn])/2) %*% (AUC.time[i] - AUC.time[i-1])
		doseToAUC[,cn] <- maxdata/AUC[,cn]
		rate.removed[,cn] <- sum(abs(diff(AUC.data[i,cn])))/(AUC.time[iend] - maxtime)
    }
    output <- rbind(AUC,doseToAUC)
    output <- rbind(output,rate.removed)
    rownames(output) <- c("AUC","doseToAUC","rate removed")
    
    ## write output files
	write.csv(inminusout.dxdt,paste(xdir,basename(x),"_clearance-[entries-exits]","-",drange,".csv",sep=""),row.names=F)
	write.csv(inminusout.pH,paste(xdir,basename(x),"_clearance-[entries-exits]-perH","-",drange,".csv",sep=""),row.names=F)
	write.csv(output,paste(xdir,basename(x),"_clearance-AUC","-",drange,".csv",sep=""),row.names=T)
}

quit()
