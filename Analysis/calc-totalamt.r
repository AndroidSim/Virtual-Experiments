#! /usr/bin/Rscript

###
## Read multiple *.csv files and calculate total amount of each
## mobile object.
##
## Time-stamp: <2021-03-15 aks>
###

argv <- commandArgs(T)

require(stats) # for statistics
tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

usage <- function() {
  print("Usage: calc-totalamt.r <dPV|dCV> dMin dMax <exp directories> ...")
  print("       dMin = min dCV or dPV distance")
  print("       dMax = max dCV or dPV distance")
  print("e.g. calc-totalamt.r dPV 0 100 data-time-format")
  print("  directories should contain files like: ")
  print("  <exp>_body-avg.csv and <exp>_celladj-dPV∈[0,100).csv")
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
		print(paste(infile,"doesn't exist for calculating total amount. exiting."))
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
		print(paste("Working on total amount of mobile objects for experiment", xdir))
	}
	
    ## read amounts in Body and all Lobular spaces
    bodyavg <- read.data.file(x, xdir, "body-avg", drange)
    celladj <- read.data.file(x, xdir, "celladj", drange)
    extra <- read.data.file(x, xdir, "extra-avg", drange)
    intra <- read.data.file(x, xdir, "mobileObject", drange)
    
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
	
    ## previous: total amount in Body + Lobuler spaces corresponds to Medium in ishc
    ## current: Body in Mouse corresponds to Media in Culture
    bodyavg.data <- bodyavg[,2:ncol(bodyavg)]
    celladj.data <- celladj[,2:ncol(celladj)]
    intra.data <- intra[,2:ncol(intra)]
    extra.data <- extra[,2:ncol(extra)]
    cnbody <- colnames(bodyavg.data)
    cncelladj <- colnames(celladj.data)
    cnintra <- colnames(intra.data)
    cnextra <- colnames(extra.data)
    cncommon <- intersect(intersect(intersect(cnbody,cncelladj),cnintra),cnextra)
    # or Reduce(intersect, list(a,b,c)) or unique(c[c%in%a[a%in%b]])
    totalamt <- as.data.frame(matrix(data = NA, nrow = length(bodyavg[,1]), ncol = length(cncommon)))
    colnames(totalamt) <- cncommon
    for (cn in cncommon) {
		temp <- bodyavg.data[,cn]+extra.data[,cn]+celladj.data[,cn]+intra.data[,cn]
		totalamt[,cn] <- temp
    }
	totalamt <- cbind(bodyavg[,1],totalamt)
	colnames(totalamt) <- c("Time",colnames(totalamt[,2:ncol(totalamt)]))
	
	# calculate total amount per vHPC
	# total amount is actually a sum of averages at different locations
	# within a distance range, and average number of vHPCs is also a sum
	# of averages within a distance range
	totalamt.pH <- cbind(totalamt[,1],totalamt[,2:ncol(totalamt)]/avgH.drange)
	colnames(totalamt.pH) <- colnames(totalamt)
	
	# take derivative
	#totalamt.time <- totalamt[,1]
    #totalamt.data <- totalamt[,2:ncol(totalamt)]
    #totalamt.dxdt <- calc.derivative(totalamt.time, totalamt)
    
    ## write output files
	write.csv(totalamt,paste(xdir,basename(x),"_totalamt","-",drange,".csv",sep=""),row.names=F)
	write.csv(totalamt.pH,paste(xdir,basename(x),"_totalamt-perH","-",drange,".csv",sep=""),row.names=F)
	#write.csv(totalamt.dxdt,paste(xdir,basename(x),"_totalamt-dxdt","-",drange,".csv",sep=""),row.names=F)
}

quit()
