#! /usr/bin/Rscript

###
##
## Use rawInput, traps/intracellular amounts, and celladj to 
## calculate the following:
##
##   Output Fraction: output/input
##   Extraction Ratio: (input - output)/input
##
## for Layer 1 (next to PV) of the Lobule graph.
##
## output(Layer) = input(Layer) - intra(band) - celladj(band)
##
###

argv <- commandArgs(T)

require(stats) # for statistics
tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

usage <- function() {
  print("Usage: layer-outext-ratios.r <dPV|dCV> dMin dMax <exp directories> ...")
  print("       dMin = min dCV or dPV distance")
  print("       dMax = max dCV or dPV distance")
  print("e.g. layer-outext-ratios.r dCV 0 100 data-time-format")
  print("  directories should contain files like <exp>_avgInput.csv,")
  print("  <exp>_extra-avg.csv, and <exp>_celladj-d[CP]V-[0-9]+.csv")
  quit()
}
if (length(argv) < 4) usage()

read.data.file <- function(x, xdir, measure, drange) {
	if (measure %in% c("body-avg","extra-avg","avgInput")) {
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

direction <- argv[1]
dMin <- as.numeric(argv[2])
dMax <- as.numeric(argv[3])
exps <- argv[-(1:3)] ## all remaining args
nxs <- length(exps)
drange <- paste(direction,"∈[",dMin,",",dMax,")",sep="")

filebases <- c("celladj", "mobileObject", "entries", "exits", "rejects", "traps")

## read in data from files in <exp>-reduced directory
## average before calculation
for (x in exps) {
	x <- paste(dirname(x), basename(x), sep="/") 
	xdir <- paste(x,"-reduced/",sep="")
	if (!file.exists(xdir)) {
		print(paste("directory",xdir,"doesn't exist. skipping."))
		next
	} else {
		print(paste("Working on exposure for experiment", xdir))
	}
	
	## read input into Layer 0
	input <- read.data.file(x, xdir, "avgInput", drange)
	input[is.na(input)] <- 0 # replace NAs with zeros
	## read cumulative traps
	traps <- read.data.file(x, xdir, "traps", drange)
    ## read amounts in extracellular and intracellular spaces
    celladj <- read.data.file(x, xdir, "celladj", drange)
    intra <- read.data.file(x, xdir, "mobileObject", drange)
    extra <- read.data.file(x, xdir, "extra-avg", drange)
    
    ## calculate output
    input.data <- input[,2:ncol(input)]
    intra.data <- intra[,2:ncol(intra)]
    extra.data <- extra[,2:ncol(extra)]
    traps.data <- traps[,2:ncol(traps)]
    celladj.data <- celladj[,2:ncol(celladj)]
    
    cninput <- colnames(input.data)
    cncelladj <- colnames(celladj.data)
    cnintra <- colnames(intra.data)
    cnextra <- colnames(extra.data)
    
    cncommon <- intersect(intersect(intersect(cninput,cncelladj),cnintra),cnextra)
    # or Reduce(intersect, list(a,b,c)) or unique(c[c%in%a[a%in%b]])
    totalamt <- as.data.frame(matrix(data = NA, nrow = length(input[,1]), ncol = length(cncommon)))
    output <- as.data.frame(matrix(data = NA, nrow = length(input[,1]), ncol = length(cncommon)))
    total.input <- as.data.frame(matrix(data = NA, nrow = length(input[,1]), ncol = length(cncommon)))
    layer.outfract <- as.data.frame(matrix(data = NA, nrow = length(input[,1]), ncol = length(cncommon)))
    layer.extratio <- as.data.frame(matrix(data = NA, nrow = length(input[,1]), ncol = length(cncommon)))
    colnames(totalamt) <- cncommon
    colnames(output) <- cncommon
    colnames(total.input) <- cncommon
    colnames(layer.outfract) <- cncommon
    colnames(layer.extratio) <- cncommon
    for (cn in cncommon) {
		totalamt[,cn] <- extra.data[,cn]+celladj.data[,cn]+intra.data[,cn]
		total.input[,cn] <- cumsum(input.data[,cn])
		output[,cn] <- total.input[,cn]-totalamt[,cn]
		#layer.extratio[,cn] <- extRatio(total.input[,cn],output[,cn]) 
		layer.outfract[,cn] <- (total.input[,cn]-totalamt[,cn])/total.input[,cn]
		layer.extratio[,cn] <- 1-layer.outfract[,cn]
    }
    
    #remaining <- totalamt[,"APAP"]+totalamt[,"G"]+totalamt[,"S"]
    #LayerOutFract <- (total.input[,"APAP"]-remaining)/total.input[,"APAP"]
    #LayerOutFract <- cbind(input[,1],LayerOutFract)
    #colnames(LayerOutFract) <- c("Time","outFract-APAP")
    #LayerOutFract[is.infinite(LayerOutFract)] <- NA
    
	totalamt <- cbind(input[,1],totalamt)
	colnames(totalamt) <- c("Time",cncommon)
	
	output <- cbind(input[,1],output)
	colnames(output) <- c("Time",cncommon)
	
	total.input <- cbind(input[,1],total.input)
	colnames(total.input) <- c("Time",cncommon)
	
	layer.outfract <- cbind(input[,1],layer.outfract)
	colnames(layer.outfract) <- c("Time",cncommon)
	layer.outfract[is.infinite(layer.outfract)] <- NA
	layer.outfract[layer.outfract < 0] <- NA
	
	layer.extratio <- cbind(input[,1],layer.extratio)
	colnames(layer.extratio) <- c("Time",cncommon)
	layer.extratio[is.infinite(layer.extratio)] <- NA
    
    ## write output files
	#write.csv(total.input,paste(xdir,basename(x),"_total-input","-",drange,".csv",sep=""),row.names=F)
	#write.csv(totalamt,paste(xdir,basename(x),"_totalamt","-",drange,".csv",sep=""),row.names=F)
	#write.csv(output,paste(xdir,basename(x),"_output","-",drange,".csv",sep=""),row.names=F)
	write.csv(layer.outfract,paste(xdir,basename(x),"_layer-outFract","-",drange,".csv",sep=""),row.names=F)
	write.csv(layer.extratio,paste(xdir,basename(x),"_layer-extRatio","-",drange,".csv",sep=""),row.names=F)
	#write.csv(LayerOutFract,paste(basename(x),"_layer-outFract-APAP","-",drange,".csv",sep=""),row.names=F)
}

# calculation before average
calcb4avg <- F
if (calcb4avg) {
	for (d in exps) {
	  d <- paste(dirname(d), basename(d), sep="/") # removes trailing slash from the argument
	  print(paste("Working on", d))

	  if (!file.exists(d)) {
		print(paste(d,"doesn't exist."))
		next
	  }

	  ## for both measures
	  for (measure in filebases) {

		pattern <- paste(measure,"-",direction,"-[0-9]+.csv.gz",sep="")
		allfiles <- list.files(path=d, pattern=pattern, recursive=T,  full.names=T)
		if (length(allfiles) <= 0) {
		  print(paste("skipping ", measure, "-", direction, sep=""))
		  next
		}
		totals <- avgByColumn(allfiles)
		fileprefix <- paste(d,"_",measure,"-",direction,sep="")
		write.csv(totals,
				  file=paste(fileprefix,".csv",sep=""),
				  row.names=F)
		distances <- unlist(strsplit(colnames(totals),":"))
		distances <- unique(distances[c(F,T)])
		Maxdist <- max(as.numeric(distances))
		if (dMin > Maxdist) {
		  print(paste("dMin, ",dMin,", is > max distance, ",Maxdist,", for direction, ",direction))
		  next
		}
		inband <- sumBandByLastTag(totals,c(dMin,dMax))
		write.csv(inband,
				  file=paste(fileprefix,"∈[",dMin,",",dMax,").csv",sep=""),
				  row.names=F)

	  } ## end   for (measure in filebases) {
	} ## end for (d in exps) {
}

quit()
