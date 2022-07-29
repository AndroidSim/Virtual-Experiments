#! /usr/bin/Rscript

###
## Read multiple *.csv files and calculate bulk avg gradient
##
## Time-stamp: <2020-03-31 aks>
###

argv <- commandArgs(T)

PLOT.SVG <- F

require(stats) # for statistics
tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })
#library(sigmoid)

usage <- function() {
	print("Usage: calc-bulk-gradient.r <gradient shape> <dPV value> <dCV value> <exp directories> ...")
	print("e.g. calc-bulk-gradient.r sigmoid 0 1 data-time-format")
	quit()
}
if (length(argv) < 4) usage()

sigmoid <- function(x) {
   1 / (1 + exp(-x))
}

##
## test for and create graphics subdirectory
##
if (!file.exists("graphics")) dir.create("graphics")

gshape <- argv[1]
vdPV <- as.numeric(argv[2])
vdCV <- as.numeric(argv[3]) 
exps <- argv[-(1:3)] ## all remaining args
nxs <- length(exps)
directions <- c("dPV", "dCV")

## read in data from files in <exp>-reduced directory
## code is from plot-hcounts.r in github repository "scripts"
report <- numeric()
xnames <- numeric()
for (x in exps) {
	x <- paste(dirname(x), basename(x), sep="/") 
	xdir <- paste(x,"-reduced/",sep="")
	if (!file.exists(xdir)) {
		print(paste("directory",xdir,"doesn't exist. skipping."))
		next
	} else {
		print(paste("Working on experiment", xdir))
	}
	bagd <- numeric()
	for (direction in directions) {
		infile <- paste(xdir,basename(x),"_Hcounts-all-",direction,".csv",sep="")
		if (!file.exists(infile)) {
			print(paste(infile,"doesn't exist. skipping."))
			next
		} else {
			print(paste("Reading file",infile))
			dat <- read.csv(infile)
		}
		
		# read all available data: each trial value, total, avg, stddev
		allHdat <- read.csv(infile, row.names=1, check.names=F)
		# keep the trial vs distance matrix with total columns appended
		allHdat <- allHdat[1:(nrow(allHdat)-2),1:(ncol(allHdat)-2)]
		colnms <- colnames(allHdat)
		colnms <- colnms[1:length(colnms)-1]
		distances <- as.numeric(colnms)
		#distances <- as.numeric(colnames(allHdat[,1:(ncol(alldatH)-1)]))
		ndpts <- length(distances)
		
		totalH <- allHdat["total","total"]
		dtotals <- allHdat["total", ]
		dtotals <- dtotals[ ,1:(ncol(dtotals)-1)]
		rownames(dtotals) <- c()
		
		## calculate prob(vHPC) or fraction(vHPC) vs distance as 
		## prob = num(vHPCs) at that distance / total num of vHPCs
		## or, in above variables, dtotals/totalH
		fractH <- dtotals/totalH
		
		## interpolate gradient from end point values and shape property
		## relevant variables from input arguments = vdPV, vdCV, gshape
		## gshape = linear, sigmoid, or exponential
		
		## The following is from function setGradients in Hepatocyte.java in islj:
		## double rps = (Double)eg.getProperty("rxnProbStart");
		## double rpf = (Double)eg.getProperty("rxnProbFinish");
		## double rp = Double.NaN;
		## String rpg = (String)eg.getProperty("rxnProbGradient");
		## if (rpg == null || rpg.contains("linear")) {
		##   rp = LinearGradient.eval(rps,rpf,0.0,(double)(dPV+dCV),(double)dPV);
		## } else if (rpg.contains("sigmoid")) {
		##   rp = SigmoidGradient.eval(rps,rpf,0.0,(double)(dPV+dCV),(double)dPV);
		## }
		
		## The following is from LinearGradient.java in bsg.util:
		## public static double eval(double refX, double refY, double valX, double valY, double x) {
		## 		double retVal = refX + x * (refY-refX)/(valY-valX);
		## linear: vG; = vdPV + distPV * (vdCV-vdPV)/((distPV+distCV)-0)
		
		## The following is from SigmoidGradient.java in bsg.util:
		## public static double eval(double refX, double refY, double intensity, double valX, double valY, double x) {
		##		double val = LinearGradient.eval(10, 0, valX, valY, x);
		##		double retVal = refX + (refY-refX)/(1.0+StrictMath.exp(intensity*(val-5)));
		
		## The following is from ExpGradient.java in bsg.util:
		## public static double eval(double refX, double refY, double intensity, double valX, double valY, double x) {
		##		double valMin = Math.min(valX,valY), valMax = Math.max(valX,valY);
		##		// if it is inverted, invert the vals 
		##		double unitPos = (refX < refY ? (x-valMin)/(valMax-valMin) : 1-((x-valMin)/(valMax-valMin)));
		##		retVal = unitPos*((StrictMath.exp(intensity*unitPos)))/(Math.pow(Math.E,intensity));
		##		// scale it
		##		double refMin = Math.min(refX,refY), refMax = Math.max(refX,refY);
		##		retVal = refMin + retVal*(refMax-refMin);    
		
		maxd <- max(distances)
		mind <- min(distances)
		if (grepl("linear", gshape)) {
			if (grepl("dPV", direction)) {
				valG <- vdPV + (distances/maxd)*(vdCV - vdPV)
			} else {
				valG <- vdCV + (distances/maxd)*(vdPV - vdCV)
			}
		} else if (grepl("sigmoid", gshape)) {
			if (grepl("dPV", direction)) {
				val <- 10 + (distances/maxd)*(0 - 10)
				valG <- vdPV + (vdCV-vdPV)/(1.0 + exp(val - 5)) 
			} else {
				#val <- 0 + (distances/maxd)*(10 - 0)
				val <- 10 + (distances/maxd)*(0 - 10)
				valG <- vdCV + (vdPV-vdCV)/(1.0 + exp(val - 5))
			}
		} else if (grepl("exponential", gshape)) {
			#E <- exp(1)
			if (vdPV < vdCV) {
				unitPos <- (distances - mind)/(maxd - mind)
			} else {
				unitPos <- 1 - ((distances - mind)/(maxd - mind))
			}
			val <- unitPos*exp(unitPos)
			minr <- min(c(vdPV,vdCV))
			maxr <- max(c(vdPV,vdCV))
			valG <- minr + val*(maxr - minr)
		} else {
		}
		
		## calculate bulk average gradient, avg = sum gradient * prob(vHPC)
		bulkavgG <- sum(valG*fractH)
		
		## combine bulk average gradient for each direction and experiment
		## into a report to be printed
		bagd <- c(bagd,bulkavgG)
		
		## plot gradient and probability of vHPC
		outfile <- paste(basename(x),"_","bulk-avg-gradient-",direction,sep="")
		outFile <- paste("graphics/",outfile, sep="")
		if (PLOT.SVG) {
			svg(paste(outFile,".svg",sep=""),10,10)
		} else {
			png(paste(outFile, ".png", sep=""), 1600, 1600)
		}
		par(mfrow=c(2,1), mar=c(5,6,4,2), cex.main=2, cex.axis=2, cex.lab=2)
        plot(distances, valG,
             main=paste(basename(x),"-Gradient-",direction,sep=""),
             xlab=direction, ylab="Gradient",
             xlim=c(mind,maxd), ylim=c(min(valG),max(valG)),
             type="l", lwd=5,
             col=2
             )
        plot(distances, fractH,
             main=paste(basename(x),"-prob(vHPC)-",direction,sep=""),
             xlab=direction, ylab="prob(vHPC)",
             xlim=c(mind,maxd), ylim=c(min(fractH),max(fractH)),
             type="l", lwd=5,
             col=3
             )
	} # end loop over directions
	## create report for experiment
	report <- rbind(report,bagd)
	xnames <- c(xnames,basename(x))
} # end loop over experiments
if (nxs > 1) {
	report <- as.matrix(report)
} else {
	report <- t(as.matrix(report))
}
colnames(report) <- c("bulk-avg-gradient-dPV","bulk-avg-gradient-dCV")
rownames(report) <- xnames
print(report)
quit()
