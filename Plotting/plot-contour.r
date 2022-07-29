#! /usr/bin/Rscript

###
## 	Make a contour plot of 3D data structured as a matrix with space (i.e. bands) vs time 
## 	vs measurements (or events). The matrix is inputted from a file as an argument.
##
## 	Time-stamp: <2020-09-11 aks>
##
###

argv <- commandArgs(T)

PLOT.SVG <- F
#PLOT.AVG <- F 
USE.FRAMES <- F
MA.WINDOW <- 181

require(stats) # for statistics
tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

library(ggplot2)
#library(lattice)
#library(gridExtra)
#library(svglite)
#library(plotly)
#library(plot3D)


usage <- function() {
    print("Usage: plot-contour.r <*-xyz.csv file of experiment1> <...experiment2> ...")
    print("  e.g. plot-contour.r Mouse_entries-Compound-dCV-xyz.csv")
    quit()
}
if (length(argv) < 1) usage()

###
## test for and create graphics subdirectory
###
if (!file.exists("graphics")) dir.create("graphics")

'%!in%' <- function(x,y)!('%in%'(x,y))

parse.read.file <- function(fname) {
	temp <- unlist(strsplit(fname,"_"))
	xname <- temp[1]
	#temp2 <- unlist(strsplit(temp[2],"-"))
	#fname.base <- paste(xname,substr(temp[2],0,regexpr('-xyz.csv',temp[2])-1),sep='-')
	mname <- substr(temp[2],0,regexpr('-xyz.csv',temp[2])-1)
	pwdf <- paste(getwd(),fname,sep="/")
	if (file.exists(pwdf)) {
		infile <- pwdf
	} else {
		redf <- paste(xname,"-reduced/",fname,sep="")
		if (file.exists(redf)) {
			infile <- redf
		} else {
			cat("\n")
			print(paste("file",fname,"in directory",paste(xname,"-reduced/",sep=""),"doesn't exist. skipping."))
			return(NA)
		}
	}
	dat <- read.csv(file=infile, check.names=FALSE, colClasses="numeric")
	return(list(xname=xname,mname=mname,dat=dat))
}

Files <- argv
nFiles <- length(Files)
xnames <- list()
mnames <- list()
fnames.base <- list()
all.dat <- list()
pndx <- 1
for (f in Files) {
	out <- parse.read.file(f)
	if (any(is.na(out))) next
	xnames[[pndx]] <- out[[1]]
	mnames[[pndx]] <- out[[2]]
	fnames.base[[pndx]] <- paste(xnames[[pndx]],mnames[[pndx]],sep="-")
	all.dat[[pndx]] <- out[[3]]
	pndx <- pndx + 1
}
nplots <- length(all.dat)
if (USE.FRAMES) {
	## determine # of plot frames
	plot.cols <- round(sqrt(nplots))
	## add a new row if we rounded up
	plot.rows <- ifelse(plot.cols >= sqrt(nplots), plot.cols, plot.cols+1)
	outFile <- paste(unique(unlist(xnames)),collapse="-")
	outFile <- paste("graphics/",outFile,"-contours",sep="")
	#par(mar=c(5,6,4,2), cex.main=2, cex.axis=2, cex.lab=2)
	if (PLOT.SVG) {
		svg(paste(outFile,".svg",sep=""),10,10)
		par(mfrow=c(plot.rows,plot.cols), cex=2, lwd=3, mar=c(5,6,4,2), cex.main=1, cex.axis=1, cex.lab=1)
	} else {
		png(paste(outFile, ".png", sep=""), 1600, 1600)
		par(mfrow=c(plot.rows,plot.cols), lwd=3, mar=c(5,6,4,2), cex.main=2, cex.axis=2, cex.lab=2)
		#par(mfrow=c(plot.rows,plot.cols), lwd=3, mar=c(1,1,1,1), cex.main=2, cex.axis=2, cex.lab=2)
	}
}
pb <- txtProgressBar(min=0,max=nplots,style=3)
setTxtProgressBar(pb,0)
for (p in c(1:nplots)) {
	xname <- xnames[[p]]
	mname <- mnames[[p]]
	fname.base <- fnames.base[[p]]
	dat <- all.dat[[p]]
	
	# x-axis = time, dPV time = dCV time
	x <- dat[,"Time"]
	# y-axis = space, each column is a band
	#y <- 1:(ncol(dat)-1)
	# z-axis = data or measurements
	z <- as.matrix(dat[,2:ncol(dat)])
	bandnames <- colnames(z)
	temp <- unlist(strsplit(bandnames,"∈"))
	direction <- unique(temp[c(T,F)])
	bands <- unique(temp[c(F,T)])
	matches <- regmatches(bands, gregexpr("[[:digit:]]+", bands))
	vbands <- as.numeric(unique(unlist(matches)))
	y <- vbands[-length(vbands)]
	
	## set margins and title, axis, and label font sizes
	##   format is c(bottom, left, top, right)
	if (!USE.FRAMES) {
		## place the right margin 1 unit for each of the characters +  units for the line
		#right.margin <- max(nchar(titles)) + max(nchar(maws)) + ifelse(plot.svg, 0, 12)
		#par(mar=c(5,6,4,right.margin), cex.main=2, cex.axis=2, cex.lab=2)
		outFile <- paste("graphics/", fname.base,"-contour", sep="")
		if (PLOT.SVG) {
			svg(paste(outFile,".svg",sep=""),10,10)
			par(lwd=3, mar=c(5,6,4,2), las=1, cex.main=2, cex.axis=1, cex.lab=1)
		} else {
			png(paste(outFile, ".png", sep=""), 1600, 1600)
			par(lwd=3, mar=c(5,10,4,2), las=2, cex.main=4, cex.axis=3, cex.lab=3)
		}
	}
	
	xlabel <- "time"
	ylabel <- direction
	#ylabels <- vbands
	nlevels <- 10
	maxz <- max(z, na.rm=T)
	minz <- min(z, na.rm=T)
	step <- (maxz-minz)/nlevels
	cols <- rev(rainbow(nlevels))
	cols[1] <- "black"
	if (maxz == 0 && minz == 0) {
		filled.contour(x,y,z, col = cols, main=fname.base, xlab=xlabel, ylab=ylabel)
	} else {
		if (USE.FRAMES) {
			options("max.contour.segments"= 300000)
			contour(x,y,z, col = cols, levels = seq(minz,maxz,step),
				main=fname.base, xlab=xlabel, ylab=ylabel)
		} else {
			filled.contour(x,y,z, col = cols, levels = seq(minz,maxz,step),
				main=fname.base, xlab=xlabel, ylab=ylabel)
		}
	}
	setTxtProgressBar(pb,getTxtProgressBar(pb)+1)
}
close(pb)
quit()
