#! /usr/bin/Rscript

###
## 	Make a 2D plot of 3D data structured as a matrix with space (i.e. bands) vs time 
## 	vs measurements (or events). The matrix is inputted from a file as an argument.
##
## 	Time-stamp: <2020-09-11 aks>
##
###

argv <- commandArgs(T)

BARPLOT <- T
PLOT.SVG <- F
#PLOT.AVG <- F 
USE.FRAMES <- T
MA.WINDOW <- 181

require(stats) # for statistics
tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

library(ggplot2)
#library(svglite)
#library(plotly)
#library(plot3D)
#library(lattice)

usage <- function() {
    print("Usage: plot-slice.r time <start time> <end time> <*-xyz.csv file of experiment1> <...experiment2> ...")
    print("	   or plot-slice.r distance/dCV/dPV <start distance> <end distance> <*-xyz.csv file of experiment1> <...experiment2> ...")
    print("  e.g. plot-slice.r time 5000 6000 Mouse_entries-Compound-dCV-xyz.csv")
    quit()
}
if (length(argv) < 4) usage()

###
## test for and create graphics subdirectory
###
if (!file.exists("graphics")) dir.create("graphics")

'%!in%' <- function(x,y)!('%in%'(x,y))

parse.read.file <- function(fname) {
	fullfn <- paste(dirname(fname), basename(fname), sep="/") 
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

if (argv[1] %in% c("distance","dCV","dPV")) {
	dist1 <- as.numeric(argv[2])
	dist2 <- as.numeric(argv[3])
	if (dist1 < dist2) {
		start.dist <- dist1
		end.dist <- dist2
	} else {
		start.dist <- dist2
		end.dist <- dist1
	}
	time.slice <- F
} else if (argv[1] %in% c("time")) {
	time1 <- as.numeric(argv[2])
	time2 <- as.numeric(argv[3])
	if (time1 < time2) {
		start.time <- time1
		end.time <- time2
	} else {
		start.time <- time2
		end.time <- time1
	}
	time.slice <- T
} else {
	print(paste(argv[1],"is not distance, dCV, dPV, or time"))
	quit()
}
Files <- argv[-(1:3)]
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
if (nplots == 0) {
	print("no data to plot")
	quit()
}
if (nplots == 1) USE.FRAMES <- F
if (USE.FRAMES) {
	## determine # of plot frames
	plot.cols <- round(sqrt(nplots))
	## add a new row if we rounded up
	plot.rows <- ifelse(plot.cols >= sqrt(nplots), plot.cols, plot.cols+1)
	fname <- c(unique(unlist(xnames)),unique(unlist(mnames)))
	outFile <- paste(fname,collapse="-")
	outFile <- paste("graphics/",outFile,"-slices",sep="")
	#par(mar=c(5,6,4,2), cex.main=2, cex.axis=2, cex.lab=2)
	if (PLOT.SVG) {
		svg(paste(outFile,".svg",sep=""),10,10)
		par(mfrow=c(plot.rows,plot.cols), cex=2, lwd=3, mar=c(5,6,4,2), cex.main=1, cex.axis=1, cex.lab=1)
	} else {
		png(paste(outFile, ".png", sep=""), 1600, 1600)
		par(mfrow=c(plot.rows,plot.cols), lwd=3, mar=c(5,6,4,2), cex.main=2, cex.axis=2, cex.lab=2)
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
	# y-axis = space or distance, each column is a band
	#y <- 1:(ncol(dat)-1)
	bandnames <- colnames(dat[,2:ncol(dat)])
	temp <- unlist(strsplit(bandnames,"∈"))
	direction <- unique(temp[c(T,F)])
	bands <- unique(temp[c(F,T)])
	matches <- regmatches(bands, gregexpr("[[:digit:]]+", bands))
	vbands <- as.numeric(unique(unlist(matches)))
	y <- vbands
	# z-axis = data or measurements
	z <- as.matrix(dat[,2:ncol(dat)])
	#pdat <- as.vector(z[which.max(x),])
	if (time.slice) {
		min.x <- min(x)
		max.x <- max(x)
		if (start.time < min.x) {
			cat("\n")
			print(paste(start.time,"is less than the min time",min.x,"for file",f))
			next
		}
		if (end.time > max.x) {
			cat("\n")
			print(paste(end.time,"is greater than the max time",max.x,"for file",f))
			next
		}
		if (start.time == end.time) {
			pdat <- as.vector(z[start.time == x,])
		} else {
			tndx <- which(start.time <= x & x <= end.time)
			pdat <- as.matrix(z[tndx,])
			pdat <- colMeans(pdat)
		}
		#temp <- as.matrix(z[(start.time <= x & x <= end.time),])
		#temp <- colMeans(temp)
	} else {
		min.y <- min(y)
		max.y <- max(y)
		if (start.dist < min.y) {
			cat("\n")
			print(paste(start.dist,"is less than the min distance",min.y,"for file",f))
			next
		}
		if (end.dist > max.y) {
			cat("\n")
			print(paste(end.dist,"is greater than the max distance",max.y,"for file",f))
			next
		}
		if (start.dist == end.dist) {
			pdat <- as.vector(z[,start.dist])
		} else {
			dndx <- which(start.dist <= y & y <= end.dist)
			#pdat <- as.matrix(z[,dndx])
			pdat <- as.matrix(z[,y[dndx]])
			pdat <- rowMeans(pdat)
		}
	}
	
	## set margins and title, axis, and label font sizes
	##   format is c(bottom, left, top, right)
	if (!USE.FRAMES) {
		## place the right margin 1 unit for each of the characters +  units for the line
		#right.margin <- max(nchar(titles)) + max(nchar(maws)) + ifelse(plot.svg, 0, 12)
		#par(mar=c(5,6,4,right.margin), cex.main=2, cex.axis=2, cex.lab=2)
		outFile <- paste("graphics/",fname.base,"-slice",sep="")
		if (PLOT.SVG) {
			svg(paste(outFile,".svg",sep=""),10,10)
			par(cex=2, lwd=3, mar=c(5,6,4,2), cex.main=1, cex.axis=1, cex.lab=1)
		} else {
			png(paste(outFile, ".png", sep=""), 1600, 1600)
			par(lwd=3, mar=c(5,6,4,2), cex.main=2, cex.axis=2, cex.lab=2)
		}
	}
	
	if (time.slice) {
		xlabel <- direction
		ylabel <- mname
		if (BARPLOT) {
			barplot(pdat, names.arg=vbands[-length(vbands)], main=xname, xlab=xlabel, ylab=ylabel)
			box()
		} else {
			plot(vbands[-length(vbands)], pdat, main=xname, xlab=xlabel, ylab=ylabel)
		}
	} else {
		xlabel <- "Time"
		ylabel <- mname
		if (BARPLOT) {
			barplot(pdat, names.arg=x, main=xname, xlab=xlabel, ylab=ylabel)
			box()
		} else {
			plot(vbands[-length(vbands)], pdat, main=xname, xlab=xlabel, ylab=ylabel)
		}
	}
	#grid()
	setTxtProgressBar(pb,getTxtProgressBar(pb)+1)
}
close(pb)
quit()
