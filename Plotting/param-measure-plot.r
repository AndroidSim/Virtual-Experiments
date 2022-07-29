#! /usr/bin/Rscript

###
## 	Make a 2D plot: parameter vs measure or measure vs measure
##  Also, make a contour plot: 2 parameters and a measure, 
##  2 measures and a parameter, or 3 measures.
##
## 	Time-stamp: <2021-04-29 aks>
##
###

argv <- commandArgs(T)

PLOT.SVG <- F
#PLOT.AVG <- F 
USE.FRAMES <- T
MA.WINDOW <- 181

require(stats) # for statistics
tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

library(ggplot2)
#library(lattice)
#library(gridExtra)
#library(svglite)
#library(plotly)
#library(plot3D)
#install.packages("rgl", dependencies = TRUE)
#library(rgl)
#library(scatterplot3d)
library(akima)

usage <- function() {
    print("Usage: param-measure-plot.r <table2 of params & measures> <.table2...> ...")
    print("  e.g. param-measure-plot.r pBind_pMetab_AUC.csv")
    quit()
}
if (length(argv) < 1) usage()

###
## test for and create graphics subdirectory
###
if (!file.exists("graphics")) dir.create("graphics")

'%!in%' <- function(x,y)!('%in%'(x,y))

parse.read.file <- function(fname) {
	## get directory name
	dname <- dirname(fname)
	## parse filenames
	fbname <- basename(fname)
	#tmp1 <- substr(files.mnames, 0, regexpr('.csv',files.mnames)-1)
	fbname <- unlist(strsplit(fbname,".csv",fixed=TRUE))
	#coords <- unlist(strsplit(fbname,"-"))
	#pwdf <- paste(getwd(),fbname,sep="/")
	if (file.exists(fname)) {
		dat <- read.csv(file=fname, check.names=FALSE, colClasses="numeric")
	} else {
		cat("\n")
		print(paste("file",fname,"in directory",dname,"doesn't exist. skipping."))
		return(NA)
	}
	#return(list(coords=coords,fbname=fbname,dat=dat))
	return(list(fbname=fbname,dat=dat))
}

D.plot <- function(pdat, use.frames) {
	coords <- colnames(pdat)
	ncoords <- length(coords)
	if (ncoords == 2) {
		## 2D plot
		# x-axis = coordinate 1
		xlabel <- coords[1]
		x <- dat[,xlabel]
		# y-axis = coordinate 2
		ylabel <- coords[2]
		y <- dat[,ylabel]
		xlim <- c(min(x),max(x))
		ylim <- c(min(y),max(y))
		plot(x, y, xlab=xlabel, ylab=ylabel, xlim=xlim, ylim=ylim, 
			#type="p",
			pch=16,
			cex=5)
		#lines(x, y, lwd=5)
	} else if (ncoords == 3) {
		## 3D plot
		#dat <- dat[order(dat[,1]),]
		
		# x-axis = coordinate 1
		xlabel <- coords[1]
		x <- dat[,xlabel]
		# y-axis = coordinate 2
		ylabel <- coords[2]
		y <- dat[,ylabel]
		# z-axis = coordinate 3
		zlabel <- coords[3]
		z <- dat[,zlabel]
		
		xlim <- c(min(x),max(x))
		ylim <- c(min(y),max(y))
		
		temp <- interp(x, y, z)
		x <- temp$x
		y <- temp$y
		z <- temp$z
		
		nlevels <- 10
		maxz <- max(z, na.rm=T)
		minz <- min(z, na.rm=T)
		step <- (maxz-minz)/nlevels
		cols <- rev(rainbow(nlevels))
		cols[1] <- "black"
		
		#persp3d(x, y, z)
		#surface3d(x, y, z)
		#image(x, y, z)
		#plot3d(x, y, z)
		#scatterplot3d(x, y, z)
		#contour(x, y, z, col = cols, levels = seq(minz,maxz,step),
		#		main=fname.base, xlab=xlabel, ylab=ylabel)
		if (use.frames) {
			plot.new()
			plot.window(xlim, ylim)
			.filled.contour(x, y, z, col = cols, levels = seq(minz,maxz,step))
		} else {
			filled.contour(x, y, z, col = cols, levels = seq(minz,maxz,step), 
				main=fname.base, xlab=xlabel, ylab=ylabel)
		}
		
	} else {
		print("number of coordinates not 2 or 3. No plotting done.")
		quit()
	}
}

files <- argv
nfiles <- length(files)
#coords <- list()
fnames.base <- list()
all.dat <- list()
pndx <- 1
for (f in files) {
	out <- parse.read.file(f)
	if (any(is.na(out))) next
	#coords[[pndx]] <- out[[1]]
	fnames.base[[pndx]] <- out[[1]]
	all.dat[[pndx]] <- out[[2]]
	pndx <- pndx + 1
}
nplots <- length(all.dat)
if (USE.FRAMES) {
	## determine # of plot frames
	plot.cols <- round(sqrt(nplots))
	## add a new row if we rounded up
	plot.rows <- ifelse(plot.cols >= sqrt(nplots), plot.cols, plot.cols+1)
	outFile <- paste(unique(unlist(fnames.base)),collapse="-")
	outFile <- paste("graphics/",outFile,"-pmplots",sep="")
	#par(mar=c(5,6,4,2), cex.main=2, cex.axis=2, cex.lab=2)
	if (PLOT.SVG) {
		svg(paste(outFile,".svg",sep=""),10,10)
		par(mfrow=c(plot.rows,plot.cols), cex=2, lwd=3, mar=c(5,6,4,2), cex.main=1, cex.axis=1, cex.lab=1)
	} else {
		png(paste(outFile, ".png", sep=""), 1600, 1600)
		#par(mfrow=c(plot.rows,plot.cols), lwd=3, mar=c(5,6,4,2), cex.main=2, cex.axis=2, cex.lab=2)
		par(mfrow=c(plot.rows,plot.cols), lwd=3, mar=c(1,1,1,1), cex.main=2, cex.axis=2, cex.lab=2)
	}
}
pb <- txtProgressBar(min=0,max=nplots,style=3)
setTxtProgressBar(pb,0)
for (p in c(1:nplots)) {
	fname.base <- fnames.base[[p]]
	dat <- all.dat[[p]]
	
	## set margins and title, axis, and label font sizes
	##   format is c(bottom, left, top, right)
	if (!USE.FRAMES) {
		outFile <- paste("graphics/", fname.base,"-pmplot", sep="")
		if (PLOT.SVG) {
			svg(paste(outFile,".svg",sep=""),10,10)
			par(lwd=3, mar=c(5,6,4,2), las=1, cex.main=2, cex.axis=1, cex.lab=1)
		} else {
			png(paste(outFile, ".png", sep=""), 1600, 1600)
			#par(lwd=3, mar=c(5,10,4,2), las=2, cex.main=4, cex.axis=3, cex.lab=3)
			par(lwd=3, mar=c(5,6,4,2), cex.main=4, cex.axis=2, cex.lab=3)
		}
	}
	
	D.plot(dat, USE.FRAMES)
	
	setTxtProgressBar(pb,getTxtProgressBar(pb)+1)
}
close(pb)
quit()
