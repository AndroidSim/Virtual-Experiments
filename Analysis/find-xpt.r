#! /usr/bin/Rscript

###
## finds the intersection point of the two curves, avg(vHPC) vs dCV
## and avg(HPC) vs dPV
##
## Time-stamp: <2019-07-11 09:46:37 gepr>
###

argv <- commandArgs(T)
AVG <- T # if AVG or SUM
PLOT.SVG <- F

require(stats) # for statistics

usage <- function() {
	print("Usage: find-midpt.r <exp directories>")
	print("  directories should contain a files experiment_Hcounts-all-dCV<dPV>.csv, ")
	quit()
}
if (length(argv) < 1) usage()

tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

NApad1stColumns <- function(first, second) {
	## insert columns of NAs that don't exist in first
	matched.columns <- match(colnames(second), colnames(first))
	for (col.index in which(is.na(matched.columns))) {
		newcol <- vector(mode="numeric",length=length(first[,1]))
		newcol[] <- NA
		total.names <- colnames(first)
		first <- cbind(first, newcol)
		colnames(first) <- c(total.names, colnames(second)[col.index])
	}
	return(first)
}

##
## test for and create graphics subdirectory
##
#if (!file.exists("graphics")) dir.create("graphics")

## the following R code for finding the intersection point of two line
## segments was obtained from the following web address:
## https://stackoverflow.com/questions/20519431/finding-point-of-intersection-in-r

## segment-segment intersection code
## http://paulbourke.net/geometry/pointlineplane/
ssi <- function(x1, x2, x3, x4, y1, y2, y3, y4){
	denom <- ((y4 - y3)*(x2 - x1) - (x4 - x3)*(y2 - y1))
	denom[abs(denom) < 1e-10] <- NA # parallel lines

	ua <- ((x4 - x3)*(y1 - y3) - (y4 - y3)*(x1 - x3)) / denom
	ub <- ((x2 - x1)*(y1 - y3) - (y2 - y1)*(x1 - x3)) / denom

	x <- x1 + ua * (x2 - x1)
	y <- y1 + ua * (y2 - y1)

	inside <- (ua >= 0) & (ua <= 1) & (ub >= 0) & (ub <= 1)
	data.frame(x = ifelse(inside, x, NA),y = ifelse(inside, y, NA))
}

exps <- argv
for (e in exps) {
	e <- paste(dirname(e), basename(e), sep="/") # removes trailing slash from the argument
	print(paste("Working on", e))

	if (!file.exists(paste(e,"-reduced/",sep=""))) {
		print(paste(e,"doesn't exist."))
		next
	}
	
	## get data from files
	vHPC.dCV <- read.csv(paste(e, "-reduced/", e, "_Hcounts-all-dCV.csv", sep=""), row.names=1, check.names=F)
	vHPC.dPV <- read.csv(paste(e, "-reduced/", e, "_Hcounts-all-dPV.csv", sep=""), row.names=1, check.names=F)
	
	rnames.dCV <- rownames(vHPC.dCV)
	rnames.dPV <- rownames(vHPC.dPV)
	cnames.dCV <- colnames(vHPC.dCV)
	cnames.dPV <- colnames(vHPC.dPV)
	
	## resolve column (distance value) differences
	#dCV <- colnames(vHPC.dCV[,1:(ncol(vHPC.dCV)-3)])
	dCV <- head(cnames.dCV,-3)
	#dPV <- colnames(vHPC.dPV[,1:(ncol(vHPC.dPV)-3)])
	dPV <- head(cnames.dPV,-3)
	d.maxl <- if (length(dCV) > length(dPV)) dCV else dPV
	
	if (AVG) {
		#vHPC.dCV <- vHPC.dCV[,1:(ncol(vHPC.dCV)-3)]
		#vHPC.dPV <- vHPC.dPV[,1:(ncol(vHPC.dPV)-3)]

		vHPC.dCV <- as.data.frame(vHPC.dCV[,1:(ncol(vHPC.dCV)-3)], row.names=rnames.dCV)
		colnames(vHPC.dCV) <- dCV
		vHPC.dPV <- as.data.frame(vHPC.dPV[,1:(ncol(vHPC.dPV)-3)], row.names=rnames.dPV)
		colnames(vHPC.dPV) <- dPV
		
		dPVdCV <- as.data.frame(NApad1stColumns(vHPC.dPV["avg",], vHPC.dCV["avg",]))
		colnames(dPVdCV) <- d.maxl
		dCVdPV <- as.data.frame(NApad1stColumns(vHPC.dCV["avg",], vHPC.dPV["avg",]))
		colnames(dCVdPV) <- d.maxl
	} else {
		vHPCtotal.dPV <- vHPC.dPV["total","total"]
		vHPCtotal.dCV <- vHPC.dCV["total","total"]
		#vHPC.dCV <- vHPC.dCV[,1:(ncol(vHPC.dCV)-3)]
		#vHPC.dPV <- vHPC.dPV[,1:(ncol(vHPC.dPV)-3)]
		
		vHPC.dCV <- as.data.frame(vHPC.dCV[,1:(ncol(vHPC.dCV)-3)], row.names=rnames.dCV)
		colnames(vHPC.dCV) <- dCV
		vHPC.dPV <- as.data.frame(vHPC.dPV[,1:(ncol(vHPC.dPV)-3)], row.names=rnames.dPV)
		colnames(vHPC.dPV) <- dPV
		
		dPVdCV <- as.data.frame(NApad1stColumns(vHPC.dPV["total",], vHPC.dCV["total",]))
		dPVdCV <- dPVdCV/vHPCtotal.dPV
		colnames(dPVdCV) <- d.maxl
		dCVdPV <- as.data.frame(NApad1stColumns(vHPC.dCV["total",], vHPC.dPV["total",]))
		dCVdPV <- dCVdPV/vHPCtotal.dCV
		colnames(dCVdPV) <- d.maxl
	}
	dat <- rbind(dPVdCV, dCVdPV)

	## define data as two curves each composed of a set of points
	## the curves most likely will not have the same number of points
	distances <- as.numeric(colnames(dat))
	## first row is "avgHPC-dPV", second row is "avgHPC-dCV"
	##matrix(avgHcts_pD,nrow=1,ncol=length(as.vector(avgHcts_pD)))
	avgdPV <- matrix(dat[1, ],nrow=1,ncol=length(dat[1, ]))
	avgdCV <- matrix(dat[2, ],nrow=1,ncol=length(dat[2, ]))
	ndpts <- length(distances)
	## change any distance with a NA to the number of the previous distance
	if (any(is.na(avgdPV))) {
		for (i in which(is.na(avgdPV))) {
			avgdPV[i] <- avgdPV[i-1]
		}
	}
	if (any(is.na(avgdCV))) {
		for (i in which(is.na(avgdCV))) {
			avgdCV[i] <- avgdCV[i-1]
		}
	}

	## start with the first line segment of one of the curves, then loop
	## over the line segments of the other curve and test for intersection
	## if no intersection, go to the next line segment of the first curve
	## continue moving along first curve until intersection point is reached
	y1 <- avgdPV
	y2 <- avgdCV
	y <- t(rbind(y1, y2))
	x1 <- distances
	x2 <- rev(distances)
	x <- t(rbind(x1,x2))
	dPVline <- cbind(x[,1],y[,1])
	dCVline <- cbind(x[,2],y[,2])

	xpt <- data.frame(x = NA,y = NA)
	notxpt <- T
	idCV <- 1:2
	while (notxpt && all(idCV <= length(distances))) {
		x3 <- dCVline[[idCV[1],1]]
		y3 <- dCVline[[idCV[1],2]]
		x4 <- dCVline[[idCV[2],1]]
		y4 <- dCVline[[idCV[2],2]]
		for (idPV in 1:(ndpts-1)) {
			## point1 = [x1,y1], point2 = [x2,y2], etc
			## point1 & 2 on the dPV line, and point3 & 4 on the dCV line
			x1 <- dPVline[[idPV,1]]
			y1 <- dPVline[[idPV,2]]
			x2 <- dPVline[[idPV+1,1]]
			y2 <- dPVline[[idPV+1,2]]
			xpt <- ssi(x1, x2, x3, x4, y1, y2, y3, y4)
			if (all(!is.na(xpt))) {
				notxpt <- F
				break
			}
		}
		idCV <- idCV+1
	}

	## output the intersection point and plot the two curves with xpt
	## The intersection point is data frame with x,y = columns
	## so dPV = x and dCV = max distance - x , and average(vHPCs) = y
	print("intersection point: dPV = x, dCV = maxd-x, avg(vHPCs) = y")
	maxd <- max(distances)
	xptreport <- cbind(xpt, maxd, maxd-xpt[1])
	colnames(xptreport) <- c("dPV", "avg(vHPCs)", "max", "dCV")
	print(xptreport)
	write.csv(xptreport,
            paste(e, "-reduced/", e, "_xpt.csv", sep=""),
            row.names=F)

	if (AVG) {
		contents <- "AVG"
	} else {
		contents <- "FRACTION"
	}
	#outfile <- paste(e,"_","xpt∈[dPV&dCV)-", contents, sep="")
	#outFile <- paste("graphics/",outfile, sep="")
	#if (PLOT.SVG) {
	#	svg(paste(outFile,".svg",sep=""),10,10)
	#	matplot(x, y, t="l", lty=1, , lwd=2, col=2:3)
	#	points(xpt, pch=21, col=1, bg=1)
	#} else {
	#	png(paste(outFile, ".png", sep=""), 1600, 1600)
	#	matplot(x, y, t="l", lty=1, , lwd=4, col=2:3)
	#	points(xpt, pch=21, col=1, bg=1, cex=3)
	#}
} # loop over experiments

quit()
