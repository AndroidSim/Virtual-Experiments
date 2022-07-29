#! /usr/bin/Rscript

###
## make a 3D plot of space (bands) vs time vs measurements
##
## Time-stamp:
##

PLOT.SVG <- F
PLOT.SURF <- F
PLOT.AVG <- F # either plot moving average or raw data, not both
MA.WINDOW <- 181
tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

argv <- commandArgs(T)

if (length(argv) < 3) {
    print("Usage: plot-xyz.r <xpt dPV> <maxd> <bandsize> <measurement> <exp>")
    print("  e.g. plot-xyz.r 33 76 5 celladj Mouse")
    print("  band files for measurement from get-bands.r must be in -reduced dir")
    quit()
}

library(stringr)
library(ggplot2)
#library(svglite)
#install.packages("plotly")
#library(plotly)
#library(plot3D)
#library(lattice)

'%!in%' <- function(x,y)!('%in%'(x,y))

###
## test for and create graphics subdirectory
###
if (!file.exists("graphics")) dir.create("graphics")

xpt <- as.numeric(argv[1])
maxd <- as.numeric(argv[2])
bandw <- as.numeric(argv[3])
measure <- argv[4]
exps <- argv[-(1:4)] ## all remaining args
directions <- c("dCV", "dPV")
measureOptions <- c("celladj", "mobileObject", "entries", "exits", "rejects", "traps","necrotic","nectrig","entries-avg-pHPC-pMC")
if (!any(str_detect(measure,measureOptions))) {
	print(paste(measure,"is not a valid measurement for plotting"))
	quit()
}	
measurevComs <- c("celladj", "mobileObject", "entries", "exits", "rejects", "traps")
measureEvents <- c("Num_Cells", "Cumu_Necrotic","Cumu_Nectrig")

for (e in exps) {
	print(paste("Working on", e))

    if (!file.exists(paste(e,"-reduced/",sep=""))) {
		print(paste(paste(e,"-reduced/",sep=""),"doesn't exist."))
        next
    }
	
	if (any(grepl(measure,c("necrotic","nectrig")))) {
		dPVtime <- vector("list")
		dCVtime <- vector("list")
	} 
	
	## read data from band files for the dPV section and put into data list
	dPVdata <- vector("list")
	dPVbands <- vector()
	timeisset <- F
	bandstart <- 0
	bandstop <- bandw
	ifile <- 1
	while (bandstop <= xpt) {
		band <- paste("-dPV","∈[",bandstart,",",bandstop,")",sep="")
		if (measure == "exposure") {
			infile <- paste(e,"_","entries-avg-pHPC-pMC",band,"-exposure",sep="")
		} else {
			infile <- paste(e,"_",measure,band,sep="")
		}
		if (!file.exists(paste(e,"-reduced/",infile,".csv",sep=""))) {
			bandstart <- bandstart + bandw
			bandstop <- bandstop + bandw
			next
		}
		dat <- read.csv(paste(e,"-reduced/",infile,".csv",sep=""),check.names=F)
		## get time and the rest of the columns
		if (any(grepl(measure,c("necrotic","nectrig")))) {
			## each data from file has different time column
			dPVtime[[ifile]] <- dat[1]
			if (grepl(measure,"necrotic")) {
				dPVdata[[ifile]] <- dat[,c("Num_Cells","Cumu_Necrotic")]
			}
			if (grepl(measure,"nectrig")) {
				dPVdata[[ifile]] <- dat[,c("Num_Cells","Cumu_Nectrig")]
			}
		} else {
			## time column, same for each data in file
			if (timeisset == F) {
				dPVtime <- dat[1]
				timeisset <- T
			}
			if (PLOT.AVG) {
				## create the moving average data set
				if (nrow(dat) < MA.WINDOW) {
					MA.WINDOW.NEW <- nrow(dat)/4
					cat("WARNING! MA Window of",MA.WINDOW,"is longer than series. Using window of",MA.WINDOW.NEW,"\n")
					MA.WINDOW <- MA.WINDOW.NEW
				}
				dat.tmp <- dat
				dat.tmp[is.na(dat.tmp)] <- 0 # replace NAs with zeros
				dat.ma <- ma.cent(dat.tmp[,2:ncol(dat.tmp)], n=MA.WINDOW)
				dat.ma <- as.data.frame(dat.ma)
				colnames(dat.ma) <- colnames(dat[2:(ncol(dat))])
				dPVdata[[ifile]] <- dat.ma
			} else {
				dPVdata[[ifile]] <- dat[2:(ncol(dat))]
			}
		}
		## the rest of the columns
		if (PLOT.AVG) {
			## create the moving average data set
			if (nrow(dat) < MA.WINDOW) {
				MA.WINDOW.NEW <- nrow(dat)/4
				cat("WARNING! MA Window of",MA.WINDOW,"is longer than series. Using window of",MA.WINDOW.NEW,"\n")
				MA.WINDOW <- MA.WINDOW.NEW
			}
			dat.tmp <- dat
			dat.tmp[is.na(dat.tmp)] <- 0 # replace NAs with zeros
			dat.ma <- ma.cent(dat.tmp[,2:ncol(dat.tmp)], n=MA.WINDOW)
			dat.ma <- as.data.frame(dat.ma)
			colnames(dat.ma) <- colnames(dat[2:(ncol(dat))])
			dPVdata[[ifile]] <- dat.ma
		} else {
			dPVdata[[ifile]] <- dat[2:(ncol(dat))]
		}
		
		dPVbands[ifile] <- band
		ifile <- ifile + 1
		bandstart <- bandstart + bandw
		bandstop <- bandstop + bandw
	}
	
	## read data from band files for the dCV section and put into data list
	dCVdata <- vector("list")
	dCVbands <- vector()
	timeisset <- F
	bandstart <- 0
	bandstop <- bandw
	ifile <- 1
	while (bandstop <= (maxd-xpt)) {
		band <- paste("-dCV","∈[",bandstart,",",bandstop,")",sep="")
		if (measure == "exposure") {
			infile <- paste(e,"_","entries-avg-pHPC-pMC",band,"-exposure",sep="")
		} else {
			infile <- paste(e,"_",measure,band,sep="")
		}
		if (!file.exists(paste(e,"-reduced/",infile,".csv",sep=""))) {
			bandstart <- bandstart + bandw
			bandstop <- bandstop + bandw
			next
		}
		dat <- read.csv(paste(e,"-reduced/",infile,".csv",sep=""),check.names=F)
		## get time and the rest of the columns
		if (any(grepl(measure,c("necrotic","nectrig")))) {
			## each data from file has different time column
			dCVtime[[ifile]] <- dat[1]
			if (grepl(measure,"necrotic")) {
				dCVdata[[ifile]] <- dat[,c("Num_Cells","Cumu_Necrotic")]
			}
			if (grepl(measure,"nectrig")) {
				dCVdata[[ifile]] <- dat[,c("Num_Cells","Cumu_Nectrig")]
			}
		} else {
			## time column, same for each data in file
			if (timeisset == F) {
				dCVtime <- dat[1]
				timeisset <- T
			}
			## the rest of the columns
			if (PLOT.AVG) {
				## create the moving average data set
				if (nrow(dat) < MA.WINDOW) {
					MA.WINDOW.NEW <- nrow(dat)/4
					cat("WARNING! MA Window of",MA.WINDOW,"is longer than series. Using window of",MA.WINDOW.NEW,"\n")
					MA.WINDOW <- MA.WINDOW.NEW
				}
				dat.tmp <- dat
				dat.tmp[is.na(dat.tmp)] <- 0 # replace NAs with zeros
				dat.ma <- ma.cent(dat.tmp[,2:ncol(dat.tmp)], n=MA.WINDOW)
				dat.ma <- as.data.frame(dat.ma)
				colnames(dat.ma) <- colnames(dat[2:(ncol(dat))])
				dCVdata[[ifile]] <- dat.ma
			} else {
				dCVdata[[ifile]] <- dat[2:(ncol(dat))]
			}
		}
		## the rest of the columns
		if (PLOT.AVG) {
			## create the moving average data set
			if (nrow(dat) < MA.WINDOW) {
				MA.WINDOW.NEW <- nrow(dat)/4
				cat("WARNING! MA Window of",MA.WINDOW,"is longer than series. Using window of",MA.WINDOW.NEW,"\n")
				MA.WINDOW <- MA.WINDOW.NEW
			}
			dat.tmp <- dat
			dat.tmp[is.na(dat.tmp)] <- 0 # replace NAs with zeros
			dat.ma <- ma.cent(dat.tmp[,2:ncol(dat.tmp)], n=MA.WINDOW)
			dat.ma <- as.data.frame(dat.ma)
			colnames(dat.ma) <- colnames(dat[2:(ncol(dat))])
			dCVdata[[ifile]] <- dat.ma
		} else {
			dCVdata[[ifile]] <- dat[2:(ncol(dat))]
		}
		
		dCVbands[ifile] <- band
		ifile <- ifile + 1
		bandstart <- bandstart + bandw
		bandstop <- bandstop + bandw
	}
	
	if (any(grepl(measure,c("necrotic","nectrig")))) {
		dPVVCd <- vector("list")
		if (length(dPVdata) > 0 && length(dCVdata) == 0) {
			minDeathTime = min(unlist(lapply(dPVtime, function(x) {min(x)})))
			maxDeathTime = max(unlist(lapply(dPVtime, function(x) {max(x)})))
			DeathTime <- minDeathTime:maxDeathTime
			columns <- colnames(dPVdata[[1]])
			dPVdat <- vector("list")
			DTdPVdata <- vector("list")
			for (band in 1:length(dPVdata)) {		
				DTdPVdata[[band]] <- data.frame(matrix(0,nrow=length(DeathTime),ncol=ncol(dPVdata[[band]])))
				DTdPVdata[[band]][match(as.matrix(dPVtime[[band]]),DeathTime),] <- dPVdata[[band]]
				colnames(DTdPVdata[[band]]) <- colnames(dPVdata[[band]])
			}
			for (i in 1:ncol(DTdPVdata[[1]])) {
				dPVdat[[i]] <- sapply(DTdPVdata, function(x) { x[[i]] })
				colnames(dPVdat[[i]]) <- paste(columns[i],dPVbands,sep="")
			}
			dtype <- "dPV"
			disttime <- DeathTime
			dPVVCd <- dPVdat
		}
		if (length(dPVdata) == 0 && length(dCVdata) > 0) {
			minDeathTime = min(unlist(lapply(dCVtime, function(x) {min(x)})))
			maxDeathTime = max(unlist(lapply(dCVtime, function(x) {max(x)})))
			DeathTime <- minDeathTime:maxDeathTime
			columns <- colnames(dCVdata[[1]])
			dCVdat <- vector("list")
			DTdCVdata <- vector("list")
			for (band in 1:length(dCVdata)) {		
				DTdCVdata[[band]] <- data.frame(matrix(0,nrow=length(DeathTime),ncol=ncol(dCVdata[[band]])))
				DTdCVdata[[band]][match(as.matrix(dCVtime[[band]]),DeathTime),] <- dCVdata[[band]]
				colnames(DTdCVdata[[band]]) <- colnames(dCVdata[[band]])
			}
			for (i in 1:ncol(DTdCVdata[[1]])) {
				dCVdat[[i]] <- sapply(DTdCVdata, function(x) { x[[i]] })
				colnames(dCVdat[[i]]) <- paste(columns[i],dCVbands,sep="")
			}
			dtype <- "dCV"
			disttime <- DeathTime
			dPVVCd <- dCVdat
		}
		if (length(dPVdata) > 0 && length(dCVdata) > 0) {
			mindPVDeathTime = min(unlist(lapply(dPVtime, function(x) {min(x)})))
			maxdPVDeathTime = max(unlist(lapply(dPVtime, function(x) {max(x)})))
			mindCVDeathTime = min(unlist(lapply(dCVtime, function(x) {min(x)})))
			maxdCVDeathTime = max(unlist(lapply(dCVtime, function(x) {max(x)})))
			if (mindPVDeathTime < mindCVDeathTime) {
				minDeathTime <- mindPVDeathTime
			} else{
				minDeathTime <- mindCVDeathTime
			}
			if (maxdPVDeathTime > maxdCVDeathTime) {
				maxDeathTime <- maxdPVDeathTime
			} else{
				maxDeathTime <- maxdPVDeathTime
			}
			DeathTime <- minDeathTime:maxDeathTime
			columns <- colnames(dPVdata[[1]])
			dPVdat <- vector("list")
			DTdPVdata <- vector("list")
			for (band in 1:length(dPVdata)) {		
				DTdPVdata[[band]] <- data.frame(matrix(0,nrow=length(DeathTime),ncol=ncol(dPVdata[[band]])))
				DTdPVdata[[band]][match(as.matrix(dPVtime[[band]]),DeathTime),] <- dPVdata[[band]]
				colnames(DTdPVdata[[band]]) <- colnames(dPVdata[[band]])
			}
			for (i in 1:ncol(DTdPVdata[[1]])) {
				dPVdat[[i]] <- sapply(DTdPVdata, function(x) { x[[i]] })
				colnames(dPVdat[[i]]) <- paste(columns[i],dPVbands,sep="")
			}
			columns <- colnames(dCVdata[[1]])
			dCVdat <- vector("list")
			DTdCVdata <- vector("list")
			for (band in 1:length(dCVdata)) {		
				DTdCVdata[[band]] <- data.frame(matrix(0,nrow=length(DeathTime),ncol=ncol(dCVdata[[band]])))
				DTdCVdata[[band]][match(as.matrix(dCVtime[[band]]),DeathTime),] <- dCVdata[[band]]
				colnames(DTdCVdata[[band]]) <- colnames(dCVdata[[band]])
			}
			for (i in 1:ncol(DTdCVdata[[1]])) {
				dCVdat[[i]] <- sapply(DTdCVdata, function(x) { x[[i]] })
				colnames(dCVdat[[i]]) <- paste(columns[i],dCVbands,sep="")
			}
			dtype <- "dPVVCd"
			disttime <- DeathTime
			if (length(dPVdat) != length(dCVdat)) {
				print("length(dPVdat) != length(dCVdat)")
				quit()
			}
			for (i in 1:length(dPVdat)) {
				dPVdf <- dPVdat[[i]]
				dCVdf <- dCVdat[[i]]
				dPVVCd[[i]] <- cbind(dPVdf,dCVdf[,rev(colnames(dCVdf))])
			}
		}
	} else {
		# for each vCompound/MITs column condense the list of data frames into one
		# data frame for each dCV and dPV section.
		if (length(dPVdata) > 0) { 
			columns <- colnames(dPVdata[[1]])
			dPVdat <- vector("list")
			for (i in 1:ncol(dPVdata[[1]])) {
				dPVdat[[i]] <- sapply(dPVdata, function(x) { x[[i]] })
				colnames(dPVdat[[i]]) <- paste(columns[i],dPVbands,sep="")
			}
		}
		if (length(dCVdata) > 0) {
			columns <- colnames(dCVdata[[1]])
			dCVdat <- vector("list")
			for (i in 1:ncol(dCVdata[[1]])) {
				dCVdat[[i]] <- sapply(dCVdata, function(x) { x[[i]] })
				colnames(dCVdat[[i]]) <- paste(columns[i],dCVbands,sep="")
			}
		}
		## reverse or flip dCV section and attach to dPV section
		## so dPV section : reverse (dCV section)
		## or dPV from 0 -> xpt : dCV from xpt -> 0
		# by now length(dPVdat) should = length(dCVdat)
		dPVVCd <- vector("list")
		if (length(dPVdata) > 0 && length(dCVdata) == 0) {
			dtype <- "dPV"
			disttime <- dPVtime
			dPVVCd <- dPVdat
		}
		if (length(dPVdata) == 0 && length(dCVdata) > 0) {
			dtype <- "dCV"
			disttime <- dCVtime
			dPVVCd <- dCVdat
		}
		if (length(dPVdata) > 0 && length(dCVdata) > 0) {
			dtype <- "dPVVCd"
			disttime <- dPVtime
			if (length(dPVdat) != length(dCVdat)) {
				print("length(dPVdat) != length(dCVdat)")
				quit()
			}
			for (i in 1:length(dPVdat)) {
				dPVdf <- dPVdat[[i]]
				dCVdf <- dCVdat[[i]]
				dPVVCd[[i]] <- cbind(dPVdf,dCVdf[,rev(colnames(dCVdf))])
			}
		}
	}
		
	## plot each data frame with the rows for time and total number of bands for space
	pb <- txtProgressBar(min=0,max=length(dPVVCd),style=3)
	setTxtProgressBar(pb,1)

	for (i in 1:length(dPVVCd)) {
		## write the data to a file
		tmp <- cbind(disttime, dPVVCd[[i]])
		colnames(tmp)[1] <- "Time"
		write.csv(x=tmp, file=paste(e, "_",measure,"-",columns[i],"-",dtype,"-xyz.csv", sep=""), row.names=FALSE)
		
		# x-axis = time, dPVtime = dCVtime
		x <- t(disttime)
		# y-axis = space, as band index, dPV section: 1:index(xpt), 
		# dCV section: index(xpt):total number of bands
		y <- 1:(length(dPVbands)+length(dCVbands))
		# z-axis = data
		z <- dPVVCd[[i]]
		
		if (PLOT.SVG) {
			if (PLOT.SURF) {
				outFile <- paste("graphics/", e, "_",measure,"-",columns[i],"-surface", sep="")
				svg(paste(outFile,".svg",sep=""),10,10)
				par(cex=2, lwd=3, mar=c(5,6,4,2), cex.main=1, cex.axis=1, cex.lab=1)
				persp(x,y,z)
				title(paste(e, "_",measure,"-",columns[i],"-surface", sep=""))
			} else {
				outFile <- paste("graphics/", e, "_",measure,"-",columns[i],"-contour", sep="")
				svg(paste(outFile,".svg",sep=""),10,10)
				par(cex=2, lwd=3, mar=c(5,6,4,2), cex.main=1, cex.axis=1, cex.lab=1)
				contour(x,y,z)
				#image(x,y,z)
				title(paste(e, "_",measure,"-",columns[i],"-contour", sep=""))
			}
			
		} else {
			if (PLOT.SURF) {
				outFile <- paste("graphics/", e, "_",measure,"-",columns[i],"-surface", sep="")
				png(paste(outFile, ".png", sep=""), 1600, 1600)
				par(cex=2, lwd=3, mar=c(5,6,4,2), cex.main=1, cex.axis=1, cex.lab=1)
				persp(x,y,z)
				title(paste(e, "_",measure,"-",columns[i],"-surface", sep=""))
			} else {
				outFile <- paste("graphics/", e, "_",measure,"-",columns[i],"-",dtype,"-contour", sep="")
				png(paste(outFile, ".png", sep=""), 1600, 1600)
				nlevels <- 10
				maxz <- max(z, na.rm=T)
				minz <- min(z, na.rm=T)
				step <- (maxz-minz)/nlevels
				ylabels <- colnames(dPVVCd[[i]])
				par(lwd=3, mar=c(5,17,4,2), las=2, cex.main=4, cex.axis=2, cex.lab=3)
				cols <- rev(rainbow(nlevels))
				cols[1] <- "black"
				if (maxz == 0 && minz == 0) {
					filled.contour(x,y,z, col = cols)
				} else {
					filled.contour(x,y,z, col = cols, levels = seq(minz,maxz,step))
				}
				axis(2, at=y, labels=ylabels)
				title(paste(e, "_",measure,"-",columns[i],"-",dtype,"-contour", sep=""))
			}	
		}
		setTxtProgressBar(pb,getTxtProgressBar(pb)+1)
	}
	close(pb)
} # loop over experiments

quit()
