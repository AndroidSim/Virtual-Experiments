#! /usr/bin/Rscript

###
##
##	Creates the data file containing a large matrix:
##		rows=time, columns=space as bands of a certain width
##		and each element some measure (e.g. amount of mobile object)
##	This data file can then be used by plot-xyz.r to make contour plots
##		with x=time, y=space, z=data.
##
##	Time-stamp: <2020-08-24 aks>
##
###

argv <- commandArgs(T)

require(stats) # for statistics
tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

usage <- function() {
  print("Usage: create-xyz.r <dPV|dCV> dMin dMax <band width> <exp directories> ...")
  print("       <dPV|dCV> = either dCV or dPV distance")
  print("       dMin = min dCV or dPV distance")
  print("       dMax = max dCV or dPV distance")
  print("       <band width> = width of distance band")
  print("e.g. create-xyz.r dCV 0 35 1 data-time-format")
  print("  directories should contain files like necrotic-[0-9]+.csv,")
  print("  entries-d[CP]V-[0-9]+.csv, and celladj-d[CP]V-[0-9]+.csv")
  print("  Note: experiment directories are the raw data-time format")
  print("  not the -reduced directory")
  quit()
}
if (length(argv) < 5) usage()

## function to create measure matrix
create.M.mat <- function(mfiles) {
	#totals <- avgByColumn(allfiles)
	#avgByColumn <- function(allfiles)
	#options(warn=1)
	tNdx <- 1
	for (t in mfiles) {
		tdat <- read.csv(t, check.names=F)
		if (tNdx == 1) {
			Time <- tdat[,1]
			tdat <- tdat[,2:ncol(tdat)]
			totals <- tdat
		} else {
			tdat <- tdat[,2:ncol(tdat)]
			totals <- pad1stColumns(totals, tdat)
			tpad <- pad1stColumns(tdat, totals)
			totals <- totals[, order(names(totals))]
			tpad <- tpad[, order(names(tpad))]
			totals <- totals + tpad
		}
		tNdx <- tNdx + 1
	} ## end of loop over trials
	if (any(grepl("rxnprod", mfiles, fixed=TRUE))) {
		distances <- sort(as.numeric(colnames(totals)))
		compounds <- "rxnprod"
	} else {
		dcnames <- unlist(strsplit(colnames(totals),":"))
		distances <- unique(dcnames[c(T,F)])
		distances <- sort(as.numeric(distances))
		compounds <- unique(dcnames[c(F,T)])
	}
	maxdist <- max(distances)
	if (dMin > maxdist) {
		print(paste("dMin",dMin,"is > max distance",maxdist,"for direction",direction))
		return(NA)
	}
	#inband <- sumBandByLastTag(totals,c(dMin,dMax))
	#sumBandByLastTag <- function(dat, band)
	ctotals <- list()
	## loop over compounds
	for (p in compounds) {
		## loop over space or bands
		bandnames <- vector()
		bandstart <- dMin
		bandstop <- bandstart + bandw
		bandx <- 1
		while (bandstop <= dMax) {
			newnames <- vector()
			## find distances within band & form column name to search
			#c(outer(distances,groups, FUN=paste))
			for (d in distances) {
				if (bandstart <= d && d < bandstop) {
					if (grepl("rxnprod", p, fixed=TRUE)) {
						temp <- d
					} else {
						temp <- paste(d,p,sep=":")
					}
					newnames <- c(newnames,temp)
				}
			}
			bandnames <- c(bandnames,paste(direction,"∈[",bandstart,",",bandstop,")",sep=""))
			cbndx <- match(newnames,colnames(totals))
			bandat <- totals[,cbndx]
			if (bandx == 1) {
				if (length(dim(bandat)) > 0) {
					ctotals[[p]] <- rowSums(bandat)
				} else {
					ctotals[[p]] <- bandat
				}
			} else {
				if (length(dim(bandat)) > 0) {
					ctotals[[p]] <- cbind(ctotals[[p]],rowSums(bandat))
				} else {
					ctotals[[p]] <- cbind(ctotals[[p]],bandat)
				}
			}
			bandstart <- bandstart + bandw
			bandstop <- bandstop + bandw
			bandx <- bandx + 1
		}
		## take average over MC trials & add back the Time column
		ctotals[[p]] <- ctotals[[p]]/length(mfiles)
		ctotals[[p]] <- cbind(Time,ctotals[[p]])
		colnames(ctotals[[p]]) <- c("Time",bandnames)		
	} ## end loop over compounds
	return(ctotals)
}

## function to create Hepatocyte (vHPC) matrix
create.H.mat <- function(hfiles) {
	## modification of hcounts-band.r
	#dMin_ord <- dMin+1
	#dMax_ord <- dMax+1
	#total <- 0
	Hband.aT <- vector()
	Trials <- vector()
	trial <- 1
	for (f in hfiles) {
		dat <- read.csv(f, colClasses="numeric", check.names=F)
		#if (dMin_ord > ncol(dat)) next
		#thismax <- ifelse(dMax_ord > ncol(dat), ncol(dat), dMax_ord-1)
		distances <- as.numeric(colnames(dat))
		maxdist <- max(distances)
		if (dMin > maxdist) {
			cat("\n")
			print(paste("dMin",dMin,"is > max distance",maxdist,"for direction",direction,"in file",f))
			#next
		}
		Trials <- c(Trials,trial)
		bandnames <- vector()
		Hband.pT <- vector()
		bandstart <- dMin
		bandstop <- bandstart + bandw
		bandx <- 1
		while (bandstop <= dMax) {
			bandnames <- c(bandnames,paste(direction,"∈[",bandstart,",",bandstop,")",sep=""))
			bandat <- dat[,bandstart <= distances & distances < bandstop]
			if (length(dim(bandat)) > 0) {
				Hband.pT <- c(Hband.pT,rowSums(bandat))
			} else {
				Hband.pT <- c(Hband.pT,bandat)
			}
			#row <- dat[1,dMin_ord:thismax]
			#rowtot <- sum(row)
			#total <- total+rowtot
			bandstart <- bandstart + bandw
			bandstop <- bandstop + bandw
			bandx <- bandx + 1
		}
		Hband.aT <- rbind(Hband.aT,Hband.pT)
		trial <- trial + 1
	} ## end of loop over files
	#avg <- total/length(files)
	Hband.aT <- matrix(Hband.aT, nrow=length(Trials), ncol=ncol(Hband.aT))
	colnames(Hband.aT) <- bandnames
	#Hband.aT <- cbind(Trials,Hband.aT)
    #colnames(Hband.aT) <- c("Trial",bandnames)
	HbaT.totals <- colSums(Hband.aT)
    HbaT.totals <- matrix(HbaT.totals,nrow=1, ncol=length(bandnames))
    colnames(HbaT.totals) <- bandnames
	if (length(Hband.aT) == 0) {
		return(NA)
	} else {
		return(list(Hband.aT=Hband.aT,HbaT.totals=HbaT.totals))
	}
}

## function to create measure per vHPC matrix
create.MpH.mat <- function(avgM.aT,mfiles,Hdat) {
	Hband.aT <- Hdat[[1]]
	HbaT.totals <- Hdat[[2]]
	MpH.pT <- list()
	Trials <- vector()
	trial <- 1
	for (f in mfiles) {
		tdat <- read.csv(f, check.names=F)
		if (trial == 1) {
			Time <- tdat[,1]
			MpH.pT[[trial]] <- tdat[,2:ncol(tdat)]
			#M.allT <- MpH.pT[[trial]]
		} else {
			MpH.pT[[trial]] <- tdat[,2:ncol(tdat)]
			#M.allT <- cbind(M.allT,MpH.pT[[trial]])
		}
		if (grepl("rxnprod", f, fixed=TRUE)) {
			distances <- sort(as.numeric(colnames(MpH.pT[[trial]])))
			compounds <- "rxnprod"
		} else {
			dcnames <- unlist(strsplit(colnames(MpH.pT[[trial]]),":"))
			distances <- unique(dcnames[c(T,F)])
			distances <- sort(as.numeric(distances))
			compounds <- unique(dcnames[c(F,T)])
			compounds <- compounds[order(compounds)]
		}
		maxdist <- max(distances)
		if (dMin > maxdist) {
			cat("\n")
			print(paste("dMin",dMin,"is > max distance",maxdist,"for direction",direction,"in file",f))
			#return(NA)
			#next
		}
		cMband.pT <- list()
		## loop over compounds
		for (p in compounds) {
			## loop over space or bands
			bandnames <- vector()
			bandstart <- dMin
			bandstop <- bandstart + bandw
			bandx <- 1
			while (bandstop <= dMax) {
				newnames <- vector()
				for (d in distances) {
					if (bandstart <= d && d < bandstop) {
						if (grepl("rxnprod", p, fixed=TRUE)) {
							temp <- d
						} else {
							temp <- paste(d,p,sep=":")
						}
						newnames <- c(newnames,temp)
					}
				}
				bandnames <- c(bandnames,paste(direction,"∈[",bandstart,",",bandstop,")",sep=""))
				cbndx <- match(newnames,colnames(MpH.pT[[trial]]))
				bandat <- MpH.pT[[trial]][,cbndx]
				if (bandx == 1) {
					if (length(dim(bandat)) > 0) {
						cMband.pT[[p]] <- rowSums(bandat)
					} else {
						cMband.pT[[p]] <- bandat
					}
				} else {
					if (length(dim(bandat)) > 0) {
						cMband.pT[[p]] <- cbind(cMband.pT[[p]],rowSums(bandat))
					} else {
						cMband.pT[[p]] <- cbind(cMband.pT[[p]],bandat)
					}
				}
				bandstart <- bandstart + bandw
				bandstop <- bandstop + bandw
				bandx <- bandx + 1
			}
			## divided each compound band vs time matrix by number of Hs
			## per each band, Hband.aT
			mat <- as.matrix(cMband.pT[[p]])
			vec <- as.vector(Hband.aT[trial,])
			# mat / rep(vec, each = nrow(mat))
			#cMband.pT[[p]] <- mat/rep(vec, each=nrow(mat))
			# mat / matrix(vec, nrow(mat), ncol(mat), byrow = TRUE)
			cMband.pT[[p]] <- mat/matrix(vec, nrow=nrow(mat), ncol=ncol(mat), byrow=T)
			# 3 different scenarios with 0 divisions: 1) 0/# = 0, 2) 0/0 = NaN, and 3) #/0 = Inf
			# for 1), do nothing, for 2), turn to 0, and for 3), turn to NaN
			# first turn NaN to 0
			cMband.pT[[p]][is.nan(cMband.pT[[p]])] <- 0
			# then turn Inf to NaN
			cMband.pT[[p]][is.infinite(cMband.pT[[p]])] <- NaN
			colnames(cMband.pT[[p]]) <- bandnames
		} ## end loop over compounds
		## form nested list with list of compound band vs time matrices
		## within the trial list
		MpH.pT[[trial]] <- cMband.pT
		Trials <- c(Trials,trial)
		trial <- trial + 1
	} ## end of loop over measure files or trials
	## average over all trials and append Time
	avgMpH.aT <- list()
	sdMpH.aT <- list()
	for (n in names(MpH.pT[[1]])) {
		tmpA <- lapply(MpH.pT, function(x) { x[[n]] })
		datavg <- tmpA[[1]]
		datstddev <- tmpA[[1]]
		for (cn in bandnames) {
			tmpB <- sapply(tmpA, function(x) { x[,cn] })
			datavg[,cn] <- apply(tmpB, 1, function(x) { mean(x, na.rm=T) })
			datstddev[,cn] <- apply(tmpB, 1, function(x) { sd(x, na.rm=T) })
		}
		datavg <- cbind(Time, datavg)
		colnames(datavg)[1] <- "Time"
		datstddev <- cbind(Time, datstddev)
		colnames(datstddev)[1] <- "Time"
		avgMpH.aT[[n]] <- datavg
		sdMpH.aT[[n]] <- datstddev
	}
	return(list(avgMpH=avgMpH.aT,sdMpH=sdMpH.aT))
}
## function to create event matrix
create.E.mat <- function(mfiles,Hdat) {
	## modification of reduce-event-data-inband.r
	cat("\n")
	E.pT <- list()
	deadTime.pT <- list()
	allEdata <- data.frame()
	Eband.aT <- vector()
	Trials <- vector()
	trial <- 1
	for (f in mfiles) {
		tdat <- read.csv(f, colClasses="numeric", check.names=F)
		tdat <- unique(tdat)
		deadTime.pT[[trial]] <- tdat["Time"] # or tdat[1]
		E.pT[[trial]] <- tdat[,c("Time","Dist_from_PV","Dist_from_CV")]
		if (direction == "dPV") deaddist <- tdat[,"Dist_from_PV"]
		else deaddist <- tdat[,"Dist_from_CV"]
		maxdeaddist <- max(deaddist)
		if (dMin > maxdeaddist) {
			cat("\n")
			print(paste("dMin",dMin,"is > max dead distance",maxdeaddist,"for direction",direction,"in file",f))
			#next
		}
		nEband.pT <- vector()
		bandnames <- vector()
		bandstart <- dMin
		bandstop <- bandstart + bandw
		bandx <- 1
		while (bandstop <= dMax) {
			bandnames <- c(bandnames,paste(direction,"∈[",bandstart,",",bandstop,")",sep=""))
			#bandat <- E.pT[[trial]][bandstart <= deaddist & deaddist < bandstop,]
			slicendx <- which(bandstart <= deaddist & deaddist < bandstop)
			## slice the data according to direction and band
			#slicendx <- which(dMin <= deaddist & deaddist < dMax)
			if (length(slicendx) == 0) {
				cat(paste("No ",measure," events ",direction,"∈[",bandstart,",",bandstop,") in ",f, "\n",sep=""))
				nEband.pT <- c(nEband.pT,0)
			} else {
				sliced <- E.pT[[trial]][slicendx,]
				nEband.pT <- c(nEband.pT,nrow(sliced))
			}
			bandstart <- bandstart + bandw
			bandstop <- bandstop + bandw
			bandx <- bandx + 1
		}
		Eband.aT <- rbind(Eband.aT,nEband.pT)
		## stack trials vertically
		allEdata <- rbind(allEdata,E.pT[[trial]])
		Trials <- c(Trials,trial)
		trial <- trial + 1
	} ## end of loop over measure files or trials
	#Eband.aT <- matrix(Eband.aT, nrow=nrow(Eband.aT), ncol=ncol(Eband.aT))
	rownames(Eband.aT) <- c()
	colnames(Eband.aT) <- bandnames
	#Eband.aT <- cbind(Trials,Eband.aT)
    #colnames(Eband.aT) <- c("Trial",bandnames)
    #sumE.aT <- apply(Eband.aT, 2, function(x) { sum(x, na.rm=T) })
	EbaT.totals <- colSums(Eband.aT)
	#avgE.aT <- apply(Eband.aT, 2, function(x) { mean(x, na.rm=T) })
	#sapply(data, mean, na.rm = T)
	#numdata<-data[sapply(data, is.numeric)]  
	#sapply(numdata, mean, na.rm = T)
	#colSds in matrixStats package
	avgE.aT <- colMeans(Eband.aT)
	#sqrt(diag(cov(data_matrix)))
	#colSd <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=sd, na.rm=na.rm)
	sdE.aT <- apply(Eband.aT, 2, function(x) { sd(x, na.rm=T) })
    #EbaT.totals <- matrix(EbaT.totals,nrow=1, ncol=length(bandnames))
    #colnames(EbaT.totals) <- bandnames
    if (any(is.na(Hdat))) {
		EpH.aT <- NA
		EpHaT.totals <- NA
    } else {
		## calculate the number of events per hepatocyte
		Hband.aT <- Hdat[[1]]
		HbaT.totals <- Hdat[[2]]
		EpH.aT <- Eband.aT/Hband.aT
		EpHaT.totals <- EbaT.totals/HbaT.totals
    }
    ## form time vs band matrix with each element containing the number
    ## of events or number of events per hepatocyte
    allEdata <- allEdata[order(allEdata[,"Time"]),]
    #indexName <- colnames(allEdata)[1]
    Eband.pTime <- vector()
	EbpH.pTime <- vector()
	if (direction == "dPV") cndd <- "Dist_from_PV" 
	else cndd <- "Dist_from_CV" 
	 ## loop through unique times
	Times <- vector()
	for (t in unique(allEdata[,"Time"])) {
		dd.pTime <- allEdata[t == allEdata[,"Time"],cndd]
		nEb.pTime <- vector()
		bandnames <- vector()
		bandstart <- dMin
		bandstop <- bandstart + bandw
		bandx <- 1
		while (bandstop <= dMax) {
			bandnames <- c(bandnames,paste(direction,"∈[",bandstart,",",bandstop,")",sep=""))
			slicendx <- which(bandstart <= dd.pTime & dd.pTime < bandstop)
			nEb.pTime <- c(nEb.pTime,length(slicendx))
			bandstart <- bandstart + bandw
			bandstop <- bandstop + bandw
			bandx <- bandx + 1
		}
		Eband.pTime <- rbind(Eband.pTime,nEb.pTime)
		Times <- c(Times,t)
	} ## end of loop over time
	if (any(is.na(Hdat))) {
		EbpH.pTime <- NA
	} else {
		## calculate the number of events per hepatocyte
		HbaT.totals <- Hdat[[2]]
		mat <- as.matrix(Eband.pTime)
		vec <- as.vector(HbaT.totals)
		# mat / matrix(vec, nrow(mat), ncol(mat), byrow = TRUE)
		EbpH.pTime <- mat/matrix(vec, nrow=nrow(mat), ncol=ncol(mat), byrow=T)
		# 3 different scenarios with 0 divisions: 1) 0/# = 0, 2) 0/0 = NaN, and 3) #/0 = Inf
		# for 1), do nothing, for 2), turn to 0, and for 3), turn to NaN
		# first turn NaN to 0
		EbpH.pTime[is.nan(EbpH.pTime)] <- 0
		# then turn Inf to NaN
		EbpH.pTime[is.infinite(EbpH.pTime)] <- NaN
	}
	rownames(Eband.pTime) <- c()
	rownames(EbpH.pTime) <- c()
	Eband.pTime <- cbind(Times,Eband.pTime)
	colnames(Eband.pTime) <- c("Time",bandnames)
	EbpH.pTime <- cbind(Times,EbpH.pTime)
	colnames(EbpH.pTime) <- c("Time",bandnames)
    return(list(EbandpTime=Eband.pTime,EbpHpTime=EbpH.pTime))
}

calc.derivative <- function(dat.time, dat) {
	dat.dxdt <- diff(as.matrix(dat[,2:ncol(dat)]))
	dat.tmp <- as.data.frame(dat.dxdt)
	dat.tmp <- cbind(dat.time[-length(dat.time)],dat.dxdt)
	colnames(dat.tmp) <- colnames(dat)
	return(dat.tmp)
}

direction <- argv[1]
dMin <- as.numeric(argv[2])
dMax <- as.numeric(argv[3])
bandw <- as.numeric(argv[4])
exps <- argv[-(1:4)] ## all remaining args
nxs <- length(exps)
#drange <- paste(direction,"∈[",dMin,",",dMax,")",sep="")

filebases <- c("celladj", "mobileObject", "entries", "exits", "rejects", "traps","necrotic","nectrig","stressed")
#measurements <- filebases
#measurements <- c("entries","celladj","mobileObject")
measurements <- c("enzymes","rxnprod","celladj","mobileObject")

## loop over experiments
for (x in exps) {
	x <- paste(dirname(x), basename(x), sep="/") # removes trailing slash from the argument
	print(paste("Working on", x))
	if (!file.exists(x)) {
		print(paste(x,"doesn't exist."))
		next
	}
	pb <- txtProgressBar(min=0,max=length(measurements),style=3)  ## progress bar
	setTxtProgressBar(pb,0); ## progress bar
	mNdx <- 1
	## loop over measures in measurements
	for (measure in measurements) {
		## read data files
		if (measure %in% c("nectrig","necrotic","stressed")) {
			mpattern <- paste(measure,"-[0-9]+.csv",sep="")
			## event measure is true
			event <- T
		} else {
			mpattern <- paste(measure,"-",direction,"-[0-9]+.csv.gz",sep="")
			event <- F
		}
		mfiles <- list.files(path=x,pattern=mpattern,recursive=T,full.names=T)
		if (length(mfiles) <= 0) {
			print(paste("skipping ", measure, "-", direction, sep=""))
			next
		}
		hpattern <- paste("hcount-",direction,"-[0-9]+.csv",sep="")
		hfiles <- list.files(path=x,pattern=hpattern,recursive=T,full.names=T)
		Hdat <- create.H.mat(hfiles)
		if (event) {
			E.aT <- create.E.mat(mfiles,Hdat)
			## write files
			cmpnames <- names(E.aT)
			for (cmp in cmpnames) {
				write.csv(x=E.aT[[cmp]], file=paste(x, "_",measure,"-",cmp,"-",direction,"-xyz.csv", sep=""), row.names=FALSE)
			}
		} else {
			avgM.aT <- create.M.mat(mfiles)
			if (any(is.na(Hdat))) {
				statMpH.aT <- NA
			} else {
				statMpH.aT <- create.MpH.mat(avgM.aT,mfiles,Hdat)
			}
			## write files
			cmpnames <- names(avgM.aT)
			for (cmp in cmpnames) {
				write.csv(x=avgM.aT[[cmp]], file=paste(x, "_",measure,"-",cmp,"-",direction,"-xyz.csv", sep=""), row.names=FALSE)
			}
			for (n in names(statMpH.aT)) {
				cmpnames <- names(statMpH.aT[[n]])
				for (cn in cmpnames) {
					write.csv(x=statMpH.aT[[n]][[cn]], file=paste(x, "_",n,"-",measure,"-",cn,"-",direction,"-xyz.csv", sep=""), row.names=FALSE)
				}
			}
		}
		setTxtProgressBar(pb,mNdx); ## progress bar
		mNdx <- mNdx + 1
	} ## end loop over measure in measurements
	close(pb) ## end progress bar
} ## end loop over experiments
quit()
