#! /usr/bin/Rscript

###
## Calculates the sums, averages, and standard deviations of Hepatocytes
## for each dCV and dPV bands over all MC trials using the raw data
## of the number of Hepatoctyes at each distance in the hcount files.
##
## Time-stamp: <2019-07-18 16:20:39 gepr>
###

argv <- commandArgs(T)

require(stats) # for statistics
tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

usage <- function() {
  print("Usage: Hcounts-avgsd.r <dCV|dPV> <dMin> <dMax> <exp directories>")
  print("  directories should contain files hcount-d[CP]V-[0-9a-z]+.csv")
  quit()
}
if (length(argv) < 4) usage()

direction <- argv[1]
dMin <- as.numeric(argv[2])
dMax <- as.numeric(argv[3])
exps <- argv[-(1:3)] ## all remaining args

Hfb <- "hcount"

for (x in exps) {
  x <- paste(dirname(x), basename(x), sep="/") # removes trailing slash from the argument
  print(paste("Working on", x))

  if (!file.exists(x)) {
    print(paste(x,"doesn't exist."))
    next
  }

  id <- 1
  xpt <- list()

  hfile <- paste(Hfb,"-",direction,"-[0-9a-z]+.csv",sep="")
  allhfiles <- list.files(path=x, pattern=hfile, recursive=T,  full.names=T)
  nHfiles <- length(allhfiles)

  if (nHfiles <= 0) {
    print(paste("no hepatocyte counts for exp ", x))
    next
  }

  nMCtrials <- nHfiles
  ## progress bar
  print(paste("distance = ",direction))
  pb <- txtProgressBar(min=0,max=nMCtrials,style=3)
  setTxtProgressBar(pb,0)

  ## loop over each MC trial
  for (i in 1:nMCtrials) {
    ## read hepatocyte count data from file
    ## first row (distances) becomes column names when read
    Hdat <- read.csv(allhfiles[i], check.names=F)
    if (i == 1) {
      Hcounts <- Hdat
    } else {
      ## combine H counts by distance from each trial
      Hcounts <- pad1stColumns(Hcounts, Hdat)
      Hpad <- pad1stColumns(Hdat, Hcounts)
      Hcounts <- rbind(Hcounts, Hpad)
    }
	
	## check if dMin is > maximum distance for this direction
    ## If true, then skip this MC trial.
    distances <- as.numeric(colnames(Hdat))
    maxdist <- max(as.numeric(distances))
    if (dMin > maxdist) {
		cat("\n")
        print(paste("dMin,",dMin,", is > max distance,",maxdist,", for direction,",direction))
        next
    }
    ## select Hcounts within band
    Hband <- Hdat[,dMin<=as.numeric(colnames(Hdat)) & as.numeric(colnames(Hdat))<dMax]
    
    ## if there is only 1 column for this band, convert it to a 1 column DF
    if (ncol(as.matrix(Hband)) == 1) {
		Hband <- as.data.frame(t(Hband))
		colnames(Hband) <- dMin
	}
	
    if (i == 1) {
      Hbandcts <- Hband
    } else {
      ## combine H counts by distance from each trial
      Hbandcts <- pad1stColumns(Hbandcts, Hband)
      Hpad <- pad1stColumns(Hband, Hbandcts)
      Hbandcts <- rbind(Hbandcts,Hpad)
    }
    setTxtProgressBar(pb,i); ## progress bar
  } ## loop over MC trials
  setTxtProgressBar(pb,i); ## progress bar
  close(pb) ## progress bar
  
  Hsums_pT <- rowSums(Hcounts)
  Hsums_pD <- colSums(Hcounts)
  sumH_pT <- sum(Hsums_pT)
  sumH_pD <- sum(Hsums_pD)
  ## these should be the same number
  if (sumH_pT != sumH_pD) {
    print(paste("sum of vHPCs over trials", sumH_pT, "!= sum of vHPCs over", direction, sumH_pD))
    q("no")
  }
  ## below is the same as rowMeans and colMeans
  ##avgMpHbandMC <- apply(MpHbandMC, c(1,2), function(x) { mean(x, na.rm=TRUE) })
  avgHcts_pT <- rowMeans(Hcounts)
  avgHcts_pD <- colMeans(Hcounts)
  avgHsums_pT <- mean(Hsums_pT)
  sdHcts_pT <- apply(Hcounts, 1, function(x) { sd(x, na.rm=TRUE) })
  sdHcts_pD <- apply(Hcounts, 2, function(x) { sd(x, na.rm=TRUE) })
  sdHsums_pT <- sd(Hsums_pT)
  ## add sums, avgs, and sds to columns, stats over MC trials
  alldistances <- colnames(Hcounts)
  stats <- c("total", "avg", "stddev")
  trialnames <- paste("Trial", 1:nMCtrials)
  ##filler <- matrix(0, nrow=3, ncol=3)
  allfiller <- matrix(c(sumH_pT, 0, 0, 0, avgHsums_pT, 0, 0, 0, sdHsums_pT), nrow=3, ncol=3)
  Hcounts <- cbind(Hcounts, Hsums_pT, avgHcts_pT, sdHcts_pT)
  pD <- rbind(as.vector(Hsums_pD), as.vector(avgHcts_pD), as.vector(sdHcts_pD))
  pDplusfill <- cbind(pD, allfiller)
  Hcounts <- rbind(as.matrix(Hcounts), pDplusfill)
  colnames(Hcounts) <- c(alldistances, stats)
  rownames(Hcounts) <- c(trialnames, stats)
  ## write hepatocyte stats to files, pD = per Distance, pT = per MC trial
  fileprefix <- paste(x, "_", "Hcounts-all-", direction, sep="")
  write.csv(Hcounts, file=paste(fileprefix, ".csv", sep=""))
  
  
  if (!exists("Hbandcts")) {
	cat("No Hs in",direction,"∈[",dMin,",",dMax,")\n")
  } else {
	Hbandsm_pT <- rowSums(Hbandcts)
	Hbandsm_pD <- colSums(Hbandcts)
	sumHband_pT <- sum(Hbandsm_pT)
	sumHband_pD <- sum(Hbandsm_pD)
	## these should be the same number
	if (sumHband_pT != sumHband_pD) {
		print(paste("sum of vHPCs over trials", sumHband_pT, "!= sum of vHPCs over", direction, sumHband_pD))
		q("no")
	}
	avgHband_pT <- rowMeans(Hbandcts)
	avgHband_pD <- colMeans(Hbandcts)
	avgHbandsm_pT <- mean(Hbandsm_pT)
	sdHband_pT <- apply(Hbandcts, 1, function(x) { sd(x, na.rm=TRUE) })
	sdHband_pD <- apply(Hbandcts, 2, function(x) { sd(x, na.rm=TRUE) })
	sdHbandsm_pT <- sd(Hbandsm_pT)
	
	banddistances <- colnames(Hbandcts)
	bandfiller <- matrix(c(sumHband_pT, 0, 0, 0, avgHbandsm_pT, 0, 0, 0, sdHbandsm_pT), nrow=3, ncol=3)
	Hbandcts <- cbind(Hbandcts, Hbandsm_pT, avgHband_pT, sdHband_pT)
	pD <- rbind(as.vector(Hbandsm_pD), as.vector(avgHband_pD), as.vector(sdHband_pD))
	pDplusfill <- cbind(pD, bandfiller)
	Hbandcts <- rbind(as.matrix(Hbandcts), pDplusfill)
	colnames(Hbandcts) <- c(banddistances,stats)
	rownames(Hbandcts) <- c(trialnames,stats)
	fileprefix <- paste(x, "_", "Hcounts-band-", direction, sep="")
    write.csv(Hbandcts, file=paste(fileprefix, "∈[", dMin, ",", dMax, ").csv", sep=""))
  }
} # loop over experiments

quit()
