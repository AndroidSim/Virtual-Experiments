#! /usr/bin/Rscript

###
## Calculates the sum of Hepatocytes for each dCV and dPV bands over all MC trials
## using the raw data of the number of Hepatoctyes at each distance in the hcount files.
##
###

argv <- commandArgs(T)

require(stats) # for statistics

usage <- function() {
  print("Usage: Measure_per_HPC-inband.r dMin dMax <exp directories>")
  print("  directories should contain files like hsolute-d[CP]V-[0-9]+.csv.gz, ")
  print (" celladj-d[CP]V-[0-9]+.csv.gz, and hcount-d[CP]V-[0-9]+.csv")
  quit()
}

if (length(argv) < 3) usage()

source("~/R/misc.r")
'%!in%' <- function(x,y)!('%in%'(x,y))

dMin <- as.numeric(argv[1])
dMax <- as.numeric(argv[2])
exps <- argv[-(1:2)] ## all remaining args
Hfb <- "hcount"

for (x in exps) {
    print(paste("Working on", x))

    if (!file.exists(x)) {
		print(paste(x,"doesn't exist."))
        next
    }

    ## for both directions
    for (dir in c("dPV","dCV")) {
		hfile <- paste(Hfb,"-",dir,"-[0-9]+.csv",sep="")
		allhfiles <- list.files(path=x, pattern=hfile, recursive=T,  full.names=T)
		nHfiles <- length(allhfiles)
			 
		if (nHfiles <= 0) {
			print(paste("no hepatocyte counts for exp ", x))
			next
		}
			 
		nMCtrials <- nHfiles
		# progress bar
		print(paste("distance = ",dir))
		pb <- txtProgressBar(min=0,max=nMCtrials,style=3)
		setTxtProgressBar(pb,0)

		# loop over each MC trial
		for (i in 1:nMCtrials) {
			# read hepatocyte count data from file
			# first row (distances) becomes column names when read
			Hdat <- read.csv(allhfiles[i], check.names=F)
            if (i == 1) {
				Hcounts <- Hdat
            } else {
				# combine H counts by distance from each trial
				Hcounts <- pad1stColumns(Hcounts, Hdat)
				Hpad <- pad1stColumns(Hdat, Hcounts)
				Hcounts <- rbind(Hcounts,Hpad)
            } 
            
			# select Hcounts within band
			Hband <- Hdat[,dMin<=as.numeric(colnames(Hdat)) & as.numeric(colnames(Hdat))<dMax]
			if (i == 1) {
				Hbandcts <- Hband
            } else {
				# combine H counts by distance from each trial
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
		Hbandsm_pT <- rowSums(Hbandcts)
		Hbandsm_pD <- colSums(Hbandcts)
		
		# below is the same as rowMeans and colMeans
		#avgMpHbandMC <- apply(MpHbandMC, c(1,2), function(x) { mean(x, na.rm=TRUE) })
		
		avgHcts_pT <- rowMeans(Hcounts) 
		avgHcts_pD <- colMeans(Hcounts)
		avgHband_pT <- rowMeans(Hbandcts)
		avgHband_pD <- colMeans(Hbandcts)
		
		sdHcts_pT <- apply(Hcounts, 1, function(x) { sd(x, na.rm=TRUE) })
		sdHcts_pD <- apply(Hcounts, 2, function(x) { sd(x, na.rm=TRUE) })
		sdHband_pT <- apply(Hbandcts, 1, function(x) { sd(x, na.rm=TRUE) })
		sdHband_pD <- apply(Hbandcts, 2, function(x) { sd(x, na.rm=TRUE) })
		
		# add sums, avgs, and sds to columns, stats over MC trials
		alldistances <- colnames(Hcounts)
		banddistances <- colnames(Hbandcts)
		stats <- c("total","avg","stddev")
		trialnames <- paste("Trial",1:nMCtrials)
		filler <- matrix(0, nrow=3, ncol=3)
		
		Hcounts <- cbind(Hcounts, Hsums_pT, avgHcts_pT, sdHcts_pT)
		pD <- rbind(as.vector(Hsums_pD), as.vector(avgHcts_pD), as.vector(sdHcts_pD))
		pDplusfill <- cbind(pD,filler)
		Hcounts <- rbind(as.matrix(Hcounts), pDplusfill)
		colnames(Hcounts) <- c(alldistances,stats)
		rownames(Hcounts) <- c(trialnames,stats)
		
		Hbandcts <- cbind(Hbandcts, Hbandsm_pT, avgHband_pT, sdHband_pT)
		pD <- rbind(as.vector(Hbandsm_pD), as.vector(avgHband_pD), as.vector(sdHband_pD))
		pDplusfill <- cbind(pD,filler)
		Hbandcts <- rbind(as.matrix(Hbandcts), pDplusfill)
		colnames(Hbandcts) <- c(banddistances,stats)
		rownames(Hbandcts) <- c(trialnames,stats)
		
		# write hepatocyte files, pD = per Distance, pT = per MC trial
		fileprefix <- paste(x,"_","Hcounts-all-",dir,sep="")
		write.csv(Hcounts, file=paste(fileprefix,".csv",sep=""))

		fileprefix <- paste(x,"_","Hcounts-band-",dir,sep="")
		write.csv(Hbandcts, file=paste(fileprefix,"∈[",dMin,",",dMax,").csv",sep=""))
    } # loop over distances
} # loop over experiments
