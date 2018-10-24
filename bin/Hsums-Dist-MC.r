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
			#cat("\n")
			Hdat <- read.csv(allhfiles[i], check.names=F)
            #print(Hdat)
            #cat("\n")
            #Hdat <- read.csv(allhfiles[i], colClasses="numeric")
            #print(Hdat)
            
			# two sets of data for each measure and H count
			# one = totals over whole lobule, other = totals within distance band
			# select measure data within band
			HbandMC <- Hdat[,dMin<=as.numeric(colnames(Hdat)) & as.numeric(colnames(Hdat))<dMax]

			# calculate the total number of Hepatocytes and within the band
			if (exists("HsumbandMC")) {
			   HsumbandMC <- cbind(HsumbandMC,rowSums(HbandMC))
			} else {
				if (length(HbandMC) > 1) {
					HsumbandMC <- rowSums(HbandMC)
				} else {
					HsumbandMC <- HbandMC
				}
			}
			if (exists("HsumallMC")) {
			   HsumallMC <- cbind(HsumallMC,rowSums(Hdat))
			} else {
			   HsumallMC <- rowSums(Hdat)
			}

			## set data from 1st file as totals, then sum with subsequent data
			if (i == 1) {
			   Htotals <- Hdat
			} else {
			   # combine H counts by distance from each trial
			   Htotals <- pad1stColumns(Htotals, Hdat)
			   Hpad <- pad1stColumns(Hdat, Htotals)
			   Htotals <- Htotals[,order(names(Htotals))]
			   Hpad <- Hpad[,order(names(Hpad))]
			   Htotals <- Htotals + Hpad
			}
			setTxtProgressBar(pb,i); ## progress bar
		} ## loop over MC trials
		setTxtProgressBar(pb,i); ## progress bar
		close(pb) ## progress bar

		# two sets of data for H count
		# one = totals over whole lobule, other = totals within distance band
		Hinband <- Htotals[,dMin<=as.numeric(colnames(Htotals)) & as.numeric(colnames(Htotals))<dMax]
		HsumallMC <- as.data.frame(as.matrix(HsumallMC))
		colnames(HsumallMC) <- paste(1:nMCtrials)
		HsumbandMC <- as.data.frame(as.matrix(HsumbandMC))
		colnames(HsumbandMC) <- paste(1:nMCtrials)

		# write hepatocyte files, pD = per Distance, pT = per MC trial
		Hband_pT <- HsumbandMC
		Htotal_pT <- HsumallMC
		Hband_pD <- Hinband
		Htotal_pD <- Htotals

		fileprefix <- paste(x,"_","Hsums-pMC-",dir,sep="")
		write.csv(cbind(Htotal_pT, total = rowSums(Htotal_pT)), file=paste(fileprefix,".csv",sep=""), row.names=F)

		fileprefix <- paste(x,"_","Hsums-pDist-",dir,sep="")
		write.csv(cbind(Htotal_pD, total = rowSums(Htotal_pD)), file=paste(fileprefix,".csv",sep=""), row.names=F)

		fileprefix <- paste(x,"_","Hsums-pMC-",dir,sep="")
		write.csv(cbind(Hband_pT, total = rowSums(Hband_pT)), file=paste(fileprefix,"∈[",dMin,",",dMax,").csv",sep=""), row.names=F)

		fileprefix <- paste(x,"_","Hsums-pDist-",dir,sep="")
		write.csv(cbind(Hband_pD, total = rowSums(Hband_pD)), file=paste(fileprefix,"∈[",dMin,",",dMax,").csv",sep=""), row.names=F)

		# remove variables for next measure distance pair
			 
		remove(HsumbandMC)
		remove(HsumallMC)
		remove(Hinband)
		remove(Htotals)
    } # loop over distances
} # loop over experiments
