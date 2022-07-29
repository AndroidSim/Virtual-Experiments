#! /usr/bin/Rscript

###
## Calculates the measure (e.g. Solute amounts) per Hepatocytes for dCV and dPV bands
## using the raw data in the mobileObject, celladj, and MITs files, as will as
## the number of Hepatoctyes at each distance in the hcount files.
##
## Also, calculates statistics over the MC trials and time, then prints them to files.
##
###

argv <- commandArgs(T)

require(stats) # for statistics

usage <- function() {
  print("Usage: dataperH-inband.r <dPV|dCV> dMin dMax <exp directories>")
  print("  directories should contain files like mobileObject-d[CP]V-[0-9]+.csv.gz, ")
  print (" celladj-d[CP]V-[0-9]+.csv.gz, and hcount-d[CP]V-[0-9]+.csv")
  quit()
}

if (length(argv) < 4) usage()

tryCatch(source("~/R/misc.r"), warning = function(w) { print("Could not source the library."); print(w$message); quit(status=-1); })

direction <- argv[1]
dMin <- as.numeric(argv[2])
dMax <- as.numeric(argv[3])
exps <- argv[-(1:3)] ## all remaining args

measurefbs <- c("enzymes","celladj", "mobileObject", "boundMObject", "entries", "exits", "rejects", "traps")
Hfb <- "hcount"

for (x in exps) {
  x <- paste(dirname(x), basename(x), sep="/") # removes trailing slash from the argument
  print(paste("Working on", x))

  if (!file.exists(x)) {
    print(paste(x,"doesn't exist."))
    next
  }

  ## for all measures
  for (measure in measurefbs) {

    mfile <- paste(measure,"-",direction,"-[0-9]+.csv.gz",sep="")
    hfile <- paste(Hfb,"-",direction,"-[0-9]+.csv",sep="")
    allmfiles <- list.files(path=x, pattern=mfile, recursive=T,  full.names=T)
    allhfiles <- list.files(path=x, pattern=hfile, recursive=T,  full.names=T)
    nMfiles <- length(allmfiles)
    nHfiles <- length(allhfiles)
    if (nMfiles <= 0) {
      print(paste("skipping ", measure, "-", direction, sep=""))
      next
    }
    if (nHfiles <= 0) {
      print(paste("no hepatocyte counts for distance ", direction," in exp ", x))
      next
    }
    if (nMfiles != nHfiles) {
      print(paste("number of measure files != number of Hepatocyte count files"))
      next
    }
    nMCtrials <- nMfiles
    ## progress bar
    print(paste("measure = ",measure,", direction = ",direction))
    pb <- txtProgressBar(min=0,max=nMCtrials,style=3)
    setTxtProgressBar(pb,0)

    Trialsused <- vector()
    ## loop over each MC trial
    for (i in 1:nMCtrials) {
      ## read measure data from file
      Mdat <- read.csv(allmfiles[i], check.names=F)
      Time <- Mdat[,1]

      ## check if dMin is > maximum distance for this direction
      ## If true, then skip this MC trial.
      distances <- unlist(strsplit(colnames(Mdat),":"))
      distances <- unique(distances[c(F,T)])
      Maxdist <- max(as.numeric(distances))
      if (dMin > Maxdist) {
        cat("\n")
        print(paste("dMin, ",dMin,", is > max distance, ",Maxdist,", for direction, ",direction))
        next
      }
      Trialsused <- c(Trialsused,i)

      ## get unique measure groups and distance between band range
      ## column name format = distance:measure
      groups <- vector()
      colmatches <- vector()
      for (cn in colnames(Mdat)[2:ncol(Mdat)]) {
        splitted <- unlist(strsplit(cn,":"))
        d <- as.numeric(splitted[1])
        g <- splitted[2]
        if (dMin <= d && d < dMax) colmatches <- c(colmatches,cn)
        groups <- c(groups,g)
      }
      groups <- unique(groups)
      groups <- groups[order(groups)]

      ## read hepatocyte count data from file
      ## first row (distances) becomes column names when read
      Hdat <- read.csv(allhfiles[i], check.names=F)

      ## two sets of data for each measure and H count
      ## one = totals over whole lobule, other = totals within distance band
      ## select measure data within band
      MbandMC <- Mdat[,colmatches]
      HbandMC <- Hdat[,dMin<=as.numeric(colnames(Hdat)) & as.numeric(colnames(Hdat))<dMax]

      ## sum data for all distances (whole lobule) and within distance band
      for (g in groups) {
        ## note the "$" in the grep() call delimits the search string
        ## (without it, a grep for "G" finds both "G" and "GSH_..."
        gsumbandMC <- MbandMC[,grep(paste(":",g,"$",sep=""),colnames(MbandMC))]
        ## if distances within band are > number of groups, need to sum
        if (ncol(MbandMC) > length(groups)) {
          gsumbandMC <- rowSums(gsumbandMC)
        }
        ## add column of group sum above to existing group sum data structure
        if (exists("allgbandMC")) {
          allgbandMC <- cbind(allgbandMC,gsumbandMC)
        } else {
          allgbandMC <- gsumbandMC
        }
      }
      if (exists("MallgbandMC")) {
        colnames(allgbandMC) <- paste(i,groups,sep=":")
        MallgbandMC <- cbind(MallgbandMC, allgbandMC)
      } else {
        MallgbandMC <- allgbandMC
        colnames(MallgbandMC) <- paste(i,groups,sep=":")
      }
      remove(allgbandMC)

      ## calculate the number of Hepatocytes within the band
      if (exists("HsumbandMC")) {
        if (length(HbandMC) > 1) {
          HsumbandMC <- cbind(HsumbandMC,rowSums(HbandMC))
        } else {
          HsumbandMC <- cbind(HsumbandMC,HbandMC)
        }
      } else {
        if (length(HbandMC) > 1) {
          HsumbandMC <- rowSums(HbandMC)
        } else {
          HsumbandMC <- HbandMC
        }
      }
      setTxtProgressBar(pb,i); ## progress bar
    } ## loop over MC trials
    setTxtProgressBar(pb,i); ## progress bar
    close(pb) ## progress bar

    ## got the following error for the below if block:
    ## Error in is.empty(Trialsused) : could not find function "is.empty"
    ## Execution halted
    ## if (is.empty(Trialsused)) {
    ##	print("Trialsused is empty")
    ##}
    ## This would be better if it worked

    if (!exists("HsumbandMC")) {
      cat("No Hs in",direction,"∈[",dMin,",",dMax,")\n")
      q()
    }

    HsumbandMC <- as.data.frame(as.matrix(HsumbandMC))
    colnames(HsumbandMC) <- Trialsused
    ##colnames(HsumbandMC) <- paste(1:nMCtrials)

    ## divide each group by the number of hepatocytes for that trial
    ## column name format = MCtrial:measure
    gpHbandMC <- MallgbandMC
    for (cn in colnames(MallgbandMC)) {
      splitted <- unlist(strsplit(cn,":"))
      t <- splitted[1]
      g <- splitted[2]
      gpHbandMC[,cn] <- MallgbandMC[,cn]/HsumbandMC[,t]
    }

### future: for loop over MC trials could be more efficient? ###

    ## calculate statistics for measure/hepatocyte over MC trials
    ## average and standard deviation
### use rowMeans function for better efficiency ? ###
    for (g in groups) {
      ## note the "$" in the grep() call delimits the search string (without it, a grep for "G" finds both "G" and "GSH_..."
      MpHbandMC <- gpHbandMC[,grep(paste(":",g,"$",sep=""),colnames(gpHbandMC))]
      ## avg and sd amount/HPC for each MC trial within the band
      if (is.vector(MpHbandMC)) {
        if (exists("avgMpHbandMC")) {
          avgMpHbandMC <- cbind(avgMpHbandMC,MpHbandMC)
        } else {
          avgMpHbandMC <- MpHbandMC
        }
        if (exists("sdMpHbandMC")) {
          sdMpHbandMC <- cbind(sdMpHbandMC,numeric(length(MpHbandMC)))
        } else {
          sdMpHbandMC <- numeric(length(MpHbandMC))
        }
      }
      if (is.matrix(MpHbandMC)) {
        if (exists("avgMpHbandMC")) {
          avgMpHbandMC <- cbind(avgMpHbandMC,apply(MpHbandMC, 1, function(x) { mean(x, na.rm=TRUE) }))
        } else {
          avgMpHbandMC <- apply(MpHbandMC, 1, function(x) { mean(x, na.rm=TRUE) })
        }
        if (exists("sdMpHbandMC")) {
          sdMpHbandMC <- cbind(sdMpHbandMC,apply(MpHbandMC, 1, function(x) { sd(x, na.rm=TRUE) }))
        } else {
          sdMpHbandMC <- apply(MpHbandMC, 1, function(x) { sd(x, na.rm=TRUE) })
        }
      } ## end if (is.matrix(MpHbandMC)) {
    } ## loop over groups

    ## label statistics and add back Time vector
    avgMpHbandMC <- cbind(Time,avgMpHbandMC)
    colnames(avgMpHbandMC) <- c("Time",groups)
    sdMpHbandMC <- cbind(Time,sdMpHbandMC)
    colnames(sdMpHbandMC) <- c("Time",groups)

    ## write statistics files
    ## average
    fileprefix <- paste(x,"_",measure,"-","avg-pHPC-pMC","-",direction,sep="")
    write.csv(avgMpHbandMC, file=paste(fileprefix,"∈[",dMin,",",dMax,").csv",sep=""), row.names=F)

    ## standard deviation
    fileprefix <- paste(x,"_",measure,"-","sd-pHPC-pMC","-",direction,sep="")
    write.csv(sdMpHbandMC, file=paste(fileprefix,"∈[",dMin,",",dMax,").csv",sep=""), row.names=F)

    ## remove variables for next measure distance pair
    remove(MallgbandMC)
    remove(HsumbandMC)
    remove(avgMpHbandMC)
    remove(sdMpHbandMC)

  } # loop over measures
} # loop over experiments
