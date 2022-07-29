#! /usr/bin/Rscript

###
## Time-stamp: <2020-01-06 12:01:14 gepr>
##
## Reduce the event data ([necrotic|nectrig|stressed]-????.csv the  files
## in the following way:
##
## sum the number of cells that necrosed at each time
## average the distances of all cells that necrosed at that time
## output:
##    time, #necrosed, avg_dpv, median_dpv, sigma_dpv, avg_dcv, median_dcv, sigma_dcv
##
###
argv <- commandArgs(TRUE)

usage <- function() {
  print("Usage: reduce-event-data-inband.r <dPV|dCV> <dMin> <dMax> <exp directories>")
  print("  directories should contain files named necrotic-[0-9]+.csv")
  quit()
}

library(stats) # for median and sd
library(stringr) # for prefix manipulation

if (length(argv) < 4) usage();

direction <- argv[1]
dMin <- as.numeric(argv[2])
dMax <- as.numeric(argv[3])
exps <- argv[-(1:3)]
prefixes <- c("necrotic","nectrig","stressed")

for (prefix in prefixes) {
  skipThisExp <- F
  for (expDir in exps) {
    expDir <- paste(dirname(expDir), basename(expDir), sep="/") # removes trailing slash from the argument
    cat(paste(prefix, " for experiment ", expDir,"\n",sep=""))
    if (!file.exists(expDir)) {
      print(paste(expDir," doesn't exist."))
      next
    }
    files <- list.files(path = expDir, pattern=paste(prefix,"-[0-9]+.csv",sep=""), recursive = TRUE)
    if (length(files) <= 0) {
      cat(paste("WARNING: No",prefix,"files in",expDir,"\n"))
      next ## prefix
    }

    eventData <- data.frame()
    # row bind (concat vertically) the data
    trial <- 1
    Trials <- 1:length(files)
    nEvent.band <- vector()
    event_dPV <- vector()
    event_dCV <- vector()
    max_dPV <- vector()
    min_dPV <- vector()
    max_dCV <- vector()
    min_dCV <- vector()
    avg_dPV <- vector()
    avg_dCV <- vector()
    median_dPV <- vector()
    median_dCV <- vector()
    stddev_dPV <- vector()
    stddev_dCV <- vector()
    for (f in files) {
      # individual trials
      tmp <- read.csv(file = paste(expDir,f,sep="/"), colClasses = "numeric")
      tmp <- unique(tmp)

      ## slice the data according to direction and band
      attach(tmp)
      if (direction == "dPV") slicer <- Dist_from_PV
      else slicer <- Dist_from_CV
      detach(tmp)
      slicendx <- which(dMin <= slicer & slicer < dMax)
      if (length(slicendx) == 0) {
        cat(paste("No ", prefix, " events ", direction, "∈[",dMin,",",dMax,") in ", f, "\n",sep=""))
        nEvent.band[trial] <- NA
        max_dPV[trial] <- NA
        min_dPV[trial] <- NA
        avg_dPV[trial] <- NA
        median_dPV[trial] <- NA
        stddev_dPV[trial] <- NA

        max_dCV[trial] <- NA
        min_dCV[trial] <- NA
        avg_dCV[trial] <- NA
        median_dCV[trial] <- NA
        stddev_dCV[trial] <- NA
      } else {
        sliced <- tmp[slicendx,]
        nEvent.band[trial] <- nrow(sliced)
        event_dPV <- sliced[,"Dist_from_PV"]
        max_dPV[trial] <- max(event_dPV)
        min_dPV[trial] <- min(event_dPV)
        avg_dPV[trial] <- mean(event_dPV)
        median_dPV[trial] <- median(event_dPV)
        stddev_dPV[trial] <- sd(event_dPV)

        event_dCV <- sliced[,"Dist_from_CV"]
        max_dCV[trial] <- max(event_dCV)
        min_dCV[trial] <- min(event_dCV)
        avg_dCV[trial] <- mean(event_dCV)
        median_dCV[trial] <- median(event_dCV)
        stddev_dCV[trial] <- sd(event_dCV)
      }

      ## stack trials vertically
      eventData <- rbind(eventData,tmp)
      trial <- trial + 1
    }

    ##
    # accumulate event statistics
    ##

    # whole experiment
    ExpStatEvent <- data.frame()
    total_exp <- sum(nEvent.band)
    avg_exp <- mean(nEvent.band)
    stddev_exp <- sd(nEvent.band)

    ## slice the data according to direction and band
    attach(eventData)
    if (direction == "dPV") slicer <- Dist_from_PV
    else slicer <- Dist_from_CV
    detach(eventData)
    slicendx <- which(dMin <= slicer & slicer < dMax)
    if (length(slicendx) == 0) {
      cat(paste("No ", prefix, " events ", direction, "∈[",dMin,",",dMax,") in entire experiment.\n",sep=""))
      skipThisExp <- T
      min_dPV_exp <- NA
      max_dPV_exp <- NA
      avg_dPV_exp <- NA
      median_dPV_exp <- NA
      stddev_dPV_exp <- NA

      min_dCV_exp <- NA
      max_dCV_exp <- NA
      avg_dCV_exp <- NA
      median_dCV_exp <- NA
      stddev_dCV_exp <- NA
    } else {
      sliced <- eventData[slicendx,]
      exp_event_dPV <- sliced[,"Dist_from_PV"]

      min_dPV_exp <- min(exp_event_dPV)
      max_dPV_exp <- max(exp_event_dPV)
      avg_dPV_exp <- mean(exp_event_dPV)
      median_dPV_exp <- median(exp_event_dPV)
      stddev_dPV_exp <- sd(exp_event_dPV)

      exp_event_dCV <- sliced[,"Dist_from_CV"]
      min_dCV_exp <- min(exp_event_dCV)
      max_dCV_exp <- max(exp_event_dCV)
      avg_dCV_exp <- mean(exp_event_dCV)
      median_dCV_exp <- median(exp_event_dCV)
      stddev_dCV_exp <- sd(exp_event_dCV)
    }

    if (skipThisExp) next()

    StatEvents <- data.frame()
    trial.mean <- rep(NA,length(Trials))
    trial.stddev <- rep(NA,length(Trials))
    StatEvents <- cbind(Trials, nEvent.band, trial.mean, trial.stddev,
                        min_dPV, max_dPV, avg_dPV, median_dPV, stddev_dPV,
                        min_dCV, max_dCV, avg_dCV, median_dCV, stddev_dCV)

    whole.exp <- c("all",total_exp, avg_exp, stddev_exp,
                   min_dPV_exp, max_dPV_exp, avg_dPV_exp, median_dPV_exp, stddev_dPV_exp,
                   min_dCV_exp, max_dCV_exp, avg_dCV_exp, median_dCV_exp, stddev_dCV_exp)

    StatEvents <- rbind(StatEvents,whole.exp)

    colnames(StatEvents) <- c("Trial","total","mean","stddev",
                              "min_dPV", "max_dPV","avg_dPV", "median_dPV", "stddev_dPV",
                              "min_dCV", "max_dCV", "avg_dCV", "median_dCV","stddev_dCV")

    write.csv(x=StatEvents, file=paste(expDir, "_events_stats_",prefix, "-", direction, "∈[",dMin,",",dMax,").csv", sep=""), row.names=FALSE)

    ## sort the whole data set based on time
    eventData <- eventData[order(eventData[,1]),]

    indexName <- colnames(eventData)[1]

    ## aggregate for command-line given direction
    ## slice the data according to direction and band
    attach(eventData)
    if (direction == "dPV") slicer <- Dist_from_PV
    else slicer <- Dist_from_CV
    detach(eventData)
    slicendx <- which(dMin <= slicer & slicer < dMax)
    if (length(slicendx) == 0) {
      cat(paste("No",prefix,"dPV∈[",dMin,",",dMax,")\n"))
    } else {
      events <- eventData[slicendx,]
      ## loop through unique times
      newEvent <- data.frame()
      for (time in unique(events[,1])) {
        attach(events)
        tmpdf <- events[events[indexName] == time,]
        detach(events)
        attach(tmpdf)
        summ <- c(time,
                  nrow(tmpdf),
                  mean(Dist_from_PV),
                  median(Dist_from_PV),
                  sd(Dist_from_PV),
                  mean(Dist_from_CV),
                  median(Dist_from_CV),
                  sd(Dist_from_CV))
        detach(tmpdf)
        ## row bind the summ data onto the bottom
        newEvent <- rbind(newEvent,summ)
      }

      ## append cumulative event and time normalized
      event.cumu <- cumsum(newEvent[[2]])
      event.norm <- event.cumu/max(event.cumu)
      time.norm <- newEvent[[1]]/max(newEvent[[1]])
      newEvent <- cbind(newEvent,event.cumu, time.norm,event.norm)

      colnames(newEvent) <- c(indexName, "Num_Cells",
                              "Mean-Dist_from_PV", "Median-Dist_from_PV", "SD-Dist_from_PV",
                              "Mean-Dist_from_CV", "Median-Dist_from_CV", "SD-Dist_from_CV",
                              paste("Cumu",str_to_title(prefix),sep="_"),
                              "Time_Norm",
                              paste("Cumu",str_to_title(prefix),"Norm",sep="_"))
      write.csv(x=newEvent, file=paste(expDir, "_", prefix,"-",direction,"∈[",dMin,",",dMax,").csv", sep=""), row.names=FALSE)
    } # end if (length(slicendx) == 0) {

  } ## end for (expDir in argv)

} ## for (prefix in prefixes)
