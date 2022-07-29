#! /usr/bin/Rscript

###
## Cumulative statistics for the whole Lobule.
##
## Time-stamp: <2019-05-13 09:47:06 gepr>
###
argv <- commandArgs(TRUE)

usage <- function() {
  print("Usage: reduce-event-data.r <exp directories>")
  print("  directories should contain files named necrotic-[0-9]+.csv")
  quit()
}

library(stats) # for median and sd
library(stringr) # for prefix manipulation

if (length(argv) <= 0) usage()

prefixes <- c("necrotic","nectrig","stressed")

for (prefix in prefixes) {
  cat(paste("Working on",prefix,"\n"))

  for (expDir in argv) {
    expDir <- paste(dirname(expDir), basename(expDir), sep="/") # removes trailing slash from the argument

    cat(paste("\t", "experiment ", expDir,"\n",sep=""))
    if (!file.exists(expDir)) {
      print(paste(expDir," doesn't exist."))
      next
    }
    files <- list.files(path = expDir, pattern=paste(prefix,"-[0-9]+.csv",sep=""), recursive = TRUE)
    if (length(files) <= 0) {
      cat(paste("WARNING: No",prefix,"files in",expDir,"\n"))
      next
    }

    eventData <- data.frame()
    # row bind (concat vertically) the data
    trial <- 1
    Trials <- 1:length(files)
    totalEvent <- vector()
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
      totalEvent[trial] <- nrow(tmp)
      if (totalEvent[trial] == 0) next
      event_dPV <- tmp[,"Dist_from_PV"]
      max_dPV[trial] <- max(event_dPV)
      min_dPV[trial] <- min(event_dPV)
      avg_dPV[trial] <- mean(event_dPV)
      median_dPV[trial] <- median(event_dPV)
      stddev_dPV[trial] <- sd(event_dPV)
      event_dCV <- tmp[,"Dist_from_CV"]
      max_dCV[trial] <- max(event_dCV)
      min_dCV[trial] <- min(event_dCV)
      avg_dCV[trial] <- mean(event_dCV)
      median_dCV[trial] <- median(event_dCV)
      stddev_dCV[trial] <- sd(event_dCV)
      eventData <- rbind(eventData,tmp)
      trial <- trial + 1
    }

    if (all(totalEvent) == 0) {
      print(paste("no",str_to_title(prefix),"events over all MC trials"))
      q()
    }

    ##
    # accumulate event statistics
    ##

    # whole experiment
    ExpStatEvent <- data.frame()
    total_exp <- sum(totalEvent)
    avg_exp <- mean(totalEvent)
    stddev_exp <- sd(totalEvent)
    exp_event_dPV <- eventData[,"Dist_from_PV"]
    min_dPV_exp <- min(exp_event_dPV)
    max_dPV_exp <- max(exp_event_dPV)
    avg_dPV_exp <- mean(exp_event_dPV)
    median_dPV_exp <- median(exp_event_dPV)
    stddev_dPV_exp <- sd(exp_event_dPV)
    exp_event_dCV <- eventData[,"Dist_from_CV"]
    min_dCV_exp <- min(exp_event_dCV)
    max_dCV_exp <- max(exp_event_dCV)
    avg_dCV_exp <- mean(exp_event_dCV)
    median_dCV_exp <- median(exp_event_dCV)
    stddev_dCV_exp <- sd(exp_event_dCV)
    
    StatEvents <- data.frame()
	trial.mean <- rep(NA,length(Trials))
    trial.stddev <- rep(NA,length(Trials))
	StatEvents <- cbind(Trials, totalEvent, trial.mean, trial.stddev,
						min_dPV, max_dPV, avg_dPV, median_dPV, stddev_dPV,
                        min_dCV, max_dCV, avg_dCV, median_dCV, stddev_dCV)
    
    whole.exp <- c("all",total_exp, avg_exp, stddev_exp, 
						min_dPV_exp, max_dPV_exp, avg_dPV_exp, median_dPV_exp, stddev_dPV_exp,
						min_dCV_exp, max_dCV_exp, avg_dCV_exp, median_dCV_exp, stddev_dCV_exp)
    
    StatEvents <- rbind(StatEvents,whole.exp)  
                        
    colnames(StatEvents) <- c("Trial","total","mean","stddev",
							"min_dPV", "max_dPV","avg_dPV", "median_dPV", "stddev_dPV",
							"min_dCV", "max_dCV", "avg_dCV", "median_dCV","stddev_dCV")
    
    write.csv(x=StatEvents, file=paste(expDir, "_events_stats_",prefix,".csv", sep=""), row.names=FALSE)

    ## sort the whole data set based on time
    eventData <- eventData[order(eventData[,1]),]

    ## loop through unique times
    indexName <- colnames(eventData)[1]
    newEvent <- data.frame()
    for (time in unique(eventData[,1])) {
      attach(eventData)
      tmpdf <- eventData[eventData[indexName] == time,]
      detach(eventData)
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
    write.csv(x=newEvent, file=paste(expDir, "_", prefix, ".csv", sep=""), row.names=FALSE)

  } ## end for (expDir in argv)

} ## for (prefix in prefixes)
