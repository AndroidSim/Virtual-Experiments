#! /usr/bin/Rscript

###
## Time-stamp: <2018-05-21 16:14:40 gepr>
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

library(stats) # for median and sd
library(stringr) # for prefix manipulation

if (length(argv) <= 0) {
    print("Usage: necrotic.r <exp directories>")
    print("  directories should contain files named necrotic-[0-9]+.csv")
    quit()
}

dMin <- as.numeric(argv[1])
dMax <- as.numeric(argv[2])
exps <- argv[-(1:2)]
prefixes <- c("necrotic","nectrig","stressed")

for (prefix in prefixes) {
  cat(paste("Working on",prefix,"\n"))

  for (expDir in exps) {
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
      tmp <- read.csv(file = paste(expDir,f,sep="/"), colClasses = "numeric")
      tmp <- unique(tmp)
      totalEvent[trial] <- nrow(tmp)
      if (totalEvent[trial] == 0) next
      event_dPV <- tmp[,"Dist_from_PV"]
      event_dPV <- event_dPV[dMin <= event_dPV & event_dPV < dMax]

      max_dPV[trial] <- max(event_dPV)
      min_dPV[trial] <- min(event_dPV)
      avg_dPV[trial] <- mean(event_dPV)
      median_dPV[trial] <- median(event_dPV)
      stddev_dPV[trial] <- sd(event_dPV)
      event_dCV <- tmp[,"Dist_from_CV"]
      event_dCV <- event_dCV[dMin <= event_dCV & event_dCV < dMax]

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

    # all trials
    TrialStatEvent <- data.frame()
    TrialStatEvent <- cbind(Trials, totalEvent, min_dPV, max_dPV, min_dCV,
                              max_dCV, avg_dPV, median_dPV, stddev_dPV, avg_dCV, median_dCV, stddev_dCV)
    colnames(TrialStatEvent) <- c("Trial", paste("total_",prefix,sep=""), "min_dPV", "max_dPV",
                                    "min_dCV", "max_dCV", "avg_dPV", "median_dPV", "stddev_dPV", "avg_dCV", "median_dCV", "stddev_dCV")

    # whole experiment
    ExpStatEvent <- data.frame()
    total_exp <- sum(totalEvent)
    avg_exp <- mean(totalEvent)
    stddev_exp <- sd(totalEvent)
    exp_event_dPV <- eventData[,"Dist_from_PV"]
    exp_event_dPV <- exp_event_dPV[dMin <= exp_event_dPV & exp_event_dPV < dMax]

    min_dPV_exp <- min(exp_event_dPV)
    max_dPV_exp <- max(exp_event_dPV)
    avg_dPV_exp <- mean(exp_event_dPV)
    median_dPV_exp <- median(exp_event_dPV)
    stddev_dPV_exp <- sd(exp_event_dPV)
    exp_event_dCV <- eventData[,"Dist_from_CV"]
    exp_event_dCV <- exp_event_dCV[dMin <= exp_event_dCV & exp_event_dCV < dMax]

    min_dCV_exp <- min(exp_event_dCV)
    max_dCV_exp <- max(exp_event_dCV)
    avg_dCV_exp <- mean(exp_event_dCV)
    median_dCV_exp <- median(exp_event_dCV)
    stddev_dCV_exp <- sd(exp_event_dCV)
    ExpStatEvent <- rbind(total_exp, avg_exp, stddev_exp, min_dPV_exp, max_dPV_exp, min_dCV_exp, max_dCV_exp,
                            avg_dPV_exp, median_dPV_exp, stddev_dPV_exp, avg_dCV_exp, median_dCV_exp, stddev_dCV_exp)
    rownames(ExpStatEvent) <- c("total","mean","stddev","min_dPV", "max_dPV",
                                  "min_dCV", "max_dCV", "avg_dPV", "median_dPV", "stddev_dPV", "avg_dCV", "median_dCV","stddev_dCV")
    colnames(ExpStatEvent) <- expDir

    # write statistics to files
    write.csv(x=ExpStatEvent, file=paste(expDir, "_stat_exp_",prefix,"∈[",dMin,",",dMax,").csv", sep=""), row.names=TRUE)
    write.csv(x=TrialStatEvent, file=paste(expDir, "_stat_MCtrials_",prefix,"∈[",dMin,",",dMax,").csv", sep=""), row.names=FALSE)

    ## sort the whole data set based on time
    eventData <- eventData[order(eventData[,1]),]

    indexName <- colnames(eventData)[1]

    ## aggregate for dCV
    dCVevents <- eventData[dMin <= eventData$Dist_from_CV & eventData$Dist_from_CV < dMax,]

    ## loop through unique times
    newEvent <- data.frame()
    if (nrow(dCVevents) <= 0) {
      cat(paste("No",prefix,"dCV ∈ [",dMin,",",dMax,")\n"))
    } else {
      for (time in unique(dCVevents[,1])) {
        attach(dCVevents)
        tmpdf <- dCVevents[dCVevents[indexName] == time,]
        detach(dCVevents)
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
      write.csv(x=newEvent, file=paste(expDir, "_", prefix,"-dCV∈[",dMin,",",dMax,").csv", sep=""), row.names=FALSE)
    } # end else

    ## aggregate for dPV
    dPVevents <- eventData[dMin <= eventData$Dist_from_PV & eventData$Dist_from_PV < dMax,]
    if (nrow(dPVevents) <= 0) {
      cat(paste("No",prefix,"dPV ∈ [",dMin,",",dMax,")\n"))
    } else {
      ## loop through unique times
      newEvent <- data.frame()
      for (time in unique(dPVevents[,1])) {
        attach(dPVevents)
        tmpdf <- dPVevents[dPVevents[indexName] == time,]
        detach(dPVevents)
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
      write.csv(x=newEvent, file=paste(expDir, "_", prefix,"-dPV∈[",dMin,",",dMax,").csv", sep=""), row.names=FALSE)
    } # end else

  } ## end for (expDir in argv)

} ## for (prefix in prefixes)
