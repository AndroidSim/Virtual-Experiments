#! /usr/bin/Rscript

###
## Average body samples from several monte carlo trials and write to a file.
##
## Time-stamp: <2018-08-23 09:35:57 gepr>
###

argv <- commandArgs(TRUE)

require(stats) # for statistics

if (length(argv) <= 0) {
    print(paste("Usage: body-avg.r <exp directories>"))
    print("  directories should contain files like body-[0-9]+.csv, ")
    print("  with a common name and common headers for columns. e.g.:")
    print("  ./2014-04-07-1645/body-00[0-9][0-9].csv")
    quit()
}


## for each experiment
for (d in argv) {

  print (paste("d = ",d))

  timeisset <- FALSE
  run <- vector()
  ## get all the enzymes files in that experiment
  files <- list.files(path = d, pattern = "(body|medium)-[0-9]+.csv", recursive = TRUE)
  compname <- unlist(strsplit(files[1],"-"))[1]

  ## for each trial file (run)
  whole <- vector("list") # list mode vector, list of matrices
  trial <- 1
  for (f in files) {

    odata <- read.csv(file = paste(d, f, sep="/"), colClasses = "numeric")
    dims <- dim(odata)
    ##print(paste("str odata = ",str(odata)))
    ## time column
    if (timeisset == FALSE) {
      time <- odata[1]
      timeisset <- TRUE
    }

    ## the rest of the columns
    whole[[trial]] <- odata[2:(length(odata))]
    trial <- trial + 1

  }

  if (length(whole) == 1) {
    dataavg <- cbind(time, whole)
  } else {
    dataavg <- whole[[1]]
    datastddev <- whole[[1]]
    for (i in 1:ncol(whole[[1]])) {
      tmp <- sapply(whole, function(x) { x[[i]] })
      dataavg[,i] <- apply(tmp, 1, function(x) { mean(x, na.rm=TRUE) })
      datastddev[,i] <- apply(tmp, 1, function(x) { sd(x, na.rm=TRUE) })
    }
    dataavg <- cbind(time, dataavg)
    datastddev <- cbind(time, datastddev)
    colnames(dataavg)[1] <- "Time"
    colnames(datastddev)[1] <- "Time"
    ## write the data to a file
    write.csv(x=datastddev, file=paste(d, "_", compname,"-stddev.csv", sep=""), row.names=FALSE)
  }
  write.csv(x=dataavg, file=paste(d, "_", compname,"-avg.csv", sep=""), row.names=FALSE)

}

