#! /usr/bin/Rscript

###
## Average body & extra samples from several monte carlo trials and write to a file.
##
## Time-stamp: <2019-05-13 09:43:58 gepr>
###

argv <- commandArgs(TRUE)

if (length(argv) <= 0) {
  print(paste("Usage: coarse-compartment.r <exp directories>"))
  print("  directories should contain files like body-[0-9]+.csv and extra[0-9]+.csv,")
  print("  with a common name and common headers for columns. e.g.:")
  print("  ./2014-04-07-1645/body-00[0-9][0-9].csv")
  quit()
}

## for each file type
## the Intro measure was added in the islj.r1169-icm branch and if those changes are merged in, the "intro" type can be included
##for (ft in c("body", "extra", "intro")) {
for (ft in c("body", "extra")) {

  ## for each experiment
  for (d in argv) {
    d <- paste(dirname(d), basename(d), sep="/") # removes trailing slash from the argument
    print (paste("Reducing", ft, "for d =",d))

    timeisset <- F
    run <- vector()

    files <- list.files(path = d, pattern = paste(ft, "-[0-9]+.csv", sep=""), recursive = T)
    if (length(files) <= 0) {
      print(paste("Warning!", ft, "files don't exist in experiment", d))
      next()
    }
    ##compname <- unlist(strsplit(files[1],"-"))[1]
    compname <- ft

    ## for each trial file (run)
    whole <- vector("list") # list mode vector, list of matrices
    trial <- 1
    for (f in files) {

      odata <- read.csv(file = paste(d, f, sep="/"), colClasses = "numeric")
      dims <- dim(odata)
      ##print(paste("str odata = ",str(odata)))
      ## time column
      if (timeisset == F) {
        time <- odata[1]
        timeisset <- T
      }

      ## the rest of the columns
      whole[[trial]] <- odata[2:(length(odata))]
      trial <- trial + 1

    } ## end for (f in files) {

    ## if there's only 1 trial, use that
    if (length(whole) == 1) {
      dataavg <- cbind(time, whole)
    } else {
      dataavg <- whole[[1]]
      datastddev <- whole[[1]]
      for (i in 1:ncol(whole[[1]])) {
        tmp <- sapply(whole, function(x) { x[[i]] })
        dataavg[,i] <- apply(tmp, 1, function(x) { mean(x, na.rm=T) })
        datastddev[,i] <- apply(tmp, 1, function(x) { sd(x, na.rm=T) })
      }
      dataavg <- cbind(time, dataavg)
      datastddev <- cbind(time, datastddev)
      colnames(dataavg)[1] <- "Time"
      colnames(datastddev)[1] <- "Time"
      ## write the data to a file
      write.csv(x=datastddev, file=paste(d, "_", compname,"-stddev.csv", sep=""), row.names=F)
    }
    write.csv(x=dataavg, file=paste(d, "_", compname,"-avg.csv", sep=""), row.names=F)

  } ## end for (d in argv) {

} ## end for (ft in c("body", "extra")) {
