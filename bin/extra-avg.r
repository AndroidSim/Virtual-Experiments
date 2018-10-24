#! /usr/bin/Rscript
####! /bin/bash setR

###
## Average extra samples from several monte carlo trials and write to a file.
##
## -gepr 2017-04-18
###

#dev.off()

argv <- commandArgs(TRUE)

if (length(argv) <= 0) {
    print(paste("Usage: extra-avg.r <exp directories>"))
    print("  Requires there be multiple files below the PWD")
    print("  with a common name and common headers for columns. e.g.:")
    print("  ./2014-04-07-1645/extra-00[0-9][0-9].csv")
    quit()
}


# for each experiment
for (d in argv) {

    print (paste("d = ",d))

    timeisset <- FALSE
    run <- vector()
    # get all the enzymes files in that experiment
    files <- list.files(path = d, pattern = "extra-[0-9]+.csv", recursive = TRUE)
    compname <- unlist(strsplit(files[1],"-"))[1]

    # for each trial file (run)
    whole <- vector("list") # list mode vector, list of matrices
    trial <- 1
    for (f in files) {

        odata <- read.csv(file = paste(d, f, sep="/"), colClasses = "numeric")
        dims <- dim(odata)
        #print(paste("str odata = ",str(odata)))
        # time column
        if (timeisset == FALSE) {
            run <- odata[1]
            timeisset <- TRUE
        }

        # the rest of the columns
        whole[[trial]] <- odata[2:(length(odata))]
        trial <- trial + 1

    }
    
    if (length(whole) == 1) {
      expmean <- whole[[1]]
    } else {
      # get the averages for each node, for each time, over all trials
      trialsum <- whole[[1]]
      for (i in 2:length(whole)-1) {
        
        #print(paste("adding in", whole[[i+1]]))
        try(trialsum+whole[[i+1]])
        if (geterrmessage() != "") q()
        trialsum <- trialsum + whole[[i+1]]
      }
      expmean <- trialsum/length(whole)
    }

    # table of outflow per node per time
    run <- cbind(run, expmean)
    # set the column names
    colnames(run) <- colnames(odata)[1:length(colnames(odata))]

    # write the data to a file
    write.csv(x=run, file=paste(d, "_", compname,".csv", sep=""), row.names=FALSE)
}

#q()
