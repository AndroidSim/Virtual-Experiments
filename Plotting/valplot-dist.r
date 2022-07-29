#! /bin/bash setR

##
# This script does averages ISLJ main output data, across all monte
# carlo trials.
#
# gepr 2013-03-25 copied and modified from islj/bin/ei.r
#

if (length(argv) <= 0) {
    print("Usage: valplot-dist.r <exp directories>")
    print("  directories should contain files named run-[0-9]+.csv")
    quit()
}
dev.off()
# for the color space max and min
minmean <-  9e10
maxmean <- -9e10

# for each experiment
for (d in argv) {

    timeisset <- FALSE
    out <- vector()
    # get all the enzymes files in that experiment
    files <- list.files(path = d, pattern = "run-[0-9]+.csv", recursive = TRUE)
    
    # for each node, for each time, sum over all files
    
    # for each enzyme file (run)
    whole <- vector("list") # list mode vector, list of matrices
    run <- 1
    for (f in files) {

        odata <- read.csv(file = paste(d, f, sep="/"))
        dims <- dim(odata)
        
        # time column
        if (timeisset == FALSE) {
            out <- odata[1]
            timeisset <- TRUE
        }
        
        # the rest of the columns
        whole[[run]] <- odata[2:(length(odata))]
        run <- run + 1
        
    }
    # get the averages for each column, for each time, over all runs
    runsum <- whole[[1]]
    for (i in 2:length(whole)-1) { 
        #print(paste("adding in", whole[[i+1]]))
        runsum <- runsum + whole[[i+1]] 
    }
    expmean <- runsum/length(whole)
    
    # table of output values per column per time
    out <- cbind(out, expmean)
    # set the column names
    colnames(out) <- colnames(odata)

    # write the data to a file
    write.csv(x=out,
        file=paste(d, "mean_output.csv", sep="_"),
        row.names=FALSE)
}

q()
