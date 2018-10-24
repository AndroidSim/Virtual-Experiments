#! /usr/bin/Rscript

##
# This script does several things to enzyme count data:
# 1) average data across all monte carlo trials, per SS node
# 2) sums over all SS nodes in each zone
# 3) takes finite diff of those sums
#
#
# -gepr 2011-06-10
#
# bkp 2012-07-10 Now works for non-constant dt\
# gepr 2012-09-10 copied from ISL & changed to work with ISLJ
# gepr 2012-09-11 clean up code and normalize by # of nodes per zone
#

argv <- commandArgs(TRUE)

if (length(argv) <= 0) {
    print("Usage: ei.r <exp directories>")
    print("  directories should contain files named enzymes_[0-9]+.csv")
    quit()
}
#dev.off()
# for the color space max and min
minmean <-  9e10
maxmean <- -9e10

# for each experiment
for (d in argv) {

    timeisset <- FALSE
    ei <- vector()
    # get all the enzymes files in that experiment
    files <- list.files(path = d, pattern = "enzymes-[0-9]+.csv", recursive = TRUE)
    
    # for each node, for each time, sum over all files
    
    # for each enzyme file (run)
    whole <- vector("list") # list mode vector, list of matrices
    run <- 1
    for (f in files) {

        edata <- read.csv(file = paste(d, f, sep="/"), colClasses = "numeric")
        dims <- dim(edata)
        
        # time column
        if (timeisset == FALSE) {
            ei <- edata[1]
            
            timeisset <- TRUE
        }
        
        # the rest of the columns
        whole[[run]] <- edata[2:(length(edata))]
        run <- run + 1
        
    }
    
    print(paste("length of whole = "),length(whole))
    if (length(whole) == 1) {
      expmean <- whole[[1]]
    } else {
      # get the averages for each node, for each time, over all runs
      runsum <- whole[[1]]
      for (i in 2:length(whole)-1) { 
        #print(paste("adding in", whole[[i+1]]))
        runsum <- runsum + whole[[i+1]] 
      }
      expmean <- runsum/length(whole)
    }
    
    # table of enzymes per node per time
    ei <- cbind(ei, expmean)
    # set the column names
    colnames(ei) <- colnames(edata)

    # average between 2 observations for each time step
    started <- FALSE
    ei.avg <- matrix(ncol=length(ei[1,]))
    avg <- vector()
    for (i in seq(from=1, to=length(ei[,1]), by=2)) {
        one <- ei[row(ei) == i]
        if (length(ei[row(ei) == i+1]) <= 0) break
        two <- ei[row(ei) == i+1]
        avg <- (one + two)/2
        if (started) { ei.avg <- rbind(ei.avg, avg) }
        else {
            ei.avg[1,] <- avg
            started <- TRUE
        }
    }
    colnames(ei.avg) <- colnames(ei)
    
    # write the data to a file
    # write.csv(x=ei.avg,
    # file=paste(d, "mean_enzyme_count_per_node.csv", sep="_"),
    # row.names=FALSE)
    
    # sum mean enzyme amounts for all nodes in each zone
    
    # specify which nodes in which zones
    nz <- vector("list") 
    
    # store indices of string "ZoneX"
    zone0 <- grep("Z_0", colnames(ei.avg))
    zone1 <- grep("Z_1", colnames(ei.avg))
    zone2 <- grep("Z_2", colnames(ei.avg))
    
    # sum the averages for each node into the zone and divide by # of nodes in that zone
    zone.sums <- cbind(ei.avg[,1], # time
    rowSums(ei.avg[,zone0])/length(zone0), # zone I
    rowSums(ei.avg[,zone1])/length(zone1), # zone II
    rowSums(ei.avg[,zone2])/length(zone2)) # zone III
    colnames(zone.sums) <- c("Time", "Zone 1", "Zone 2", "Zone 3")
    
    # write the sums to a file
    write.csv(x=zone.sums,
    file=paste(d, "mean_enzyme_count_per_zone.csv", sep="_"),
    row.names=FALSE)
    
    # finite diffs of the means per zone
    dt <- diff(zone.sums[,1]) # works for non-constant dt
    zone.dxdt <- cbind(zone.sums[-1,1],
    diff(zone.sums[,2])/dt, # zone I
    diff(zone.sums[,3])/dt, # zone II
    diff(zone.sums[,4])/dt) # zone III
    colnames(zone.dxdt) <- colnames(zone.sums)
    
    write.csv(x=zone.dxdt,
    file=paste(d, "mean_enzyme_count_per_zone-dxdt.csv", sep="_"),
    row.names=FALSE)
}

q()
