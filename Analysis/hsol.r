#! /usr/bin/Rscript

##
# This script sums metabolite data across all monte carlo trials, per
# cell per time
#
# -gepr 2013-05-15
# -aks 2017-07-25

argv <- commandArgs(TRUE)

if (length(argv) <= 0) {
    print("Usage: hsol.r <exp directories>")
    print("  directories should contain files named hsolute_zone_1_2-[0-9]+.csv")
    quit()
}
# for the color space max and min
minmean <-  9e10
maxmean <- -9e10

avgPerZone <- function(path, fileNameRoot, extractZone) {

    timeisset <- FALSE
    ei <- vector()

    # get all the zone 1&2 reaction product files in that experiment
    files <- list.files(path = path,
                        pattern = paste(fileNameRoot,"_",extractZone,"-[0-9]+.csv",sep=""),
                        recursive = TRUE)

    # for each node, for each time, sum over all files

    # for each file
    zoneData <- vector()
    count <- 1
    for (f in files) {

        rxndata <- read.csv(file = paste(path, f, sep="/"), colClasses = "numeric")
        dims <- dim(rxndata)
        cat("Read ", dims, " from file ", f, "\n")

        # time column
        if (timeisset == FALSE) {
            rxn <- rxndata[1]
            timeisset <- TRUE
        }

        rxnnames <- list()
        for (c in 2:length(rxndata))
          rxnnames[c-1] <- unlist(strsplit(unlist(strsplit(colnames(rxndata)[c],'X'))[2],'\\.'))[6]
        umn <- unique(rxnnames)
        numnames <- length(umn)
        
        # the rest of the columns into respective zones
        for (colNdx in seq(2,length(rxndata),numnames)) {
          if (unlist(strsplit(unlist(strsplit(colnames(rxndata)[colNdx],'X'))[2],'\\.'))[1] == extractZone) {
            for (mn in 1:numnames) {
              if (count == 1) {
                zoneData <- cbind(zoneData, rxndata[,(colNdx+mn-1)])
              } else zoneData[,mn] <- zoneData[,mn] + rxndata[,(colNdx+mn-1)]
            }
            count <- count + 1
          } 
        } # end colNdx loop
        
    }
    time <- rxn
    zData <- zoneData
    hepCount <- count
    
    zoneData <- zoneData/count

    print("Binding time, zoneData")
    
    # table of reaction product per cell for each time
    rxn <- cbind(rxn, zoneData)
    # set the column names
    zoneNames <- c()
    for (hn in umn) {
      zoneNames <- c(zoneNames, paste("Z",extractZone, hn))
    }
    colnames(rxn) <- c("Time", zoneNames)
	colnames(zData) <- umn

    print("writing data to the file")
    
    # write the rxn sums data to a file
    write.csv(x=rxn,
    file=paste(path, "_", fileNameRoot, "-",extractZone,".csv", sep=""),
    row.names=FALSE)
    
    output <- list("time" = time, "zData" = zData, "hepCount" = hepCount)
    return(output)
}

# for each experiment
for (expDir in argv) {
    z0 <- avgPerZone(expDir, "hsolute_zone", 0)
    z1 <- avgPerZone(expDir, "hsolute_zone", 1)
    z2 <- avgPerZone(expDir, "hsolute_zone", 2)
    
    time <- z0[[1]]
    
    z0data <- z0[[2]]
    z0hepcount <- z0[[3]]
    print(paste("z0hepcount = ",z0hepcount))
    
    z1data <- z1[[2]]
    z1hepcount <- z1[[3]]
    print(paste("z1hepcount = ",z1hepcount))
    
    z2data <- z2[[2]]
    z2hepcount <- z2[[3]]
    print(paste("z2hepcount = ",z2hepcount))
    
    total_data <- z0data + z1data + z2data
    total_count <- z0hepcount + z1hepcount + z2hepcount
    print(paste("total hepcount = ",total_count))
    
    totalperH <- total_data/total_count
    
    # create data structure to output to file
    data2file <- cbind("Time" = time, totalperH)
    # write data to file
    write.csv(x=data2file,file=paste(expDir, "_hsolute_total.csv", sep=""),row.names=FALSE)
}

q()
