#! /usr/bin/Rscript
#! /bin/bash setR
argv <- commandArgs(TRUE)
#dev.off()

###
## Sum the columns with the 1st element of the column name is within the band
## specified by [dMin,dMax).  Column elements are separated with ":".
##
## gepr 2017-03-27
###


usage <- function() {
    print("Usage: eg-inband.r dMin dMax exp_file.csv")
    print("  exp_file.csv should look like:")
    print("    Time, 0:Phase1, 0:Phase2, ...")
    print("  Should work for either dPV or dCV data produced by the ")
    print("  reduction scripts.")
    quit()
}

if (length(argv) < 3) usage()
dMin <- as.numeric(argv[1])
dMax <- as.numeric(argv[2])
fn <- argv[3]

source("~/R/misc.r")

dat <- read.csv(fn,check.names=F)
sums <- sumBandByLastTag(dat,c(dMin,dMax))

outfilebase <- substr(fn,0,regexpr('_',fn)-1)
write.csv(sums,paste(outfilebase,"_enzymes∈[",dMin,",",dMax,").csv",sep=""),row.names=F)

#q()
