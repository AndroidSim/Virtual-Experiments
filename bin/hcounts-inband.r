#! /usr/bin/Rscript
###
## Slices out the hcount-????.csv datat within the given band, and prints the
## average number of Hepatocytes in that band.
##
## Time-stamp: <2018-06-04 14:06:05 gepr>
###
argv <- commandArgs(TRUE)

if (length(argv) < 4 || !(argv[1] == "dCV" || argv[1] == "dPV")) {
  print("Usage: hcounts-inband.r <dCV or dPV> min max exp1 exp2 ...")
  q()
}

direction <- argv[1]
dmin <- as.numeric(argv[2])
dmax <- as.numeric(argv[3])
exps <- argv[-(1:3)]

for (exp in exps) {
  files <- list.files(path = exp, pattern=paste("hcount-",direction,"-[0-9]+.csv",sep=""), recursive=T)
  total <- 0
  for (f in files) {
    fn <- paste(exp,"/",f,sep="")
    d <- read.csv(fn, colClasses="numeric")
    if (dmin > ncol(d)) next
    thismax <- ifelse(dmax > ncol(d), ncol(d), dmax)
    row <- d[1,(dmin+1):thismax]
    rowtot <- sum(row)
    print(cbind(row[1,], rowtot))
    total <- total+rowtot
  }
  avg <- total/length(files)
  print(paste(exp,": avg =",avg))
}

